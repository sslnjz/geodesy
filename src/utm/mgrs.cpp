/**********************************************************************************
*  MIT License                                                                    *
*                                                                                 *
*  Copyright (c) 2021 Binbin Song <ssln.jzs@gmail.com>                            *
*                                                                                 *
*  Geodesy tools for conversions between (historical) datums                      *
*  (c) Chris Veness 2005-2019                                                     *
*  www.movable-type.co.uk/scripts/latlong-convert-coords.html                     *
*  www.movable-type.co.uk/scripts/geodesy-library.html#latlon-ellipsoidal-datum   *
*                                                                                 *
*  Permission is hereby granted, free of charge, to any person obtaining a copy   *
*  of this software and associated documentation files (the "Software"), to deal  *
*  in the Software without restriction, including without limitation the rights   *
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      *
*  copies of the Software, and to permit persons to whom the Software is          *
*  furnished to do so, subject to the following conditions:                       *
*                                                                                 *
*  The above copyright notice and this permission notice shall be included in all *
*  copies or substantial portions of the Software.                                *
*                                                                                 *
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       *
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    *
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         *
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  *
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  *
*  SOFTWARE.                                                                      *
***********************************************************************************/

#include "mgrs.h"

#include "latlon_utm.h"
#include "strutil.h"
#include "utm_mgrs.h"

#include <cctype>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

using namespace geodesy;

namespace
{
constexpr double oneHundredKmMetres = 100e3;
constexpr double rowCycleMetres = 2000e3;

[[nodiscard]] bool containsLetter(const char* letters, char letter)
{
   return std::string(letters).find(letter) != std::string::npos;
}

[[nodiscard]] bool isDigitString(const std::string& text)
{
   if (text.empty())
      return false;

   for (char ch : text)
   {
      if (!std::isdigit(static_cast<unsigned char>(ch)))
         return false;
   }
   return true;
}

[[nodiscard]] int parseZone(const std::string& text)
{
   if (!isDigitString(text))
      throw std::invalid_argument("invalid MGRS zone: expected two digits");

   const int zone = std::stoi(text);
   if (!(1 <= zone && zone <= 60))
      throw std::invalid_argument("invalid MGRS zone: expected integer in range [1, 60]");

   return zone;
}

void validateGridZoneDesignator(int zone, char band)
{
   // Svalbard band X omits zones 32, 34, and 36; those longitudes are covered
   // by the widened neighbouring zones used by UTM and MGRS.
   if (band == 'X' && (zone == 32 || zone == 34 || zone == 36))
      throw std::invalid_argument("invalid MGRS grid zone designator");
}

[[nodiscard]] double parseOffsetMetres(const std::string& text, const char* fieldName)
{
   if (!isDigitString(text))
      throw std::invalid_argument(std::string("invalid MGRS ") + fieldName + ": expected digits");

   return std::stod(text);
}

[[nodiscard]] std::string padOffset(double offsetMetres, unsigned int digits)
{
   const double divisor = std::pow(10.0, 5.0 - digits / 2.0);
   const auto truncated = static_cast<unsigned int>(std::floor(offsetMetres / divisor));

   std::ostringstream stream;
   stream << std::setw(static_cast<int>(digits / 2)) << std::setfill('0') << truncated;
   return stream.str();
}

[[nodiscard]] std::string normaliseUnseparatedReference(const std::string& text)
{
   if (text.size() < 5)
      throw std::invalid_argument("invalid MGRS grid reference");

   const std::string eastingNorthing = text.substr(5);
   if (eastingNorthing.size() % 2 != 0)
      throw std::invalid_argument("invalid MGRS grid reference: easting/northing precision mismatch");

   return text.substr(0, 3) + " " + text.substr(3, 2) + " "
      + eastingNorthing.substr(0, eastingNorthing.size() / 2) + " "
      + eastingNorthing.substr(eastingNorthing.size() / 2);
}
}

Mgrs::Mgrs(int zone, char band, char e100k, char n100k, double easting, double northing, Datum datum)
   : m_zone(zone)
   , m_band(band)
   , m_e100k(e100k)
   , m_n100k(n100k)
   , m_easting(easting)
   , m_northing(northing)
   , m_datum(datum)
{
   if (!(1 <= zone && zone <= 60))
      throw std::invalid_argument("invalid MGRS zone: expected integer in range [1, 60]");
   if (!containsLetter(latBands, band))
      throw std::invalid_argument("invalid MGRS latitude band");
   validateGridZoneDesignator(zone, band);
   if (!containsLetter(e100kLetters[(zone - 1) % 3], e100k))
      throw std::invalid_argument("invalid MGRS 100km easting square letter");
   if (!containsLetter(n100kLetters[(zone - 1) % 2], n100k))
      throw std::invalid_argument("invalid MGRS 100km northing square letter");
   if (!std::isfinite(easting) || !(0.0 <= easting && easting < oneHundredKmMetres))
      throw std::invalid_argument("invalid MGRS easting: expected finite metres in [0, 100000)");
   if (!std::isfinite(northing) || !(0.0 <= northing && northing < oneHundredKmMetres))
      throw std::invalid_argument("invalid MGRS northing: expected finite metres in [0, 100000)");
   if (!datum.ellipsoid)
      throw std::invalid_argument("unrecognized datum");
}

UtmMgrs Mgrs::toUtm() const
{
   const Utm::Hemisphere hemisphere = m_band >= 'N' ? Utm::Hemisphere::N : Utm::Hemisphere::S;

   // The first 100 km square letter maps to full UTM easting. The +1 offset mirrors
   // the MGRS false easting pattern, whose first valid column begins at 100 km.
   const auto columnIndex = std::string(e100kLetters[(m_zone - 1) % 3]).find(m_e100k) + 1;
   const auto e100kNum = static_cast<double>(columnIndex) * oneHundredKmMetres;

   const auto rowIndex = std::string(n100kLetters[(m_zone - 1) % 2]).find(m_n100k);
   const auto n100kNum = static_cast<double>(rowIndex) * oneHundredKmMetres;

   const double latBand = static_cast<double>(std::string(latBands).find(m_band)) * 8.0 - 80.0;
   const auto bandBottomNorthing = std::floor(LatLonUtm(latBand, 0.0).toUtm().northing() / oneHundredKmMetres)
      * oneHundredKmMetres;

   // MGRS northing row letters repeat every 2,000 km. Add complete row cycles
   // until the square lies within or just above the latitude band bottom.
   double rowCycleNorthing = 0.0;
   while (rowCycleNorthing + n100kNum + m_northing < bandBottomNorthing)
      rowCycleNorthing += rowCycleMetres;

   return UtmMgrs(m_zone, hemisphere, e100kNum + m_easting,
      rowCycleNorthing + n100kNum + m_northing, m_datum, std::nullopt, std::nullopt, false);
}

std::string Mgrs::toString(unsigned int digits) const
{
   if (digits > 10 || digits % 2 != 0)
      throw std::invalid_argument("invalid MGRS precision: expected even digits in range [0, 10]");

   std::ostringstream zoneStream;
   zoneStream << std::setw(2) << std::setfill('0') << m_zone;

   const std::string prefix = zoneStream.str() + m_band + " " + m_e100k + m_n100k;
   if (digits == 0)
      return prefix;

   return prefix + " " + padOffset(m_easting, digits) + " " + padOffset(m_northing, digits);
}

Mgrs Mgrs::parse(const std::string& mgrsGridRef)
{
   const std::string stripped = strutil::strip(mgrsGridRef);
   if (stripped.empty())
      throw std::invalid_argument("invalid MGRS grid reference");

   const bool hasWhitespace = stripped.find_first_of(" \t\r\n") != std::string::npos;
   const std::string normalised = hasWhitespace ? stripped : normaliseUnseparatedReference(stripped);

   std::istringstream stream(normalised);
   std::vector<std::string> fields;
   std::string field;
   while (stream >> field)
      fields.push_back(field);

   if (fields.size() != 4)
      throw std::invalid_argument("invalid MGRS grid reference: expected GZD square easting northing");

   if (fields[0].size() != 3 || fields[1].size() != 2)
      throw std::invalid_argument("invalid MGRS grid reference");

   const int zone = parseZone(fields[0].substr(0, 2));
   const char band = fields[0][2];
   const char e100k = fields[1][0];
   const char n100k = fields[1][1];

   std::string easting = fields[2];
   std::string northing = fields[3];
   if (easting.size() != northing.size() || easting.size() > 5)
      throw std::invalid_argument("invalid MGRS grid reference: easting/northing precision mismatch");

   // Short grid references name larger squares; right-padding expresses their
   // south-west corner in metre offsets within the 100 km square.
   easting = easting + std::string(5 - easting.size(), '0');
   northing = northing + std::string(5 - northing.size(), '0');

   return Mgrs(zone, band, e100k, n100k, parseOffsetMetres(easting, "easting"),
      parseOffsetMetres(northing, "northing"));
}
