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
#include "dms.h"

#include "algorithm.h"
#include "strutil.h"

#include <cmath>
#include <iomanip>
#include <regex>
#include <sstream>
#include <vector>

namespace
{
    int defaultDecimalPlaces(geodesy::Dms::eFormat format)
    {
        switch (format)
        {
        case geodesy::Dms::D:
            return 4;
        case geodesy::Dms::DM:
            return 2;
        case geodesy::Dms::DMS:
            return 0;
        case geodesy::Dms::N:
            return 4;
        }

        return 4;
    }

    double decimalScale(int decimalPlaces)
    {
        return std::pow(10.0, decimalPlaces);
    }

    bool parseNumberToken(const std::string& token, double& value)
    {
        try
        {
            std::size_t consumed = 0;
            value = std::stod(token, &consumed);
            return consumed == token.size() && std::isfinite(value);
        }
        catch (const std::exception&)
        {
            return false;
        }
    }

    std::string formatFixed(double value, int decimalPlaces)
    {
        std::ostringstream stream;
        stream << std::fixed << std::setprecision(decimalPlaces) << value;
        return stream.str();
    }

    std::string padLeft(std::string value, std::size_t width)
    {
        if (value.size() >= width)
            return value;

        return std::string(width - value.size(), '0') + value;
    }

    std::string formatDegrees(double degrees, int decimalPlaces)
    {
        return padLeft(formatFixed(degrees, decimalPlaces), decimalPlaces == 0 ? 3 : 4 + decimalPlaces);
    }

    std::string formatMinutes(double minutes, int decimalPlaces)
    {
        return padLeft(formatFixed(minutes, decimalPlaces), decimalPlaces == 0 ? 2 : 3 + decimalPlaces);
    }

    std::string formatSeconds(double seconds, int decimalPlaces)
    {
        return padLeft(formatFixed(seconds, decimalPlaces), decimalPlaces == 0 ? 2 : 3 + decimalPlaces);
    }
}

using namespace geodesy;

std::string& Dms::mutableSeparator()
{
   // Keep the current C++ API default empty while allowing callers to opt into another separator.
   static std::string separator = "";
   return separator;
}

std::string Dms::separator()
{
   return mutableSeparator();
}

void Dms::setSeparator(const std::string& sep)
{
   mutableSeparator() = sep;
}

double Dms::parse(double degrees)
{
   return degrees;
}

double Dms::parse(const char* dms)
{
   if (dms == nullptr)
      return NAN;

   return parse(std::string(dms));
}

double Dms::parse(const std::string& dms)
{
   const std::string trimmed = strutil::strip(dms);
   if (trimmed.empty())
      return NAN;

   double decimalDegrees = NAN;
   if (std::regex_match(trimmed, std::regex(R"(^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$)")))
   {
      return parseNumberToken(trimmed, decimalDegrees) ? decimalDegrees : NAN;
   }

   // Strip the leading sign and trailing hemisphere before splitting degree, minute, and second fields.
   std::string unsignedDms = std::regex_replace(trimmed, std::regex(R"(^-)"), "");
   unsignedDms = std::regex_replace(unsignedDms, std::regex(R"(\s*[NSEW]\s*$)", std::regex_constants::icase), "");

   std::vector<std::string> dmsParts = strutil::split_regex(unsignedDms, "[^0-9.,]+");
   while (!dmsParts.empty() && dmsParts.back().empty())
      dmsParts.pop_back();
   while (!dmsParts.empty() && dmsParts.front().empty())
      dmsParts.erase(dmsParts.begin());

   if (dmsParts.empty() || dmsParts.size() > 3)
      return NAN;

   std::vector<double> numericParts;
   numericParts.reserve(dmsParts.size());
   for (const auto& part : dmsParts)
   {
      double value = NAN;
      if (!parseNumberToken(part, value))
         return NAN;
      numericParts.push_back(value);
   }

   // Interpret split fields using the reference d/m/s rules: d, d+m/60, or d+m/60+s/3600.
   switch (numericParts.size())
   {
   case 3:
      decimalDegrees = numericParts[0] + numericParts[1] / 60.0 + numericParts[2] / 3600.0;
      break;
   case 2:
      decimalDegrees = numericParts[0] + numericParts[1] / 60.0;
      break;
   case 1:
      decimalDegrees = numericParts[0];
      break;
   default:
      return NAN;
   }

   if (std::regex_search(trimmed, std::regex(R"(^-|\s*[WS]\s*$)", std::regex_constants::icase)))
      decimalDegrees = -decimalDegrees;

   return decimalDegrees;
}

std::string Dms::toDms(double deg, eFormat format, std::optional<int> dp)
{
   // give up here if we can't make a number from deg
   if (std::isnan(deg) || !std::isfinite(deg))
      return "";

   const int decimalPlaces = dp.value_or(defaultDecimalPlaces(format));
   if (decimalPlaces < 0)
      return "";

   if (format == N)
      return formatFixed(deg, decimalPlaces);

   // (unsigned result ready for appending compass dir'n)
   deg = std::fabs(deg);

   std::string dms;
   switch (format)
   {
   case D:
      {
         dms = formatDegrees(deg, decimalPlaces) + "°";
      }
      break;
   case DM:
      {
         const double scale = decimalScale(decimalPlaces);
         const double roundedMinutes = std::round(deg * 60.0 * scale) / scale;
         const auto degrees = static_cast<int>(std::floor(roundedMinutes / 60.0));
         const double minutes = roundedMinutes - degrees * 60.0;
         dms = padLeft(std::to_string(degrees), 3) + "°" + mutableSeparator()
            + formatMinutes(minutes, decimalPlaces) + "′";
      }
      break;
   case DMS:
      {
         const double scale = decimalScale(decimalPlaces);
         const double roundedSeconds = std::round(deg * 3600.0 * scale) / scale;
         const auto degrees = static_cast<int>(std::floor(roundedSeconds / 3600.0));
         const auto minutes = static_cast<int>(std::floor((roundedSeconds - degrees * 3600.0) / 60.0));
         const double seconds = roundedSeconds - degrees * 3600.0 - minutes * 60.0;
         dms = padLeft(std::to_string(degrees), 3) + "°" + mutableSeparator()
            + padLeft(std::to_string(minutes), 2) + "′" + mutableSeparator()
            + formatSeconds(seconds, decimalPlaces) + "″";
      }
      break;
   case N:
      break;
   }
   return dms;
}

std::string Dms::toLat(double deg, eFormat format, std::optional<int> dp)
{
   const std::string lat = toDms(wrap90(deg), format, dp);
   return lat.empty() ? "-" : lat.substr(1) + mutableSeparator() + (deg < 0 ? "S" : "N");
}

std::string Dms::toLon(double deg, eFormat format, std::optional<int> dp)
{
   const std::string lon = toDms(wrap180(deg), format, dp);
   return lon.empty() ? "-" : lon + mutableSeparator() + (deg < 0 ? "W" : "E");
}

std::string Dms::toBearing(double deg, eFormat format, std::optional<int> dp)
{
   const std::string brng = toDms(wrap360(deg), format, dp);
   // just in case rounding took us up to 360∼!
   if (brng.empty())
      return "-";

   return brng.rfind("360", 0) == 0 ? "0" + brng.substr(3) : brng;
}

std::string Dms::compassPoint(double bearing, int precision)
{
   if (precision < 1 || precision > 3)
      throw std::range_error("invalid precision" + std::to_string(precision));

   // note precision could be extended to 4 for quarter-winds (eg NbNW), but I think they are little used
   bearing = wrap360(bearing); // normalize to range 0..360∼
   const std::string cardinals[] = {
         "N", "NNE", "NE", "ENE",
         "E", "ESE", "SE", "SSE",
         "S", "SSW", "SW", "WSW",
         "W", "WNW", "NW", "NNW" };

   // no of compass points at req＊d precision (1=>4, 2=>8, 3=>16)
   const int n = 4 * static_cast<int>(std::exp2(precision - 1)); 
   return cardinals[static_cast<int>(std::round(bearing * n / 360)) % n * 16 / n];
}

double Dms::wrap360(double degrees)
{
   return geodesy::wrap360(degrees);
}

double Dms::wrap180(double degrees)
{
   return geodesy::wrap180(degrees);
}

double Dms::wrap90(double degrees)
{
   return geodesy::wrap90(degrees);
}
