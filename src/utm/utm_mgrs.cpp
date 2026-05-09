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

#include "utm_mgrs.h"

#include "algorithm.h"
#include "latlon_utm.h"
#include "mgrs.h"

#include <cmath>
#include <stdexcept>

using namespace geodesy;

namespace
{
constexpr double oneHundredKmMetres = 100e3;
}

UtmMgrs::UtmMgrs(int zone, Hemisphere h, double easting, double northing, std::optional<Datum> datum,
   std::optional<double> convergence, std::optional<double> scale, bool verifyEN)
   : Utm(zone, h, easting, northing, datum, convergence, scale, verifyEN)
{
}

Mgrs UtmMgrs::toMgrs() const
{
   const auto latLon = toLatLon();
   if (!(-80.0 <= latLon.lat() && latLon.lat() <= 84.0))
      throw std::invalid_argument("latitude outside MGRS limits");

   // Latitude bands are 8 degrees high, with X repeated to cover 80..84 degrees north.
   const auto bandIndex = static_cast<size_t>(std::floor(latLon.lat() / 8.0 + 10.0));
   const char band = latBands[bandIndex];

   const int col = static_cast<int>(std::floor(m_easting / oneHundredKmMetres));
   if (!(1 <= col && col <= 8))
      throw std::invalid_argument("invalid UTM easting for MGRS conversion");
   const char e100k = e100kLetters[(m_zone - 1) % 3][col - 1];

   const int row = static_cast<int>(std::fmod(std::floor(m_northing / oneHundredKmMetres), 20.0));
   const char n100k = n100kLetters[(m_zone - 1) % 2][row];

   // Only the metre offsets inside the 100 km square are retained in the MGRS reference.
   const double easting = std::stod(geodesy::toFixed(std::fmod(m_easting, oneHundredKmMetres), 6));
   const double northing = std::stod(geodesy::toFixed(std::fmod(m_northing, oneHundredKmMetres), 6));

   return Mgrs(m_zone, band, e100k, n100k, easting, northing, m_datum.value());
}
