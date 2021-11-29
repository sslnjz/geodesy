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

#include "latlon.h"
#include "dms.h"

#include <cmath>
#include <sstream>
#include <iomanip>

using namespace geodesy;

LatLon::LatLon(double lat, double lon)
   : m_lat(Dms::wrap90(lat))
   , m_lon(Dms::wrap180(lon))
{

}

double LatLon::lat() const
{
   return m_lat;
}

double LatLon::latitude() const
{
   return m_lat;
}

double LatLon::lon() const
{
   return m_lon;
}

double LatLon::lng() const
{
   return m_lon;
}

double LatLon::longitude() const
{
   return m_lon;
}

bool geodesy::operator==(const LatLon& p1, const LatLon& p2)
{
   return std::fabs(p1.m_lat - p2.m_lat) <= std::numeric_limits<double>::epsilon() &&
          std::fabs(p1.m_lon - p2.m_lon) <= std::numeric_limits<double>::epsilon();
}

bool geodesy::operator!=(const LatLon& p1, const LatLon& p2)
{
   return (!(p1 == p2));
}

void LatLon::setLat(double lat)
{
   m_lat = Dms::wrap90(lat);
}

void LatLon::setLatitude(double lat)
{
   m_lat = Dms::wrap90(lat);
}

void LatLon::setLon(double lon)
{
   m_lon = Dms::wrap180(lon);
}

void LatLon::setLng(double lon)
{
   m_lon = Dms::wrap180(lon);
}

void LatLon::setLongitude(double lon)
{
   m_lon = Dms::wrap180(lon);
}

std::string LatLon::toString(Dms::eFormat e,  std::optional<int> dp) const
{
   // note: explicitly set dp to undefined for passing through to toLat/toLon
   if(Dms::N == e && std::nullopt == dp)
      dp = 4;
   return Dms::toLat(m_lat, e, dp) + ", " + Dms::toLon(m_lon, e, dp);
}

std::string LatLon::toGeoJSON() const
{
   std::stringstream ss;
   ss << std::setprecision(9);
   ss << "{ type: \"Point\", coordinates : [";
   ss << m_lon;
   ss << ",";
   ss << m_lat;
   ss << "] }";

   return ss.str();
}
