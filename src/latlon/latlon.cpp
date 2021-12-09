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

#include <cmath>
#include <sstream>
#include <iomanip>

#include "dms.h"
#include "strutil.h"


using namespace geodesy;

LatLon::LatLon() noexcept
   : m_lat(0), m_lon(0)
{
}

LatLon::LatLon(double lat, double lon)
{
   if (std::isnan(lat)) throw std::invalid_argument("invalid lat’");
   if (std::isnan(lon)) throw std::invalid_argument("invalid lon’");

   m_lat = Dms::wrap90(lat);
   m_lon = Dms::wrap180(lon);
}

LatLon::LatLon(const std::string& lat, const std::string& lon)
{
   if (std::isnan(Dms::parse(lat))) throw std::invalid_argument("invalid lat’");
   if (std::isnan(Dms::parse(lon))) throw std::invalid_argument("invalid lon’");

   m_lat = Dms::wrap90(Dms::parse(lat));
   m_lon = Dms::wrap180(Dms::parse(lon));
}

LatLon::LatLon(const LatLon& rhs) noexcept
{
   if(this != &rhs)
   {
      m_lat = rhs.m_lat;
      m_lon = rhs.m_lon;
   }
}
LatLon::LatLon(LatLon&& rhs) noexcept
{
   if (this != &rhs)
   {
      m_lat = std::exchange(rhs.m_lat, 0);
      m_lon = std::exchange(rhs.m_lon, 0);
   }
}

LatLon::~LatLon() = default;


LatLon& LatLon::operator=(const LatLon& rhs) noexcept
{
   if (this != &rhs)
   {
      m_lat = rhs.m_lat;
      m_lon = rhs.m_lon;
   }

   return *this;
}
LatLon& LatLon::operator=(LatLon&& rhs) noexcept
{
   if (this != &rhs)
   {
      m_lat = std::exchange(rhs.m_lat, 0);
      m_lon = std::exchange(rhs.m_lon, 0);
   }

   return *this;
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
   if (std::fabs(p1.m_lat - p2.m_lat) > std::numeric_limits<double>::epsilon()) return false;
   if (std::fabs(p1.m_lon - p2.m_lon) > std::numeric_limits<double>::epsilon()) return false;
   return true;
}

bool geodesy::operator!=(const LatLon& p1, const LatLon& p2)
{
   return (!(p1 == p2));
}

LatLon LatLon::parse(double lat, double lon)
{
   return { lat, lon };
}

LatLon LatLon::parse(const std::string& dms)
{
   if(const auto parts = strutil::split(dms, L','); parts.size() == 2)
   {
      return { parts[0], parts[1] };
   }

   throw std::invalid_argument("invalid point");
}

LatLon LatLon::parse(const std::string& lat, const std::string& lon)
{
   return { lat, lon };
}

std::string LatLon::toString(Dms::eFormat e,  std::optional<int> dp) const
{
   // note: explicitly set dp to undefined for passing through to toLat/toLon
   if (e == Dms::N) 
   { // signed numeric degrees
      if (dp == std::nullopt) dp = 4;
      std::stringstream ss;
      ss << std::fixed << std::setprecision(*dp);
      ss << m_lat << ", " << m_lon;
      return ss.str();
   }

   return Dms::toLat(m_lat, e, dp) + ", " + Dms::toLon(m_lon, e, dp);
}

std::string LatLon::toGeoJSON() const
{
   std::stringstream ss;
   ss << std::setprecision(9);
   ss << "{ type: \"Point\", coordinates: [";
   ss << m_lon << ", " << m_lat;
   ss << "] }";

   return ss.str();
}

bool LatLon::equals(const LatLon& point) const
{
   return *this == point;
}
