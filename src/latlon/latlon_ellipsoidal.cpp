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
#include "latlon_ellipsoidal.h"
#include "cartesian.h"

#include <sstream>
#include <iomanip>
#include <utility>

using namespace geodesy;

LatLonEllipsoidal::LatLonEllipsoidal()
    : LatLon()
    , m_epoch(std::nullopt)
    , m_datum(std::nullopt)
    , m_referenceFrame(std::nullopt)
    , m_height(0.0)
{
}

LatLonEllipsoidal::LatLonEllipsoidal(double lat, double lon, double height,
   std::optional<Datum> datum, 
   std::optional<ReferenceFrame> reference, 
   std::optional<std::string> epoch)
   : LatLon(lat, lon)
   , m_epoch(std::move(epoch))
   , m_datum(datum)
   , m_referenceFrame(std::move(reference))
   , m_height(height)
{
}

void LatLonEllipsoidal::setHeight(double height)
{
    m_height = height;
}

void LatLonEllipsoidal::setDatum(const Datum &datum)
{
    m_datum = { datum.ellipsoid, datum.transforms };
}

double LatLonEllipsoidal::height() const
{
   return m_height;
}

Datum LatLonEllipsoidal::datum() const
{
   return *m_datum;
}

Ellipsoids LatLonEllipsoidal::ellipsoids()
{
   return g_ellipsoids;
}

Datums LatLonEllipsoidal::datums()
{
   return g_datums;
}

LatLonEllipsoidal LatLonEllipsoidal::parse(double lat, double lon, double height)
{
   const auto latlon = LatLon::parse(lat, lon);
   return LatLonEllipsoidal(latlon.lat(), latlon.lon(), height);
}

LatLonEllipsoidal LatLonEllipsoidal::parse(const std::string& dms, double height)
{
   const auto latlon = LatLon::parse(dms);
   return LatLonEllipsoidal(latlon.lat(), latlon.lon(), height);
}

LatLonEllipsoidal LatLonEllipsoidal::parse(const std::string& lat, const std::string& lon, double height)
{
   const auto latlon = LatLon::parse(lat, lon);
   return LatLonEllipsoidal(latlon.lat(), latlon.lon(), height);
}

LatLonEllipsoidal LatLonEllipsoidal::parse(const std::string& lat, const std::string& lon, std::string height)
{
   const auto latlon = LatLon::parse(lat, lon);

   double h = 0;
   try
   {
      h = std::stod(height);
   }
   catch (const std::exception& e)
   {
      throw e;
   }

   return LatLonEllipsoidal(latlon.lat(), latlon.lon(), h);
}

bool LatLonEllipsoidal::equals(const LatLonEllipsoidal& point) const
{
   return *this == point;
}

std::string LatLonEllipsoidal::toString(Dms::eFormat format, std::optional<int> dp, std::optional<int> dph) const
{
   std::string height;
   if (dph != std::nullopt) {
      height = (m_height >= 0 ? " +" : " ") + toFixed(m_height, *dph) + "m";
   }

   if (format == Dms::N) 
   { // signed numeric degrees
      if (dp == std::nullopt) dp = 4;
      return toFixed(m_lat, *dp) + ", " + toFixed(m_lon, *dp) + (dph ? height : "");
   }
   return Dms::toLat(m_lat, format, dp) + ", " + Dms::toLon(m_lon, format, dp) + (dph ? height : "");
}



Cartesian LatLonEllipsoidal::toCartesian() const
{
   // x = (ν+h)⋅cosφ⋅cosλ, y = (ν+h)⋅cosφ⋅sinλ, z = (ν⋅(1-e²)+h)⋅sinφ
   // where ν = a/√(1−e²⋅sinφ⋅sinφ), e² = (a²-b²)/a² or (better conditioned) 2⋅f-f²
   const Ellipsoid ellipsoid = m_datum
                     ? m_datum->ellipsoid
                     : m_referenceFrame ? m_referenceFrame->ellipsoid : ellipsoids().WGS84;

   const auto phi = toRadians(m_lat);
   const auto lambda = toRadians(m_lon);
   const auto h = m_height;
   const auto a = ellipsoid.a, f = ellipsoid.f;

   const auto sinphi = std::sin(phi), cosphi = std::cos(phi);
   const auto sinlambda = std::sin(lambda), coslambda = std::cos(lambda);

   const auto eSq = 2*f - f*f;                      // 1st eccentricity squared ≡ (a²-b²)/a²
   const auto nu = a / std::sqrt(1 - eSq*sinphi*sinphi); // radius of curvature in prime vertical

   const auto x = (nu+h) * cosphi * coslambda;
   const auto y = (nu+h) * cosphi * sinlambda;
   const auto z = (nu*(1-eSq)+h) * sinphi;

   return { x, y, z };
}

bool LatLonEllipsoidal::operator==(const LatLonEllipsoidal& point) const
{
   if (std::fabs(m_lat - point.m_lat) > std::numeric_limits<double>::epsilon()) return false;
   if (std::fabs(m_lon - point.m_lon) > std::numeric_limits<double>::epsilon()) return false;
   if (std::fabs(m_height - point.m_height) > std::numeric_limits<double>::epsilon()) return false;
   if (m_epoch != point.m_epoch) return false;
   if (m_datum != point.m_datum) return false;
   if (m_referenceFrame != point.m_referenceFrame) return false;

   return true;
}
