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

using namespace geodesy;

LatLonEllipsoidal::LatLonEllipsoidal()
    : m_epoch(std::nullopt)
    , m_datum(std::nullopt)
    , m_referenceFrame(std::nullopt)
    , m_lat(0.0)
    , m_lon(0.0)
    , m_height(0.0)
{
}

LatLonEllipsoidal::~LatLonEllipsoidal()
= default;

LatLonEllipsoidal::LatLonEllipsoidal(double lat, double lon, double height,
                                     std::optional<Datum> datum, std::optional<ReferenceFrame> reference, std::optional<float> epoch)
    : m_epoch(epoch)
    , m_datum(datum)
    , m_referenceFrame(reference)
    , m_lat(lat)
    , m_lon(lon)
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

double LatLonEllipsoidal::lat() const
{
   return m_lat;
}

double LatLonEllipsoidal::latitude() const
{
   return m_lat;
}

void LatLonEllipsoidal::setLat(double lat)
{
   m_lat = Dms::wrap90(lat);
}

void LatLonEllipsoidal::setLatitude(double lat)
{
   m_lat = Dms::wrap90(lat);
}

double LatLonEllipsoidal::lon() const
{
   return m_lon;
}

double LatLonEllipsoidal::lng() const
{
   return m_lon;
}

double LatLonEllipsoidal::longitude() const
{
   return m_lon;
}

void LatLonEllipsoidal::setLon(const double lon)
{
   m_lon = Dms::wrap180(lon);
}

void LatLonEllipsoidal::setLng(const double lon)
{
   m_lon = Dms::wrap180(lon);
}

void LatLonEllipsoidal::setLongitude(const double lon)
{
   m_lon = Dms::wrap180(lon);
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

bool LatLonEllipsoidal::equals(const LatLonEllipsoidal& point) const
{
   return *this == point;
}

Cartesian LatLonEllipsoidal::toCartesian() const
{
   // x = (ν+h)⋅cosφ⋅cosλ, y = (ν+h)⋅cosφ⋅sinλ, z = (ν⋅(1-e²)+h)⋅sinφ
   // where ν = a/√(1−e²⋅sinφ⋅sinφ), e² = (a²-b²)/a² or (better conditioned) 2⋅f-f²
   const Ellipsoid ellipsoid = m_datum
                     ? m_datum->ellipsoid
                     : m_referenceFrame ? m_referenceFrame->ellipsoid : ellipsoids().WGS84;

   const auto φ = toRadians(m_lat);
   const auto λ = toRadians(m_lon);
   const auto h = m_height;
   const auto a = ellipsoid.a, f = ellipsoid.f;

   const auto sinφ = std::sin(φ), cosφ = std::cos(φ);
   const auto sinλ = std::sin(λ), cosλ = std::cos(λ);

   const auto eSq = 2*f - f*f;                      // 1st eccentricity squared ≡ (a²-b²)/a²
   const auto ν = a / std::sqrt(1 - eSq*sinφ*sinφ); // radius of curvature in prime vertical

   const auto x = (ν+h) * cosφ * cosλ;
   const auto y = (ν+h) * cosφ * sinλ;
   const auto z = (ν*(1-eSq)+h) * sinφ;

   return { x, y, z };
}

std::string LatLonEllipsoidal::toString(Dms::eFormat format, std::optional<int> dph) const
{
   std::stringstream hwss;
   hwss << (m_height >= 0 ? L" +" : L" ");
   hwss << std::fixed << std::setprecision(dph.value_or(0)) << m_height << L"m";

   std::stringstream llwss;
   if (format == Dms::N)
   {
      // signed numeric degrees
      llwss << std::fixed << std::setprecision(4);
      llwss << m_lat << ",";
      llwss << m_lon;
      llwss << hwss.str();
      return llwss.str();
   }

   llwss << Dms::toLatitude(m_lat, format) << ", ";
   llwss << Dms::toLatitude(m_lon, format);
   llwss << (bool(dph) ? "" : hwss.str());

   return llwss.str();
}

