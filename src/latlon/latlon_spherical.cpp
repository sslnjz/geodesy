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
#include "latlon_spherical.h"

#include "vector3d.h"
#include "algorithm.h"

using namespace geodesy;

LatLonSpherical::LatLonSpherical(double lat, double lon)
   : m_lat(Dms::wrap90(lat))
   , m_lon(Dms::wrap180(lon))
{

}

double LatLonSpherical::lat() const
{
   return m_lat;
}

double LatLonSpherical::latitude() const
{
   return m_lat;
}

double LatLonSpherical::lon() const
{
   return m_lon;
}

double LatLonSpherical::lng() const
{
   return m_lon;
}

double LatLonSpherical::longitude() const
{
   return m_lon;
}

double LatLonSpherical::getMetresToKm()
{
   return 1.000 / 1000.000;
}

double LatLonSpherical::getMetresToMiles()
{
   return 1.000 / 1609.344;
}

double LatLonSpherical::getMetresToNauticalMiles()
{
   return 1.000 / 1852.000;
}

bool geodesy::operator==(const LatLonSpherical& p1, const LatLonSpherical& p2)
{
   return std::fabs(p1.m_lat - p2.m_lat) < std::numeric_limits<double>::epsilon() &&
      std::fabs(p1.m_lon - p2.m_lon) < std::numeric_limits<double>::epsilon();
}

bool geodesy::operator!=(const LatLonSpherical& p1, const LatLonSpherical& p2)
{
   return (!(p1 == p2));
}

void LatLonSpherical::setLat(double lat)
{
   m_lat = Dms::wrap90(lat);
}

void LatLonSpherical::setLatitude(double lat)
{
   m_lat = Dms::wrap90(lat);
}

void LatLonSpherical::setLon(double lon)
{
   m_lon = Dms::wrap180(lon);
}

void LatLonSpherical::setLng(double lon)
{
   m_lon = Dms::wrap180(lon);
}

void LatLonSpherical::setLongitude(double lon)
{
   m_lon = Dms::wrap180(lon);
}

double LatLonSpherical::distanceTo(const LatLonSpherical& point, double radius) const
{
   // a = sin²(Δφ/2) + cos(φ1)⋅cos(φ2)⋅sin²(Δλ/2)
   // δ = 2·atan2(√(a), √(1−a))
   // see mathforum.org/library/drmath/view/51879.html for derivation
   const auto φ1 = toRadians(m_lat), λ1 = toRadians(m_lon);
   const auto φ2 = toRadians(point.lat()), λ2 = toRadians(point.lon());
   const auto Δφ = φ2 - φ1;
   const auto Δλ = λ2 - λ1;
   const auto a = std::sin(Δφ / 2) * std::sin(Δφ / 2) + std::cos(φ1) * std::cos(φ2) * std::sin(Δλ / 2) *
      std::sin(Δλ / 2);
   const double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
   const double d = radius * c;
   return d;
}

double LatLonSpherical::initialBearingTo(const LatLonSpherical& point) const
{
   if (point == *this)
      return std::numeric_limits<double>::infinity() * 0.0; // coincident points
   // tanθ = sinΔλ⋅cosφ2 / cosφ1⋅sinφ2 − sinφ1⋅cosφ2⋅cosΔλ
   // see mathforum.org/library/drmath/view/55417.html for derivation
   const auto φ1 = toRadians(m_lat);
   const auto φ2 = toRadians(point.m_lat);
   const auto Δλ = toRadians(point.m_lon - m_lon);
   const auto x = std::cos(φ1) * std::sin(φ2) - std::sin(φ1) * std::cos(φ2) * std::cos(Δλ);
   const auto y = std::sin(Δλ) * std::cos(φ2);
   const auto θ = std::atan2(y, x);
   const auto bearing = toDegrees(θ);
   return Dms::wrap360(bearing);
}

double LatLonSpherical::finalBearingTo(const LatLonSpherical& point) const
{
   // get initial bearing from destination point to this point & reverse it by adding 180°
   const double bearing = point.initialBearingTo(*this) + 180;
   return Dms::wrap360(bearing);
}

LatLonSpherical LatLonSpherical::midpointTo(const LatLonSpherical& point) const
{
   // φm = atan2( sinφ1 + sinφ2, √( (cosφ1 + cosφ2⋅cosΔλ)² + cos²φ2⋅sin²Δλ ) )
   // λm = λ1 + atan2(cosφ2⋅sinΔλ, cosφ1 + cosφ2⋅cosΔλ)
   // midpoint is sum of vectors to two points: mathforum.org/library/drmath/view/51822.html
   const auto φ1 = toRadians(m_lat);
   const auto λ1 = toRadians(m_lon);
   const auto φ2 = toRadians(point.m_lat);
   const auto Δλ = toRadians((point.m_lon - m_lon));

   // get cartesian coordinates for the two points
   const vector3d A = {std::cos(φ1), 0, std::sin(φ1)}; // place point A on prime meridian y=0
   const vector3d B = {std::cos(φ2) * std::cos(Δλ), std::cos(φ2) * std::sin(Δλ), std::sin(φ2)};
   // vector to midpoint is sum of vectors to two points (no need to normalize)
   const vector3d C = A + B;
   const auto φm = std::atan2(C.z(), std::sqrt(C.x() * C.x() + C.y() * C.y()));
   const auto λm = λ1 + std::atan2(C.y(), C.x());
   const auto lat = toDegrees(φm);
   const auto lon = toDegrees(λm);
   return { lat, lon };
}

LatLonSpherical LatLonSpherical::intermediatePointTo(const LatLonSpherical& point, double fraction) const
{
   if (*this == point)
   {
      return {m_lat, m_lon}; // coincident points
   }
   const auto φ1 = toRadians(m_lat), λ1 = toRadians(m_lon);
   const auto φ2 = toRadians(point.m_lat), λ2 = toRadians(point.m_lon);
   // distance between points
   const auto Δφ = φ2 - φ1;
   const auto Δλ = λ2 - λ1;
   const auto a = std::sin(Δφ / 2) * std::sin(Δφ / 2)
      + std::cos(φ1) * std::cos(φ2) * std::sin(Δλ / 2) * std::sin(Δλ / 2);
   const auto δ = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
   const auto A = std::sin((1 - fraction) * δ) / std::sin(δ);
   const auto B = std::sin(fraction * δ) / std::sin(δ);
   const auto x = A * std::cos(φ1) * std::cos(λ1) + B * std::cos(φ2) * std::cos(λ2);
   const auto y = A * std::cos(φ1) * std::sin(λ1) + B * std::cos(φ2) * std::sin(λ2);
   const auto z = A * std::sin(φ1) + B * std::sin(φ2);
   const auto φ3 = std::atan2(z, std::sqrt(x * x + y * y));
   const auto λ3 = std::atan2(y, x);
   return { toDegrees(φ3), toDegrees(λ3) };
}

LatLonSpherical LatLonSpherical::destinationPoint(double distance, double bearing, double radius) const
{
   // sinφ2 = sinφ1⋅cosδ + cosφ1⋅sinδ⋅cosθ
   // tanΔλ = sinθ⋅sinδ⋅cosφ1 / cosδ−sinφ1⋅sinφ2
   // see mathforum.org/library/drmath/view/52049.html for derivation
   const auto δ = distance / radius; // angular distance in radians
   const auto θ = toRadians(bearing);
   const auto φ1 = toRadians(m_lat), λ1 = toRadians(m_lon);
   const auto sinφ2 = std::sin(φ1) * std::cos(δ) + std::cos(φ1) * std::sin(δ) * std::cos(θ);
   const auto φ2 = std::asin(sinφ2);
   const auto y =std::sin(θ) * std::sin(δ) * std::cos(φ1);
   const auto x =std::cos(δ) - std::sin(φ1) * sinφ2;
   const auto λ2 = λ1 + std::atan2(y, x);
   return { toDegrees(φ2), toDegrees(λ2) };
}

LatLonSpherical LatLonSpherical::intersection(const LatLonSpherical& p1, double brng1, const LatLonSpherical& p2, double brng2)
{
   // see www.edwilliams.org/avform.htm#Intersection
   const auto φ1 = toRadians(p1.lat()), λ1 = toRadians(p1.lon());
   const auto φ2 = toRadians(p2.lat()), λ2 = toRadians(p2.lon());
   const auto θ13 = toRadians(brng1), θ23 = toRadians(brng2);
   const auto Δφ = φ2 - φ1, Δλ = λ2 - λ1;
   // angular distance p1-p2
   const auto δ12 = 2 * std::asin(std::sqrt(std::sin(Δφ / 2) * std::sin(Δφ / 2)
      + std::cos(φ1) * std::cos(φ2) * std::sin(Δλ / 2) * std::sin(Δλ / 2)));
   if (std::abs(δ12) < std::numeric_limits<double>::epsilon()) 
   {
      return { p1.lat(), p1.lon() }; // coincident points
   }

   // initial/final bearings between points
   const auto cosθa = (std::sin(φ2) - std::sin(φ1) * std::cos(δ12)) / (std::sin(δ12) * std::cos(φ1));
   const auto cosθb = (std::sin(φ1) - std::sin(φ2) * std::cos(δ12)) / (std::sin(δ12) * std::cos(φ2));
   const auto θa = std::acos(std::fmin(std::fmax(cosθa, -1), 1)); // protect against rounding errors
   const auto θb = std::acos(std::fmin(std::fmax(cosθb, -1), 1)); // protect against rounding errors
   const auto θ12 = std::sin(λ2 - λ1) > 0 ? θa : 2 * π - θa;
   const auto θ21 = std::sin(λ2 - λ1) > 0 ? 2 * π - θb : θb;
   const auto α1 = θ13 - θ12; // angle 2-1-3
   const auto α2 = θ21 - θ23; // angle 1-2-3

   if (std::sin(α1) < std::numeric_limits<double>::epsilon()
      && std::sin(α2) < std::numeric_limits<double>::epsilon()) 
   {
      throw std::runtime_error("infinite intersections"); // infinite intersections
   }

   if (std::sin(α1) * std::sin(α2) < 0) 
   {
      throw std::runtime_error("ambiguous intersection (antipodal?)"); // ambiguous intersection (antipodal?)
   }

   const auto cosα3 = -std::cos(α1) * std::cos(α2) + std::sin(α1) * std::sin(α2) * std::cos(δ12);
   const auto δ13 = std::atan2(std::sin(δ12) * std::sin(α1) * std::sin(α2), std::cos(α2) + std::cos(α1) * cosα3);
   const auto φ3 = std::asin(std::sin(φ1) * std::cos(δ13) + std::cos(φ1) * std::sin(δ13) * std::cos(θ13));
   const auto Δλ13 = std::atan2(std::sin(θ13) * std::sin(δ13) * std::cos(φ1), std::cos(δ13) - std::sin(φ1) * std::sin(φ3));
   const auto λ3 = λ1 + Δλ13;
   return { toDegrees(φ3), toDegrees(λ3) };
}

double LatLonSpherical::crossTrackDistanceTo(const LatLonSpherical& pathStart, const LatLonSpherical& pathEnd,
                                             double radius) const
{
   const auto R = radius;
   const auto δ13 = pathStart.distanceTo(*this, R) / R;
   const auto θ13 = toRadians(pathStart.initialBearingTo(*this));
   const auto θ12 = toRadians(pathStart.initialBearingTo(pathEnd));
   const auto δxt = std::asin(std::sin(δ13) * std::sin(θ13 - θ12));
   return δxt * R;
}

double LatLonSpherical::alongTrackDistanceTo(const LatLonSpherical& pathStart, const LatLonSpherical& pathEnd,
                                             double radius) const
{
   const auto R = radius;
   const auto δ13 = pathStart.distanceTo(*this, R) / R;
   const auto θ13 = toRadians(pathStart.initialBearingTo(*this));
   const auto θ12 = toRadians(pathStart.initialBearingTo(pathEnd));
   const auto δxt = std::asin(std::sin(δ13) * std::sin(θ13 - θ12));
   const auto δat = std::acos(std::cos(δ13) / std::abs(std::cos(δxt)));
   return δat * sign(std::cos(θ12 - θ13)) * R;
}

double LatLonSpherical::maxLatitude(double bearing) const
{
   const auto θ = toRadians(bearing);
   const auto φ = toRadians(m_lat);
   const auto φMax = std::acos(std::abs(std::sin(θ) * std::cos(φ)));
   return toDegrees(φMax);
}

std::pair<double, double> LatLonSpherical::crossingParallels(const LatLonSpherical& point1,
                                                             const LatLonSpherical& point2, double latitude)
{
   if (point1 == (point2)) 
   {
      return {}; // coincident points
   }
   const auto φ  = toRadians(latitude);
   const auto φ1 = toRadians(point1.lat());
   const auto λ1 = toRadians(point1.lon());
   const auto φ2 = toRadians(point2.lat());
   const auto λ2 = toRadians(point2.lon());
   const auto Δλ = λ2 - λ1;
   const auto x = std::sin(φ1) * std::cos(φ2) * std::cos(φ) * std::sin(Δλ);
   const auto y = std::sin(φ1) * std::cos(φ2) * std::cos(φ) * std::cos(Δλ) - std::cos(φ1) * std::sin(φ2) * std::cos(φ);
   const double z = std::cos(φ1) * std::cos(φ2) * std::sin(φ) * std::sin(Δλ);

   if (z * z > x * x + y * y) 
   {
      return {}; // great circle doesn't reach latitude
   }

   const auto λm = std::atan2(-y, x); // longitude at max latitude
   const auto Δλi = std::acos(z / std::sqrt(x * x + y * y)); // Δλ from λm to intersection points
   const auto λi1 = λ1 + λm - Δλi;
   const auto λi2 = λ1 + λm + Δλi;
   const auto lon1 = toDegrees(λi1);
   const auto lon2 = toDegrees(λi2);
   return { Dms::wrap180(lon1), Dms::wrap180(lon2) };
}

double LatLonSpherical::rhumbDistanceTo(const LatLonSpherical& point, double radius) const
{
   // see www.edwilliams.org/avform.htm#Rhumb
   const auto R = radius;
   const auto φ1 = toRadians(m_lat);
   const auto φ2 = toRadians(point.lat());
   const auto Δφ = φ2 - φ1;
   auto Δλ = toRadians(std::abs(point.lon() - m_lon));
   // if dLon over 180° take shorter rhumb line across the anti-meridian:
   if (std::abs(Δλ) > π) 
      Δλ = Δλ > 0 ? -(2 * π - Δλ) : (2 * π + Δλ);
   // on Mercator projection, longitude distances shrink by latitude; q is the 'stretch factor'
   // q becomes ill-conditioned along E-W line (0/0); use empirical tolerance to avoid it
   const auto Δψ = std::log(std::tan(φ2 / 2 + π / 4) / std::tan(φ1 / 2 + π / 4));
   const auto q = std::abs(Δψ) > 10e-12 ? Δφ / Δψ : std::cos(φ1);
   // distance is pythagoras on 'stretched' Mercator projection, √(Δφ² + q²·Δλ²)
   const auto δ = std::sqrt(Δφ * Δφ + q * q * Δλ * Δλ); // angular distance in radians
   const auto d = δ * R;
   return d;
}

double LatLonSpherical::rhumbBearingTo(const LatLonSpherical& point) const
{
   if (*this == point) 
      return std::numeric_limits<double>::infinity() * 0.0; // coincident points
   const auto φ1 = toRadians(m_lat);
   const auto φ2 = toRadians(point.lat());
   auto Δλ = toRadians((point.lon() - m_lon));
   // if dLon over 180° take shorter rhumb line across the anti-meridian:
   if (std::abs(Δλ) > π)
      Δλ = Δλ > 0 ? -(2 * π - Δλ) : (2 * π + Δλ);
   const auto Δψ = std::log(std::tan(φ2 / 2 + π / 4) / std::tan(φ1 / 2 + π / 4));
   const auto θ = std::atan2(Δλ, Δψ);
   const auto bearing = toDegrees(θ);
   return Dms::wrap360(bearing);
}

LatLonSpherical LatLonSpherical::rhumbDestinationPoint(double distance, double bearing, double radius) const
{
   const auto φ1 = toRadians(m_lat);
   const auto λ1 = toRadians(m_lon);
   const auto θ = toRadians(bearing);
   const auto δ = distance / radius; // angular distance in radians
   const auto Δφ = δ * std::cos(θ);
   double φ2 = φ1 + Δφ;
   // check for some daft bugger going past the pole, normalise latitude if so
   if (std::abs(φ2) > π / 2) 
   {
      φ2 = φ2 > 0 ? π - φ2 : -π - φ2;
   }
   const auto Δψ = std::log(std::tan(φ2 / 2 + π / 4) / std::tan(φ1 / 2 + π / 4));
   const auto q = std::abs(Δψ) > 10e-12 ? Δφ / Δψ : std::cos(φ1); // E-W course becomes ill-conditioned with 0/0
   const auto Δλ = δ * std::sin(θ) / q;
   const auto λ2 = λ1 + Δλ;
   return { toDegrees(φ2), toDegrees(λ2) };
}

LatLonSpherical LatLonSpherical::rhumbMidpointTo(const LatLonSpherical& point) const
{
   // see mathforum.org/kb/message.jspa?messageID=148837
   const auto φ1 = toRadians(m_lat);
   double λ1 = toRadians(m_lon);
   const auto φ2 = toRadians(point.lat());
   const auto λ2 = toRadians(point.lon());
   if (std::abs(λ2 - λ1) > π)
      λ1 += 2 * π; // crossing anti-meridian
   const auto φ3 = (φ1 + φ2) / 2;
   const auto f1 = std::tan(π / 4 + φ1 / 2);
   const auto f2 = std::tan(π / 4 + φ2 / 2);
   const auto f3 = std::tan(π / 4 + φ3 / 2);
   double λ3 = ((λ2 - λ1) * std::log(f3) + λ1 * std::log(f2) - λ2 * std::log(f1)) / std::log(f2 / f1);
   if (!std::isfinite(λ3)) λ3 = (λ1 + λ2) / 2; // parallel of latitude
   const auto lat = toDegrees(φ3);
   const auto lon = toDegrees(λ3);
   return {lat, lon};
}

double LatLonSpherical::areaOf(std::vector<LatLonSpherical>& polygon, double radius)
{
   if (polygon.empty()) 
      return std::numeric_limits<double>::infinity() * 0.0;

   // uses method due to Karney: osgeo-org.1560.x6.nabble.com/Area-of-a-spherical-polygon-td3841625.html;
   // for each edge of the polygon, tan(E/2) = tan(Δλ/2)·(tan(φ₁/2)+tan(φ₂/2)) / (1+tan(φ₁/2)·tan(φ₂/2))
   // where E is the spherical excess of the trapezium obtained by extending the edge to the equator
   // (Karney's method is probably more efficient than the more widely known L’Huilier’s Theorem)
   const auto R = radius;
   // close polygon so that last point equals first point
   const auto closed = polygon[0] == (polygon[polygon.size() - 1]);
   if (!closed) {
      polygon.push_back(polygon[0]);
   }
   const size_t nVertices = polygon.size() - 1;
   auto S = 0.0; // spherical excess in steradians
   for (size_t v = 0; v < nVertices; ++v) 
   {
      const auto φ1 = toRadians(polygon[v].lat());
      const auto φ2 = toRadians(polygon[v + 1].lat());
      const auto Δλ = toRadians((polygon[v + 1].lon() - polygon[v].lon()));
      const auto E = 2 * std::atan2(std::tan(Δλ / 2) * (std::tan(φ1 / 2) + std::tan(φ2 / 2)), 1 + std::tan(φ1 / 2) * std::tan(φ2 / 2));
      S += E;
   }

   if (auto isPoleEnclosedBy = 
      [](const std::vector<LatLonSpherical>& p){
         //TODO: any better test than this
         double ΣΔ = 0.0;
         double prevBrng = p[0].initialBearingTo(p[1]);
         for (size_t v = 0; v < p.size() - 1; v++)
         {
            const auto initBrng = p[v].initialBearingTo(p[v + 1]);
            const auto finalBrng = p[v].finalBearingTo(p[v + 1]);
            ΣΔ += std::fmod((initBrng - prevBrng + 540), 360) - 180;
            ΣΔ += std::fmod((finalBrng - initBrng + 540), 360) - 180;
            prevBrng = finalBrng;
         }
         const auto initBrng = p[0].initialBearingTo(p[1]);
         ΣΔ += std::fmod((initBrng - prevBrng + 540), 360) - 180;
         // TODO: fix (intermittant) edge crossing pole - eg (85,90), (85,0), (85,-90)
         return (std::abs(ΣΔ) < 90);
      }; isPoleEnclosedBy(polygon)) 
   {
      S = std::fabs(S) - 2 * π;
   }

   const auto A = std::abs(S * R * R); // area in units of R
   if (!closed) 
   {
      polygon.pop_back(); // restore polygon to pristine condition
   }
   return A;
}

std::string LatLonSpherical::toString(Dms::eFormat e) const
{
   // note: explicitly set dp to undefined for passing through to toLat/toLon
   const std::string lat = Dms::toLatitude(m_lat, e);
   const std::string lon = Dms::toLongitude(m_lon, e);
   return lat + "," + lon;
}

std::string LatLonSpherical::toGeoJSON() const
{
   return std::string("{ type: \"Point\", coordinates : [")
      + std::to_string(m_lon) + std::string(",")
      + std::to_string(m_lat) + std::string("] }");
}
