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
   const auto phi1 = toRadians(m_lat), lambda1 = toRadians(m_lon);
   const auto phi2 = toRadians(point.lat()), lambda2 = toRadians(point.lon());
   const auto DELTAphi = phi2 - phi1;
   const auto DELTAlambda = lambda2 - lambda1;
   const auto a = std::sin(DELTAphi / 2) * std::sin(DELTAphi / 2) + std::cos(phi1) * std::cos(phi2) * std::sin(DELTAlambda / 2) *
      std::sin(DELTAlambda / 2);
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
   const auto phi1 = toRadians(m_lat);
   const auto phi2 = toRadians(point.m_lat);
   const auto DELTAlambda = toRadians(point.m_lon - m_lon);
   const auto x = std::cos(phi1) * std::sin(phi2) - std::sin(phi1) * std::cos(phi2) * std::cos(DELTAlambda);
   const auto y = std::sin(DELTAlambda) * std::cos(phi2);
   const auto theta = std::atan2(y, x);
   const auto bearing = toDegrees(theta);
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
   const auto phi1 = toRadians(m_lat);
   const auto lambda1 = toRadians(m_lon);
   const auto phi2 = toRadians(point.m_lat);
   const auto DELTAlambda = toRadians((point.m_lon - m_lon));

   // get cartesian coordinates for the two points
   const vector3d A = {std::cos(phi1), 0, std::sin(phi1)}; // place point A on prime meridian y=0
   const vector3d B = {std::cos(phi2) * std::cos(DELTAlambda), std::cos(phi2) * std::sin(DELTAlambda), std::sin(phi2)};
   // vector to midpoint is sum of vectors to two points (no need to normalize)
   const vector3d C = A + B;
   const auto phim = std::atan2(C.z(), std::sqrt(C.x() * C.x() + C.y() * C.y()));
   const auto lambdam = lambda1 + std::atan2(C.y(), C.x());
   const auto lat = toDegrees(phim);
   const auto lon = toDegrees(lambdam);
   return { lat, lon };
}

LatLonSpherical LatLonSpherical::intermediatePointTo(const LatLonSpherical& point, double fraction) const
{
   if (*this == point)
   {
      return {m_lat, m_lon}; // coincident points
   }
   const auto phi1 = toRadians(m_lat), lambda1 = toRadians(m_lon);
   const auto phi2 = toRadians(point.m_lat), lambda2 = toRadians(point.m_lon);
   // distance between points
   const auto DELTAphi = phi2 - phi1;
   const auto DELTAlambda = lambda2 - lambda1;
   const auto a = std::sin(DELTAphi / 2) * std::sin(DELTAphi / 2)
      + std::cos(phi1) * std::cos(phi2) * std::sin(DELTAlambda / 2) * std::sin(DELTAlambda / 2);
   const auto delta = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
   const auto A = std::sin((1 - fraction) * delta) / std::sin(delta);
   const auto B = std::sin(fraction * delta) / std::sin(delta);
   const auto x = A * std::cos(phi1) * std::cos(lambda1) + B * std::cos(phi2) * std::cos(lambda2);
   const auto y = A * std::cos(phi1) * std::sin(lambda1) + B * std::cos(phi2) * std::sin(lambda2);
   const auto z = A * std::sin(phi1) + B * std::sin(phi2);
   const auto phi3 = std::atan2(z, std::sqrt(x * x + y * y));
   const auto lambda3 = std::atan2(y, x);
   return { toDegrees(phi3), toDegrees(lambda3) };
}

LatLonSpherical LatLonSpherical::destinationPoint(double distance, double bearing, double radius) const
{
   // sinφ2 = sinφ1⋅cosδ + cosφ1⋅sinδ⋅cosθ
   // tanΔλ = sinθ⋅sinδ⋅cosφ1 / cosδ−sinφ1⋅sinφ2
   // see mathforum.org/library/drmath/view/52049.html for derivation
   const auto delta = distance / radius; // angular distance in radians
   const auto theta = toRadians(bearing);
   const auto phi1 = toRadians(m_lat), lambda1 = toRadians(m_lon);
   const auto sinphi2 = std::sin(phi1) * std::cos(delta) + std::cos(phi1) * std::sin(delta) * std::cos(theta);
   const auto phi2 = std::asin(sinphi2);
   const auto y =std::sin(theta) * std::sin(delta) * std::cos(phi1);
   const auto x =std::cos(delta) - std::sin(phi1) * sinphi2;
   const auto lambda2 = lambda1 + std::atan2(y, x);
   return { toDegrees(phi2), toDegrees(lambda2) };
}

LatLonSpherical LatLonSpherical::intersection(const LatLonSpherical& p1, double brng1, const LatLonSpherical& p2, double brng2)
{
   // see www.edwilliams.org/avform.htm#Intersection
   const auto phi1 = toRadians(p1.lat()), lambda1 = toRadians(p1.lon());
   const auto phi2 = toRadians(p2.lat()), lambda2 = toRadians(p2.lon());
   const auto theta13 = toRadians(brng1), theta23 = toRadians(brng2);
   const auto DELTAphi = phi2 - phi1, DELTAlambda = lambda2 - lambda1;
   // angular distance p1-p2
   const auto delta12 = 2 * std::asin(std::sqrt(std::sin(DELTAphi / 2) * std::sin(DELTAphi / 2)
      + std::cos(phi1) * std::cos(phi2) * std::sin(DELTAlambda / 2) * std::sin(DELTAlambda / 2)));
   if (std::abs(delta12) < std::numeric_limits<double>::epsilon()) 
   {
      return { p1.lat(), p1.lon() }; // coincident points
   }

   // initial/final bearings between points
   const auto costhetaa = (std::sin(phi2) - std::sin(phi1) * std::cos(delta12)) / (std::sin(delta12) * std::cos(phi1));
   const auto costhetab = (std::sin(phi1) - std::sin(phi2) * std::cos(delta12)) / (std::sin(delta12) * std::cos(phi2));
   const auto thetaa = std::acos(std::fmin(std::fmax(costhetaa, -1), 1)); // protect against rounding errors
   const auto thetab = std::acos(std::fmin(std::fmax(costhetab, -1), 1)); // protect against rounding errors
   const auto theta12 = std::sin(lambda2 - lambda1) > 0 ? thetaa : 2 * pi - thetaa;
   const auto theta21 = std::sin(lambda2 - lambda1) > 0 ? 2 * pi - thetab : thetab;
   const auto alpha1 = theta13 - theta12; // angle 2-1-3
   const auto alpha2 = theta21 - theta23; // angle 1-2-3

   if (std::sin(alpha1) < std::numeric_limits<double>::epsilon()
      && std::sin(alpha2) < std::numeric_limits<double>::epsilon()) 
   {
      throw std::runtime_error("infinite intersections"); // infinite intersections
   }

   if (std::sin(alpha1) * std::sin(alpha2) < 0) 
   {
      throw std::runtime_error("ambiguous intersection (antipodal?)"); // ambiguous intersection (antipodal?)
   }

   const auto cosalpha3 = -std::cos(alpha1) * std::cos(alpha2) + std::sin(alpha1) * std::sin(alpha2) * std::cos(delta12);
   const auto delta13 = std::atan2(std::sin(delta12) * std::sin(alpha1) * std::sin(alpha2), std::cos(alpha2) + std::cos(alpha1) * cosalpha3);
   const auto phi3 = std::asin(std::sin(phi1) * std::cos(delta13) + std::cos(phi1) * std::sin(delta13) * std::cos(theta13));
   const auto DELTAlambda13 = std::atan2(std::sin(theta13) * std::sin(delta13) * std::cos(phi1), std::cos(delta13) - std::sin(phi1) * std::sin(phi3));
   const auto lambda3 = lambda1 + DELTAlambda13;
   return { toDegrees(phi3), toDegrees(lambda3) };
}

double LatLonSpherical::crossTrackDistanceTo(const LatLonSpherical& pathStart, const LatLonSpherical& pathEnd,
                                             double radius) const
{
   const auto R = radius;
   const auto delta13 = pathStart.distanceTo(*this, R) / R;
   const auto theta13 = toRadians(pathStart.initialBearingTo(*this));
   const auto theta12 = toRadians(pathStart.initialBearingTo(pathEnd));
   const auto deltaxt = std::asin(std::sin(delta13) * std::sin(theta13 - theta12));
   return deltaxt * R;
}

double LatLonSpherical::alongTrackDistanceTo(const LatLonSpherical& pathStart, const LatLonSpherical& pathEnd,
                                             double radius) const
{
   const auto R = radius;
   const auto delta13 = pathStart.distanceTo(*this, R) / R;
   const auto theta13 = toRadians(pathStart.initialBearingTo(*this));
   const auto theta12 = toRadians(pathStart.initialBearingTo(pathEnd));
   const auto deltaxt = std::asin(std::sin(delta13) * std::sin(theta13 - theta12));
   const auto deltaat = std::acos(std::cos(delta13) / std::abs(std::cos(deltaxt)));
   return deltaat * sign(std::cos(theta12 - theta13)) * R;
}

double LatLonSpherical::maxLatitude(double bearing) const
{
   const auto theta = toRadians(bearing);
   const auto phi = toRadians(m_lat);
   const auto phiMax = std::acos(std::abs(std::sin(theta) * std::cos(phi)));
   return toDegrees(phiMax);
}

std::pair<double, double> LatLonSpherical::crossingParallels(const LatLonSpherical& point1,
                                                             const LatLonSpherical& point2, double latitude)
{
   if (point1 == (point2)) 
   {
      return {}; // coincident points
   }
   const auto phi = toRadians(latitude);
   const auto phi1 = toRadians(point1.lat());
   const auto lambda1 = toRadians(point1.lon());
   const auto phi2 = toRadians(point2.lat());
   const auto lambda2 = toRadians(point2.lon());
   const auto DELTAlambda = lambda2 - lambda1;
   const auto x = std::sin(phi1) * std::cos(phi2) * std::cos(phi) * std::sin(DELTAlambda);
   const auto y = std::sin(phi1) * std::cos(phi2) * std::cos(phi) * std::cos(DELTAlambda) - std::cos(phi1) * std::sin(phi2) * std::cos(phi);
   const double z = std::cos(phi1) * std::cos(phi2) * std::sin(phi) * std::sin(DELTAlambda);

   if (z * z > x * x + y * y) 
   {
      return {}; // great circle doesn't reach latitude
   }

   const auto lambdam = std::atan2(-y, x); // longitude at max latitude
   const auto DELTAlambdai = std::acos(z / std::sqrt(x * x + y * y)); // Δλ from λm to intersection points
   const auto lambdai1 = lambda1 + lambdam - DELTAlambdai;
   const auto lambdai2 = lambda1 + lambdam + DELTAlambdai;
   const auto lon1 = toDegrees(lambdai1);
   const auto lon2 = toDegrees(lambdai2);
   return { Dms::wrap180(lon1), Dms::wrap180(lon2) };
}

double LatLonSpherical::rhumbDistanceTo(const LatLonSpherical& point, double radius) const
{
   // see www.edwilliams.org/avform.htm#Rhumb
   const auto R = radius;
   const auto phi1 = toRadians(m_lat);
   const auto phi2 = toRadians(point.lat());
   const auto DELTAphi = phi2 - phi1;
   auto DELTAlambda = toRadians(std::abs(point.lon() - m_lon));
   // if dLon over 180° take shorter rhumb line across the anti-meridian:
   if (std::abs(DELTAlambda) > pi) 
      DELTAlambda = DELTAlambda > 0 ? -(2 * pi - DELTAlambda) : (2 * pi + DELTAlambda);
   // on Mercator projection, longitude distances shrink by latitude; q is the 'stretch factor'
   // q becomes ill-conditioned along E-W line (0/0); use empirical tolerance to avoid it
   const auto DELTApsi = std::log(std::tan(phi2 / 2 + pi / 4) / std::tan(phi1 / 2 + pi / 4));
   const auto q = std::abs(DELTApsi) > 10e-12 ? DELTAphi / DELTApsi : std::cos(phi1);
   // distance is pythagoras on 'stretched' Mercator projection, √(Δφ² + q²·Δλ²)
   const auto delta = std::sqrt(DELTAphi * DELTAphi + q * q * DELTAlambda * DELTAlambda); // angular distance in radians
   const auto d = delta * R;
   return d;
}

double LatLonSpherical::rhumbBearingTo(const LatLonSpherical& point) const
{
   if (*this == point) 
      return std::numeric_limits<double>::infinity() * 0.0; // coincident points
   const auto phi1 = toRadians(m_lat);
   const auto phi2 = toRadians(point.lat());
   auto DELTAlambda = toRadians((point.lon() - m_lon));
   // if dLon over 180° take shorter rhumb line across the anti-meridian:
   if (std::abs(DELTAlambda) > pi)
      DELTAlambda = DELTAlambda > 0 ? -(2 * pi - DELTAlambda) : (2 * pi + DELTAlambda);
   const auto DELTAψ = std::log(std::tan(phi2 / 2 + pi / 4) / std::tan(phi1 / 2 + pi / 4));
   const auto theta = std::atan2(DELTAlambda, DELTAψ);
   const auto bearing = toDegrees(theta);
   return Dms::wrap360(bearing);
}

LatLonSpherical LatLonSpherical::rhumbDestinationPoint(double distance, double bearing, double radius) const
{
   const auto phi1 = toRadians(m_lat);
   const auto lambda1 = toRadians(m_lon);
   const auto theta = toRadians(bearing);
   const auto delta = distance / radius; // angular distance in radians
   const auto DELTAphi = delta * std::cos(theta);
   double phi2 = phi1 + DELTAphi;
   // check for some daft bugger going past the pole, normalise latitude if so
   if (std::abs(phi2) > pi / 2) 
   {
      phi2 = phi2 > 0 ? pi - phi2 : -pi - phi2;
   }
   const auto DELTApsi = std::log(std::tan(phi2 / 2 + pi / 4) / std::tan(phi1 / 2 + pi / 4));
   const auto q = std::abs(DELTApsi) > 10e-12 ? DELTAphi / DELTApsi : std::cos(phi1); // E-W course becomes ill-conditioned with 0/0
   const auto DELTAlambda = delta * std::sin(theta) / q;
   const auto lambda2 = lambda1 + DELTAlambda;
   return { toDegrees(phi2), toDegrees(lambda2) };
}

LatLonSpherical LatLonSpherical::rhumbMidpointTo(const LatLonSpherical& point) const
{
   // see mathforum.org/kb/message.jspa?messageID=148837
   const auto phi1 = toRadians(m_lat);
   double lambda1 = toRadians(m_lon);
   const auto phi2 = toRadians(point.lat());
   const auto lambda2 = toRadians(point.lon());
   if (std::abs(lambda2 - lambda1) > pi)
      lambda1 += 2 * pi; // crossing anti-meridian
   const auto phi3 = (phi1 + phi2) / 2;
   const auto f1 = std::tan(pi / 4 + phi1 / 2);
   const auto f2 = std::tan(pi / 4 + phi2 / 2);
   const auto f3 = std::tan(pi / 4 + phi3 / 2);
   double lambda3 = ((lambda2 - lambda1) * std::log(f3) + lambda1 * std::log(f2) - lambda2 * std::log(f1)) / std::log(f2 / f1);
   if (!std::isfinite(lambda3)) lambda3 = (lambda1 + lambda2) / 2; // parallel of latitude
   const auto lat = toDegrees(phi3);
   const auto lon = toDegrees(lambda3);
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
      const auto phi1 = toRadians(polygon[v].lat());
      const auto phi2 = toRadians(polygon[v + 1].lat());
      const auto DELTAlambda = toRadians((polygon[v + 1].lon() - polygon[v].lon()));
      const auto E = 2 * std::atan2(std::tan(DELTAlambda / 2) * (std::tan(phi1 / 2) + std::tan(phi2 / 2)), 1 + std::tan(phi1 / 2) * std::tan(phi2 / 2));
      S += E;
   }

   if (auto isPoleEnclosedBy = 
      [](const std::vector<LatLonSpherical>& p){
         //TODO: any better test than this
         double ΣDELTA = 0.0;
         double prevBrng = p[0].initialBearingTo(p[1]);
         for (size_t v = 0; v < p.size() - 1; v++)
         {
            const auto initBrng = p[v].initialBearingTo(p[v + 1]);
            const auto finalBrng = p[v].finalBearingTo(p[v + 1]);
            ΣDELTA += std::fmod((initBrng - prevBrng + 540), 360) - 180;
            ΣDELTA += std::fmod((finalBrng - initBrng + 540), 360) - 180;
            prevBrng = finalBrng;
         }
         const auto initBrng = p[0].initialBearingTo(p[1]);
         ΣDELTA += std::fmod((initBrng - prevBrng + 540), 360) - 180;
         // TODO: fix (intermittant) edge crossing pole - eg (85,90), (85,0), (85,-90)
         return (std::abs(ΣDELTA) < 90);
      }; isPoleEnclosedBy(polygon)) 
   {
      S = std::fabs(S) - 2 * pi;
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
   const std::string lat = Dms::toLat(m_lat, e);
   const std::string lon = Dms::toLon(m_lon, e);
   return lat + "," + lon;
}

std::string LatLonSpherical::toGeoJSON() const
{
   return std::string("{ type: \"Point\", coordinates : [")
      + std::to_string(m_lon) + std::string(",")
      + std::to_string(m_lat) + std::string("] }");
}
