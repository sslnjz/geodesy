#include "latlon_spherical.h"

#include "vector3d.hpp"
#include "algorithm.hpp"

using namespace geodesy;

LatLonSpherical::LatLonSpherical(double lat, double lon)
   : m_lat(Dms::wrap90(lat))
   , m_lon(Dms::wrap180(lon))
{

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
   if (isnan(radius)) 
      throw std::invalid_argument("invalid radius: " + std::to_string(radius));

   // a = sin²(Δφ/2) + cos(φ1)⋅cos(φ2)⋅sin²(Δλ/2)
   // δ = 2·atan2(√(a), √(1−a))
   // see mathforum.org/library/drmath/view/51879.html for derivation
   const double R = radius;
   const double φ1 = toRadians(m_lat), λ1 = toRadians(m_lon);
   const double φ2 = toRadians(point.lat()), λ2 = toRadians(point.lon());
   const double Δφ = φ2 - φ1;
   const double Δλ = λ2 - λ1;
   const double a = std::sin(Δφ / 2) * std::sin(Δφ / 2) + std::cos(φ1) * std::cos(φ2) * std::sin(Δλ / 2) *
      std::sin(Δλ / 2);
   const double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
   const double d = R * c;
   return d;
}

double LatLonSpherical::initialBearingTo(const LatLonSpherical& point) const
{
   if (point == *this)
      return std::numeric_limits<double>::infinity() * 0.0; // coincident points
   // tanθ = sinΔλ⋅cosφ2 / cosφ1⋅sinφ2 − sinφ1⋅cosφ2⋅cosΔλ
   // see mathforum.org/library/drmath/view/55417.html for derivation
   const double φ1 = toRadians(m_lat);
   const double φ2 = toRadians(point.m_lat);
   const double Δλ = toRadians(point.m_lon - m_lon);
   const double x = std::cos(φ1) * std::sin(φ2) - std::sin(φ1) * std::cos(φ2) * std::cos(Δλ);
   const double y = std::sin(Δλ) * std::cos(φ2);
   const double θ = std::atan2(y, x);
   const double bearing = toDegrees(θ);
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
   const double φ1 = toRadians(m_lat);
   const double λ1 = toRadians(m_lon);
   const double φ2 = toRadians(point.m_lat);
   const double Δλ = toRadians((point.m_lon - m_lon));

   // get cartesian coordinates for the two points
   const vector3d A = {std::cos(φ1), 0, std::sin(φ1)}; // place point A on prime meridian y=0
   const vector3d B = {std::cos(φ2) * std::cos(Δλ), std::cos(φ2) * std::sin(Δλ), std::sin(φ2)};
   // vector to midpoint is sum of vectors to two points (no need to normalize)
   const vector3d C = A + B;
   const double φm = std::atan2(C.z(), std::sqrt(C.x() * C.x() + C.y() * C.y()));
   const double λm = λ1 + std::atan2(C.y(), C.x());
   const double lat = toDegrees(φm);
   const double lon = toDegrees(λm);
   return { lat, lon };
}

LatLonSpherical LatLonSpherical::intermediatePointTo(const LatLonSpherical& point, double fraction) const
{
   if (*this == point)
   {
      return {m_lat, m_lon}; // coincident points
   }
   const double φ1 = toRadians(m_lat), λ1 = toRadians(m_lon);
   const double φ2 = toRadians(point.m_lat), λ2 = toRadians(point.m_lon);
   // distance between points
   const double Δφ = φ2 - φ1;
   const double Δλ = λ2 - λ1;
   const double a = std::sin(Δφ / 2) * std::sin(Δφ / 2)
      + std::cos(φ1) * std::cos(φ2) * std::sin(Δλ / 2) * std::sin(Δλ / 2);
   const double δ = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
   const double A = std::sin((1 - fraction) * δ) / std::sin(δ);
   const double B = std::sin(fraction * δ) / std::sin(δ);
   const double x = A * std::cos(φ1) * std::cos(λ1) + B * std::cos(φ2) * std::cos(λ2);
   const double y = A * std::cos(φ1) * std::sin(λ1) + B * std::cos(φ2) * std::sin(λ2);
   const double z = A * std::sin(φ1) + B * std::sin(φ2);
   const double φ3 = std::atan2(z, std::sqrt(x * x + y * y));
   const double λ3 = std::atan2(y, x);
   return {toDegrees(φ3), toDegrees(λ3)};
}

LatLonSpherical LatLonSpherical::destinationPoint(double distance, double bearing, double radius) const
{
   // sinφ2 = sinφ1⋅cosδ + cosφ1⋅sinδ⋅cosθ
   // tanΔλ = sinθ⋅sinδ⋅cosφ1 / cosδ−sinφ1⋅sinφ2
   // see mathforum.org/library/drmath/view/52049.html for derivation
   const double δ = distance / radius; // angular distance in radians
   const double θ = toRadians(bearing);
   const double φ1 = toRadians(m_lat), λ1 = toRadians(m_lon);
   const double sinφ2 = std::sin(φ1) * std::cos(δ) + std::cos(φ1) * std::sin(δ) * std::cos(θ);
   const double φ2 = std::asin(sinφ2);
   const double y =std::sin(θ) * std::sin(δ) * std::cos(φ1);
   const double x =std::cos(δ) - std::sin(φ1) * sinφ2;
   const double λ2 = λ1 + std::atan2(y, x);
   return { toDegrees(φ2), toDegrees(λ2) };
}

LatLonSpherical LatLonSpherical::intersection(const LatLonSpherical& p1, double brng1, const LatLonSpherical& p2, double brng2)
{
   if (isnan(brng1)) throw std::invalid_argument("invalid brng1 ‘${ brng1 }’");
   if (isnan(brng2)) throw std::invalid_argument("invalid brng2 ‘${ brng2 }’");

   // see www.edwilliams.org/avform.htm#Intersection
   const double φ1 = toRadians(p1.lat()), λ1 = toRadians(p1.lon());
   const double φ2 = toRadians(p2.lat()), λ2 = toRadians(p2.lon());
   const double θ13 = toRadians(brng1), θ23 = toRadians(brng2);
   const double Δφ = φ2 - φ1, Δλ = λ2 - λ1;
   // angular distance p1-p2
   const double δ12 = 2 * std::asin(std::sqrt(std::sin(Δφ / 2) * std::sin(Δφ / 2)
      + std::cos(φ1) * std::cos(φ2) * std::sin(Δλ / 2) * std::sin(Δλ / 2)));
   if (std::abs(δ12) < std::numeric_limits<double>::epsilon()) {
      return { p1.lat(), p1.lon() }; // coincident points
   }

   // initial/final bearings between points
   const double cosθa = (std::sin(φ2) - std::sin(φ1) * std::cos(δ12)) / (std::sin(δ12) * std::cos(φ1));
   const double cosθb = (std::sin(φ1) - std::sin(φ2) * std::cos(δ12)) / (std::sin(δ12) * std::cos(φ2));
   const double θa = std::acos(std::fmin(std::fmax(cosθa, -1), 1)); // protect against rounding errors
   const double θb = std::acos(std::fmin(std::fmax(cosθb, -1), 1)); // protect against rounding errors
   const double θ12 = std::sin(λ2 - λ1) > 0 ? θa : 2 * π - θa;
   const double θ21 = std::sin(λ2 - λ1) > 0 ? 2 * π - θb : θb;
   const double α1 = θ13 - θ12; // angle 2-1-3
   const double α2 = θ21 - θ23; // angle 1-2-3

   if (std::sin(α1) < std::numeric_limits<double>::epsilon()
      && std::sin(α2) < std::numeric_limits<double>::epsilon()) {
      throw std::runtime_error("infinite intersections"); // infinite intersections
   }

   if (std::sin(α1) * std::sin(α2) < 0) {
      throw std::runtime_error("ambiguous intersection (antipodal?)"); // ambiguous intersection (antipodal?)
   }

   const double cosα3 = -std::cos(α1) * std::cos(α2) + std::sin(α1) * std::sin(α2) * std::cos(δ12);
   const double δ13 = std::atan2(std::sin(δ12) * std::sin(α1) * std::sin(α2), std::cos(α2) + std::cos(α1) * cosα3);
   const double φ3 = std::asin(std::sin(φ1) * std::cos(δ13) + std::cos(φ1) * std::sin(δ13) * std::cos(θ13));
   const double Δλ13 = std::atan2(std::sin(θ13) * std::sin(δ13) * std::cos(φ1), std::cos(δ13) - std::sin(φ1) * std::sin(φ3));
   const double λ3 = λ1 + Δλ13;
   return { toDegrees(φ3), toDegrees(λ3) };
}

double LatLonSpherical::crossTrackDistanceTo(const LatLonSpherical& pathStart, const LatLonSpherical& pathEnd,
                                             double radius) const
{
   const double R = radius;
   const double δ13 = pathStart.distanceTo(*this, R) / R;
   const double θ13 = toRadians(pathStart.initialBearingTo(*this));
   const double θ12 = toRadians(pathStart.initialBearingTo(pathEnd));
   const double δxt = std::asin(std::sin(δ13) * std::sin(θ13 - θ12));
   return δxt * R;
}

double LatLonSpherical::alongTrackDistanceTo(const LatLonSpherical& pathStart, const LatLonSpherical& pathEnd,
                                             double radius) const
{
   const double R = radius;
   const double δ13 = pathStart.distanceTo(*this, R) / R;
   const double θ13 = toRadians(pathStart.initialBearingTo(*this));
   const double θ12 = toRadians(pathStart.initialBearingTo(pathEnd));
   const double δxt = std::asin(std::sin(δ13) * std::sin(θ13 - θ12));
   const double δat = std::acos(std::cos(δ13) / std::abs(std::cos(δxt)));
   return δat * sign(std::cos(θ12 - θ13)) * R;
}

double LatLonSpherical::maxLatitude(double bearing) const
{
   const double θ = toRadians(bearing);
   const double φ = toRadians(m_lat);
   const double φMax = std::acos(std::abs(std::sin(θ) * std::cos(φ)));
   return toDegrees(φMax);
}

std::pair<double, double> LatLonSpherical::crossingParallels(const LatLonSpherical& point1,
                                                             const LatLonSpherical& point2, double latitude)
{
   if (point1 == (point2)) {
      return {}; // coincident points
   }
   const double φ = toRadians(latitude);
   const double φ1 = toRadians(point1.lat());
   const double λ1 = toRadians(point1.lon());
   const double φ2 = toRadians(point2.lat());
   const double λ2 = toRadians(point2.lon());
   const double Δλ = λ2 - λ1;
   const double x = std::sin(φ1) * std::cos(φ2) * std::cos(φ) * std::sin(Δλ);
   const double y = std::sin(φ1) * std::cos(φ2) * std::cos(φ) * std::cos(Δλ) - std::cos(φ1) * std::sin(φ2) *
      std::cos(φ);
   const double z = std::cos(φ1) * std::cos(φ2) * std::sin(φ) * std::sin(Δλ);

   if (z * z > x * x + y * y) {
      return {}; // great circle doesn't reach latitude
   }

   const double λm = std::atan2(-y, x); // longitude at max latitude
   const double Δλi = std::acos(z / std::sqrt(x * x + y * y)); // Δλ from λm to intersection points
   const double λi1 = λ1 + λm - Δλi;
   const double λi2 = λ1 + λm + Δλi;
   const double lon1 = toDegrees(λi1);
   const double lon2 = toDegrees(λi2);
   return {Dms::wrap180(lon1), Dms::wrap180(lon2)};
}

double LatLonSpherical::rhumbDistanceTo(const LatLonSpherical& point, double radius) const
{
   // see www.edwilliams.org/avform.htm#Rhumb
   const double R = radius;
   const double φ1 = toRadians(m_lat);
   const double φ2 = toRadians(point.lat());
   const double Δφ = φ2 - φ1;
   double Δλ = toRadians(std::abs(point.lon() - m_lon));
   // if dLon over 180° take shorter rhumb line across the anti-meridian:
   if (std::abs(Δλ) > π) 
      Δλ = Δλ > 0 ? -(2 * π - Δλ) : (2 * π + Δλ);
   // on Mercator projection, longitude distances shrink by latitude; q is the 'stretch factor'
   // q becomes ill-conditioned along E-W line (0/0); use empirical tolerance to avoid it
   const double Δψ = std::log(std::tan(φ2 / 2 + π / 4) / std::tan(φ1 / 2 + π / 4));
   const double q = std::abs(Δψ) > 10e-12 ? Δφ / Δψ : std::cos(φ1);
   // distance is pythagoras on 'stretched' Mercator projection, √(Δφ² + q²·Δλ²)
   const double δ = std::sqrt(Δφ * Δφ + q * q * Δλ * Δλ); // angular distance in radians
   const double d = δ * R;
   return d;
}

double LatLonSpherical::rhumbBearingTo(const LatLonSpherical& point) const
{
   if (*this == point) 
      return std::numeric_limits<double>::infinity() * 0.0; // coincident points
   const double φ1 = toRadians(m_lat);
   const double φ2 = toRadians(point.lat());
   double Δλ = toRadians((point.lon() - m_lon));
   // if dLon over 180° take shorter rhumb line across the anti-meridian:
   if (std::abs(Δλ) > π)
      Δλ = Δλ > 0 ? -(2 * π - Δλ) : (2 * π + Δλ);
   const double Δψ = std::log(std::tan(φ2 / 2 + π / 4) / std::tan(φ1 / 2 + π / 4));
   const double θ = std::atan2(Δλ, Δψ);
   const double bearing = toDegrees(θ);
   return Dms::wrap360(bearing);
}

LatLonSpherical LatLonSpherical::rhumbDestinationPoint(double distance, double bearing, double radius) const
{
   const double φ1 = toRadians(m_lat);
   const double λ1 = toRadians(m_lon);
   const double θ = toRadians(bearing);
   const double δ = distance / radius; // angular distance in radians
   const double Δφ = δ * std::cos(θ);
   double φ2 = φ1 + Δφ;
   // check for some daft bugger going past the pole, normalise latitude if so
   if (std::abs(φ2) > π / 2)
      φ2 = φ2 > 0 ? π - φ2 : -π - φ2;
   const double Δψ = std::log(std::tan(φ2 / 2 + π / 4) / std::tan(φ1 / 2 + π / 4));
   const double q = std::abs(Δψ) > 10e-12 ? Δφ / Δψ : std::cos(φ1); // E-W course becomes ill-conditioned with 0/0
   const double Δλ = δ * std::sin(θ) / q;
   const double λ2 = λ1 + Δλ;
   return { toDegrees(φ2), toDegrees(λ2) };
}

LatLonSpherical LatLonSpherical::rhumbMidpointTo(const LatLonSpherical& point) const
{
   // see mathforum.org/kb/message.jspa?messageID=148837
   const double φ1 = toRadians(m_lat);
   double λ1 = toRadians(m_lon);
   const double φ2 = toRadians(point.lat());
   const double λ2 = toRadians(point.lon());
   if (std::abs(λ2 - λ1) > π)
      λ1 += 2 * π; // crossing anti-meridian
   const double φ3 = (φ1 + φ2) / 2;
   const double f1 = std::tan(π / 4 + φ1 / 2);
   const double f2 = std::tan(π / 4 + φ2 / 2);
   const double f3 = std::tan(π / 4 + φ3 / 2);
   double λ3 = ((λ2 - λ1) * std::log(f3) + λ1 * std::log(f2) - λ2 * std::log(f1)) / std::log(f2 / f1);
   if (!std::isfinite(λ3)) λ3 = (λ1 + λ2) / 2; // parallel of latitude
   const double lat = toDegrees(φ3);
   const double lon = toDegrees(λ3);
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
   const double R = radius;
   // close polygon so that last point equals first point
   const bool closed = polygon[0] == (polygon[polygon.size() - 1]);
   if (!closed) 
      polygon.push_back(polygon[0]);
   const size_t nVertices = polygon.size() - 1;
   double S = 0; // spherical excess in steradians
   for (int v = 0; v < nVertices; v++) {
      const double φ1 = toRadians(polygon[v].lat());
      const double φ2 = toRadians(polygon[v + 1].lat());
      const double Δλ = toRadians((polygon[v + 1].lon() - polygon[v].lon()));
      const double E = 2 * std::atan2(std::tan(Δλ / 2) * (std::tan(φ1 / 2) + std::tan(φ2 / 2)), 1 + std::tan(φ1 / 2) * std::tan(φ2 / 2));
      S += E;
   }

   // returns whether polygon encloses pole: sum of course deltas around pole is 0° rather than
   // normal ±360°: blog.element84.com/determining-if-a-spherical-polygon-contains-a-pole.html
   auto isPoleEnclosedBy = [](const std::vector<LatLonSpherical>& p)
   {
      //TODO: any better test than this
      double ΣΔ = 0;
      double prevBrng = p[0].initialBearingTo(p[1]);
      for (size_t v = 0; v < p.size() - 1; v++)
      {
         const double initBrng = p[v].initialBearingTo(p[v + 1]);
         const double finalBrng = p[v].finalBearingTo(p[v + 1]);
         ΣΔ += std::fmod((initBrng - prevBrng + 540), 360) - 180;
         ΣΔ += std::fmod((finalBrng - initBrng + 540), 360) - 180;
         prevBrng = finalBrng;
      }
      const double initBrng = p[0].initialBearingTo(p[1]);
      ΣΔ += std::fmod((initBrng - prevBrng + 540), 360) - 180;
      // TODO: fix (intermittant) edge crossing pole - eg (85,90), (85,0), (85,-90)
      return (std::abs(ΣΔ) < 90);
   };

   if (isPoleEnclosedBy(polygon)) {
      S = std::abs(S) - 2 * π;
   }
   const double A = std::abs(S * R * R); // area in units of R
   if (!closed) 
      polygon.pop_back(); // restore polygon to pristine condition
   return A;
}

std::wstring LatLonSpherical::toString(Dms::eFormat e) const
{
   // note: explicitly set dp to undefined for passing through to toLat/toLon
   const std::wstring lat = Dms::toLatitude(m_lat, e);
   const std::wstring lon = Dms::toLongitude(m_lon, e);
   return lat + L"," + lon;
}

std::wstring LatLonSpherical::toGeoJSON() const
{
   return std::wstring(L"{ type: \"Point\", coordinates : [")
      + std::to_wstring(m_lon) + std::wstring(L",")
      + std::to_wstring(m_lat) + std::wstring(L"] }");
}
