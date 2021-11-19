#include "cartesian.h"

using namespace geodesy;

Cartesian::Cartesian()
   : vector3d()
{

}

Cartesian::Cartesian(double x, double y, double z)
   : vector3d(x, y, z)
{

}

LatLonEllipsoidal Cartesian::toLatLonEllipsoidal(Ellipsoid ellipsoid) const
{
   // note ellipsoid is available as a parameter for when toLatLon gets subclassed to
   // Ellipsoidal_Datum / Ellipsoidal_Referenceframe.
   //const Cartesian cts = *this;
   //const Ellipsoid els = ellipsoid;
   auto [a, b, f] = ellipsoid;

   const double e2 = 2 * f - f * f; // 1st eccentricity squared ≡ (a²−b²)/a²
   const double ε2 = e2 / (1 - e2); // 2nd eccentricity squared ≡ (a²−b²)/b²
   const double p = std::sqrt(x() * x() + y() * y()); // distance from minor axis
   const double R = std::sqrt(p * p + z() * z()); // polar radius

   // parametric latitude (Bowring eqn.17, replacing tanβ = z·a / p·b)
   const double tanβ = (b * z()) / (a * p) * (1 + ε2 * b / R);
   const double sinβ = tanβ / std::sqrt(1 + tanβ * tanβ);
   const double cosβ = sinβ / tanβ;

   // geodetic latitude (Bowring eqn.18: tanφ = z+ε²⋅b⋅sin³β / p−e²⋅cos³β)
   const double φ = std::isnan(cosβ)
      ? 0
      : std::atan2(z() + ε2 * b * sinβ * sinβ * sinβ, p - e2 * a * cosβ * cosβ * cosβ);

   // longitude
   const double λ = std::atan2(y(), x());

   // height above ellipsoid (Bowring eqn.7)
   const double sinφ = std::sin(φ), cosφ = std::cos(φ);
   const double ν = a / std::sqrt(1 - e2 * sinφ * sinφ); // length of the normal terminated by the minor axis
   const double h = p * cosφ + z() * sinφ - (a * a / ν);

   return { toDegrees(φ), toDegrees(λ), h };
}