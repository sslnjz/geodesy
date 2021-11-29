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

LatLonEllipsoidal Cartesian::toLatLon(Ellipsoid ellipsoid) const
{
   // note ellipsoid is available as a parameter for when toLatLon gets subclassed to
   // Ellipsoidal_Datum / Ellipsoidal_Referenceframe.
   //const Cartesian cts = *this;
   //const Ellipsoid els = ellipsoid;
   auto [a, b, f] = ellipsoid;

   const double e2 = 2 * f - f * f; // 1st eccentricity squared ≡ (a²−b²)/a²
   const double epsilon2 = e2 / (1 - e2); // 2nd eccentricity squared ≡ (a²−b²)/b²
   const double p = std::sqrt(x() * x() + y() * y()); // distance from minor axis
   const double R = std::sqrt(p * p + z() * z()); // polar radius

   // parametric latitude (Bowring eqn.17, replacing tanβ = z·a / p·b)
   const double tanbeta = (b * z()) / (a * p) * (1 + epsilon2 * b / R);
   const double sinbeta = tanbeta / std::sqrt(1 + tanbeta * tanbeta);
   const double cosbeta = sinbeta / tanbeta;

   // geodetic latitude (Bowring eqn.18: tanφ = z+ε²⋅b⋅sin³β / p−e²⋅cos³β)
   const double phi = std::isnan(cosbeta)
      ? 0
      : std::atan2(z() + epsilon2 * b * sinbeta * sinbeta * sinbeta, p - e2 * a * cosbeta * cosbeta * cosbeta);

   // longitude
   const double lambda = std::atan2(y(), x());

   // height above ellipsoid (Bowring eqn.7)
   const double sinphi = std::sin(phi), cosphi = std::cos(phi);
   const double nu = a / std::sqrt(1 - e2 * sinphi * sinphi); // length of the normal terminated by the minor axis
   const double h = p * cosphi + z() * sinphi - (a * a / nu);

   return { toDegrees(phi), toDegrees(lambda), h };
}