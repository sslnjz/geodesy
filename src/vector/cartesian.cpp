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
#include "cartesian.h"

#include <cmath>
#include <stdexcept>

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
   auto [a, b, f] = ellipsoid;

   if (!std::isfinite(x()) || !std::isfinite(y()) || !std::isfinite(z()))
   {
      throw std::invalid_argument("invalid cartesian coordinate");
   }

   if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(f)
      || a <= 0 || b <= 0 || b > a || f <= 0 || f >= 1)
   {
      throw std::invalid_argument("invalid ellipsoid");
   }

   const double e2 = 2 * f - f * f; // 1st eccentricity squared ≡ (a²−b²)/a²
   const double epsilon2 = e2 / (1 - e2); // 2nd eccentricity squared ≡ (a²−b²)/b²
   const double p = std::sqrt(x() * x() + y() * y()); // distance from minor axis
   const double R = std::sqrt(p * p + z() * z()); // polar radius

   if (R == 0.0)
   {
      throw std::domain_error("cartesian origin has no geodetic latitude or longitude");
   }

   if (p == 0.0)
   {
      // Longitude is undefined on the minor axis; the reference convention is 0°.
      // Height is measured along the ellipsoid normal, so polar height is |z|-b.
      const double lat = z() < 0 ? -90.0 : 90.0;
      const double h = std::fabs(z()) - b;
      return { lat, 0.0, h };
   }

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
