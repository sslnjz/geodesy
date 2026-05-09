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

#include "latlon_osgridref.h"

#include "algorithm.h"
#include "osgridref.h"

#include <cmath>

namespace geodesy
{

LatLonOsGridRef::LatLonOsGridRef(double lat, double lon, double height, std::optional<Datum> datum)
      : LatLonEllipsoidalDatum(lat, lon, height, datum)
{
}

OsGridRef LatLonOsGridRef::toOsGrid() const
{
   // The National Grid projection is defined on OSGB36, so WGS84 inputs first pass through Helmert conversion.
   const LatLonOsGridRef point = this->datum() == datums().OSGB36
                         ? *this
                         : convertDatum(datums().OSGB36);

   const double phi = toRadians(point.lat());
   const double lambda = toRadians(point.lon());

   const auto [ a, b, f ] = nationalGrid.ellipsoid;
   const double phi0 = toRadians(nationalGrid.trueOrigin.lat);
   const double lambda0 = toRadians(nationalGrid.trueOrigin.lon);
   const double E0 = -nationalGrid.falseOrigin.easting;
   const double N0 = -nationalGrid.falseOrigin.northing;
   const double F0 = nationalGrid.scaleFactor;

   const double e2 = 1.0 - (b * b) / (a * a);
   const double n = (a - b) / (a + b);
   const double n2 = n * n;
   const double n3 = n2 * n;

   const double cosPhi = std::cos(phi);
   const double sinPhi = std::sin(phi);
   const double nu = a * F0 / std::sqrt(1.0 - e2 * sinPhi * sinPhi);
   const double rho = a * F0 * (1.0 - e2) / std::pow(1.0 - e2 * sinPhi * sinPhi, 1.5);
   const double eta2 = nu / rho - 1.0;

   // Meridional arc from true origin to latitude on the Airy 1830 ellipsoid.
   const double Ma = (1.0 + n + (5.0 / 4.0) * n2 + (5.0 / 4.0) * n3) * (phi - phi0);
   const double Mb = (3.0 * n + 3.0 * n2 + (21.0 / 8.0) * n3)
      * std::sin(phi - phi0) * std::cos(phi + phi0);
   const double Mc = ((15.0 / 8.0) * n2 + (15.0 / 8.0) * n3)
      * std::sin(2.0 * (phi - phi0)) * std::cos(2.0 * (phi + phi0));
   const double Md = (35.0 / 24.0) * n3
      * std::sin(3.0 * (phi - phi0)) * std::cos(3.0 * (phi + phi0));
   const double M = b * F0 * (Ma - Mb + Mc - Md);

   const double cos3Phi = cosPhi * cosPhi * cosPhi;
   const double cos5Phi = cos3Phi * cosPhi * cosPhi;
   const double tanPhi = std::tan(phi);
   const double tan2Phi = tanPhi * tanPhi;
   const double tan4Phi = tan2Phi * tan2Phi;

   const double I = M + N0;
   const double II = (nu / 2.0) * sinPhi * cosPhi;
   const double III = (nu / 24.0) * sinPhi * cos3Phi * (5.0 - tan2Phi + 9.0 * eta2);
   const double IIIA = (nu / 720.0) * sinPhi * cos5Phi * (61.0 - 58.0 * tan2Phi + tan4Phi);
   const double IV = nu * cosPhi;
   const double V = (nu / 6.0) * cos3Phi * (nu / rho - tan2Phi);
   const double VI = (nu / 120.0) * cos5Phi
      * (5.0 - 18.0 * tan2Phi + tan4Phi + 14.0 * eta2 - 58.0 * tan2Phi * eta2);

   const double deltaLambda = lambda - lambda0;
   const double deltaLambda2 = deltaLambda * deltaLambda;
   const double deltaLambda3 = deltaLambda2 * deltaLambda;
   const double deltaLambda4 = deltaLambda3 * deltaLambda;
   const double deltaLambda5 = deltaLambda4 * deltaLambda;
   const double deltaLambda6 = deltaLambda5 * deltaLambda;

   double northing = I + II * deltaLambda2 + III * deltaLambda4 + IIIA * deltaLambda6;
   double easting = E0 + IV * deltaLambda + V * deltaLambda3 + VI * deltaLambda5;

   // Keep the reference implementation's millimetre rounding before constructing the value type.
   northing = std::stod(toFixed(northing, 3));
   easting = std::stod(toFixed(easting, 3));

   return OsGridRef(easting, northing);
}

LatLonOsGridRef LatLonOsGridRef::convertDatum(const Datum& toDatum) const
{
   const LatLonEllipsoidalDatum osgbED = LatLonEllipsoidalDatum::convertDatum(toDatum);
   return LatLonOsGridRef(osgbED.lat(), osgbED.lon(), osgbED.height(), osgbED.datum());
}

}
