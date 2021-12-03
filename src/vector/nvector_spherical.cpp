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
#include "nvector_spherical.h"
#include "latlon_nvector_spherical.h"

using namespace geodesy;

NvectorSpherical::NvectorSpherical(double x, double y, double z)
   : vector3d( vector3d(x, y, z).unit())
{

}

vector3d NvectorSpherical::greatCircle(double bearing) const
{
   const double theta = toRadians(bearing);

   const vector3d N = vector3d(0, 0, 1); // n-vector representing north pole
   const vector3d e = N.cross(*this); // easting
   const vector3d n = this->cross(e); // northing
   const vector3d eʹ = e * (std::cos(theta) / e.length());
   const vector3d nʹ = n * (std::sin(theta) / n.length());
   const vector3d c = nʹ - (eʹ);

   return c;
}

LatLonNvectorSpherical NvectorSpherical::toLatLon() const
{
   // tanφ = z / √(x²+y²), tanλ = y / x (same as ellipsoidal calculation)
   const double phi = std::atan2(z(), std::sqrt(x() * x() + y() * y()));
   const double lambda = std::atan2(y(), x());

   return { toDegrees(phi), toDegrees(lambda) };
}
