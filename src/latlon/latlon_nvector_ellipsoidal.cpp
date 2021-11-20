
/**********************************************************************************
*  MIT License                                                                    *
*                                                                                 *
*  Copyright (c) 2021 Binbin Song <ssln.jzs@gmail.com>                       *
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

#include "latlon_nvector_ellipsoidal.h"
#include "nvector_ellipsoidal.h"
#include "ned.h"

using namespace geodesy;

LatLonNvectorEllipsoidal::LatLonNvectorEllipsoidal(double lat, double lon, double h, Datum datum)
    : LatLonEllipsoidal(lat, lon, h, datum)
{
}

Cartesian LatLonNvectorEllipsoidal::toCartesian()
{
    const Cartesian c = LatLonEllipsoidal::toCartesian();  // c is 'Cartesian'

    // return Cartesian_Nvector to have toNvector() available as method of exported LatLon
    return NvectorCartesian(c.x(), c.y(), c.z());
}

Ned LatLonNvectorEllipsoidal::deltaTo(const LatLonEllipsoidal &point)
{
    // get delta in cartesian frame
    const auto c1 = toCartesian();
    const auto c2 = point.toCartesian();
    const auto δc = c2 - c1;
    // get local (n-vector) coordinate frame
    const auto n1 = toNvector();
    const auto a = vector3d(0, 0, 1); // axis vector pointing to 90°N
    const auto d = n1.negate();           // down (pointing opposite to n-vector)
    const auto e = a.cross(n1).unit();    // east (pointing perpendicular to the plane)
    const auto n = e.cross(d);            // north (by right hand rule)
    // rotation matrix is built from n-vector coordinate frame axes (using row vectors)
    const double r[3][3] = {
            {n.x(), n.y(), n.z()},
            {e.x(), e.y(), e.z()},
            {d.x(), d.y(), d.z()}
    };
    // apply rotation to δc to get delta in n-vector reference frame
    const auto δn = Cartesian(
            r[0][0]*δc.x() + r[0][1]*δc.y() + r[0][2]*δc.z(),
            r[1][0]*δc.x() + r[1][1]*δc.y() + r[1][2]*δc.z(),
            r[2][0]*δc.x() + r[2][1]*δc.y() + r[2][2]*δc.z()
    );
    return Ned(δn.x(), δn.y(), δn.z());
}

LatLonEllipsoidal LatLonNvectorEllipsoidal::destinationPoint(const Ned &delta)
{
    // convert North-East-Down delta to standard x/y/z vector in coordinate frame of n-vector
    const auto δn = vector3d(delta.north(), delta.east(), delta.down());

    // get local (n-vector) coordinate frame
    const auto n1 = toNvector();
    const auto a = vector3d(0, 0, 1); // axis vector pointing to 90°N
    const auto d = n1.negate();           // down (pointing opposite to n-vector)
    const auto e = a.cross(n1).unit();    // east (pointing perpendicular to the plane)
    const auto n = e.cross(d);            // north (by right hand rule)

    // rotation matrix is built from n-vector coordinate frame axes (using column vectors)
    const double r[3][3] = {
            {n.x(), e.x(), d.x()},
            {n.y(), e.y(), d.y()},
            {n.z(), e.z(), d.z()}
    };

    // apply rotation to δn to get delta in cartesian (ECEF) coordinate reference frame
    const auto δc = Cartesian(
            r[0][0]*δn.x() + r[0][1]*δn.y() + r[0][2]*δn.z(),
            r[1][0]*δn.x() + r[1][1]*δn.y() + r[1][2]*δn.z(),
            r[2][0]*δn.x() + r[2][1]*δn.y() + r[2][2]*δn.z()
    );

    // apply (cartesian) delta to c1 to obtain destination point as cartesian coordinate
    const auto c1 = toCartesian();              // convert this LatLon to Cartesian
    const auto v2 = c1 + δc;                     // the plus() gives us a plain vector,..
    const auto c2 = Cartesian(v2.x(), v2.y(), v2.z()); // ... need to convert it to Cartesian to get LatLon

    // return destination cartesian coordinate as latitude/longitude
    return c2.toLatLon();
}

NvectorEllipsoidal LatLonNvectorEllipsoidal::toNvector()
{ // note: replicated in LatLonNvectorSpherical
    const auto φ = toRadians(lat());
    const auto λ = toRadians(lon());

    const auto sinφ = std::sin(φ), cosφ = std::cos(φ);
    const auto sinλ = std::sin(λ), cosλ = std::cos(λ);

    // right-handed vector: x -> 0°E,0°N; y -> 90°E,0°N, z -> 90°N
    const auto x = cosφ * cosλ;
    const auto y = cosφ * sinλ;
    const auto z = sinφ;

    return NvectorEllipsoidal(x, y, z, height(), datum());
}
