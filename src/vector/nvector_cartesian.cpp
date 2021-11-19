
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

#include "nvector_cartesian.h"
#include "ellipsoids.h"
#include "latlon_ellipsoidal.h"
#include "nvector_ellipsoidal.h"

using namespace geodesy;

NvectorCartesian::NvectorCartesian()
{

}

NvectorCartesian::NvectorCartesian(double x, double y, double z)
    : Cartesian(x, y, z)
{

}

NvectorEllipsoidal NvectorCartesian::toNvector(Datum datum) {
    const auto [ a, b, f ] = datum.ellipsoid;

    const auto e2 = 2*f - f*f; // e² = 1st eccentricity squared ≡ (a²-b²)/a²
    const auto e4 = e2*e2;     // e⁴

    const auto p = (x()*x() + y()*y()) / (a*a);
    const auto q = z()*z() * (1-e2) / (a*a);
    const auto r = (p + q - e4) / 6;
    const auto s = (e4*p*q) / (4*r*r*r);
    const auto t = std::cbrt(1 + s + std::sqrt(2*s+s*s));
    const auto u = r * (1 + t + 1/t);
    const auto v = std::sqrt(u*u + e4*q);
    const auto w = e2 * (u + v - q) / (2*v);
    const auto k = std::sqrt(u + v + w*w) - w;
    const auto d = k * std::sqrt(x()*x() + y()*y()) / (k + e2);

    const auto tmp = 1 / std::sqrt(d*d + z()*z());
    const auto xʹ = tmp * k/(k+e2) * x();
    const auto yʹ = tmp * k/(k+e2) * y();
    const auto zʹ = tmp * z();
    const auto h = (k + e2 - 1)/k * std::sqrt(d*d + z()*z());

    return NvectorEllipsoidal(xʹ, yʹ, zʹ, h, datum);
}
