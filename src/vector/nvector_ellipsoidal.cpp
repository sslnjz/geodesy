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

#include "nvector_ellipsoidal.h"

using namespace geodesy;

NvectorEllipsoidal::NvectorEllipsoidal(double x, double y, double z, double h, Datum datum)
        : vector3d(vector3d(x, y, z).unit()), m_h(h), m_datum(datum)
{

}

std::string NvectorEllipsoidal::toString(int dp, std::optional<int> dph)
{
    std::stringstream ssh;
    if(m_h >= 0)
    {
        ssh << "+" << std::setprecision(dph.value_or(0)) << m_h << "m";
    }
    std::stringstream ss;
    ss << "[" << std::setprecision(dp) << x() << ", " << y() << ", " << z() << "]"
        << (bool(dph) ? ssh.str() : "");
    return ss.str();
}

NvectorCartesian NvectorEllipsoidal::toCartesian()
{
    const auto [ a, b, f ] = m_datum.ellipsoid;

    const auto m = (1-f) * (1-f); // (1−f)² = b²/a²
    const auto n = b / std::sqrt(x()*x()/m + y()*y()/m + z()*z());

    const auto xʹ = n * x() / m + x() * m_h;
    const auto yʹ = n * y() / m + y() * m_h;
    const auto zʹ = n * z()     + z() * m_h;

    return NvectorCartesian(xʹ, yʹ, zʹ);
}

LatLonNvectorEllipsoidal NvectorEllipsoidal::toLatLon()
{
    // tanφ = z / √(x²+y²), tanλ = y / x (same as spherical calculation)

    const auto phi = std::atan2(z(), std::sqrt(x()*x() + y()*y()));
    const auto lambda = std::atan2(y(), x());
    return LatLonNvectorEllipsoidal(toDegrees(phi), toDegrees(lambda), m_h, m_datum);
}
