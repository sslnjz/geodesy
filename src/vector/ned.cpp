
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

#include "ned.h"

#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>

#include "dms.h"
#include "algorithm.h"

using namespace geodesy;

Ned::Ned(double n, double e, double d)
    : m_north(n), m_east(e), m_down(d)
{

}

double Ned::length()
{
    return std::sqrt(m_north*m_north + m_east*m_east + m_down*m_down);
}

double Ned::bearing()
{
    const auto θ = std::atan2(m_east, m_north);
    return Dms::wrap360(toDegrees(θ)); // normalise to range 0..360°
}

double Ned::elevation()
{
    const auto α = std::asin(m_down/length());
    return toDegrees(-α);
}

Ned Ned::fromDistanceBearingElevation(double dist, double brng, double elev)
{
    const auto θ = toRadians(brng);
    const auto α = toRadians(elev);

    const auto sinθ = std::sin(θ), cosθ = std::cos(θ);
    const auto sinα = std::sin(α), cosα = std::cos(α);
    const auto n = cosθ * dist*cosα;
    const auto e = sinθ * dist*cosα;
    const auto d = -sinα * dist;
    return Ned(n, e, d);
}

std::string Ned::toString(int dp)
{
    std::stringstream ss;
    ss << std::setprecision(dp) << "[N:"
       << m_north << ",E:"
       << m_east << ",D:"
       << m_down << "]";
    return ss.str();
}
