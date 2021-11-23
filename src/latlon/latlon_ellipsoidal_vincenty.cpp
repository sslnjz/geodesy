
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

#include "latlon_ellipsoidal_vincenty.h"

using namespace geodesy;

struct geodesy::VI
{
    double distance;
    double initialBearing;
    double finalBearing;
    int    iterations;
};

struct geodesy::VD
{
    LatLonEllipsoidalVincenty point;
    double finalBearing;
    int    iterations;
};

LatLonEllipsoidalVincenty::LatLonEllipsoidalVincenty(double lat, double lon, double height,
    std::optional<Datum> datum,
    std::optional<ReferenceFrame> reference,
    std::optional<float> epoch)
    : LatLonEllipsoidal(lat, lon, height, datum, reference, epoch)
{

}

LatLonEllipsoidalVincenty::LatLonEllipsoidalVincenty()
    : LatLonEllipsoidal()
{

}

double LatLonEllipsoidalVincenty::distanceTo(const LatLonEllipsoidal &point)
{
    try
    {
        const auto dist = inverse(point).distance;
        return dist; // round to 1mm precision
    }
    catch (const std::domain_error& e)
    {
        return NAN;
    }
    catch (const std::exception& e)
    {
        throw e;
    }
}

VI LatLonEllipsoidalVincenty::inverse(const LatLonEllipsoidal &point)
{
    if (height() != 0 || point.height() != 0)
    {
        throw std::range_error("point must be on the surface of the ellipsoid");
    }

    const LatLonEllipsoidal p1 = *this, p2 = point;
    const auto φ1 = toRadians(p1.lat()), λ1 = toRadians(p1.lon());
    const auto φ2 = toRadians(p2.lat()), λ2 = toRadians(p2.lon());

    // allow alternative ellipsoid to be specified
    const auto ellipsoid = m_datum.has_value() ? m_datum.value().ellipsoid : LatLonEllipsoidal::ellipsoids().WGS84;
    const auto [ a, b, f ] = ellipsoid;

    const auto L = λ2 - λ1; // L = difference in longitude, U = reduced latitude, defined by tan U = (1-f)·tanφ.
    const auto tanU1 = (1-f) * std::tan(φ1), cosU1 = 1 / std::sqrt((1 + tanU1*tanU1)), sinU1 = tanU1 * cosU1;
    const auto tanU2 = (1-f) * std::tan(φ2), cosU2 = 1 / std::sqrt((1 + tanU2*tanU2)), sinU2 = tanU2 * cosU2;

    const auto antipodal = std::abs(L) > π/2 ||std::abs(φ2-φ1) > π/2;

    double λ = L, sinλ = 0.0, cosλ = 0.0; // λ = difference in longitude on an auxiliary sphere
    double σ = antipodal ? π : 0, sinσ = 0.0, cosσ = antipodal ? -1.0 : 1.0, sinSqσ = 0.0; // σ = angular distance P₁ P₂ on the sphere
    double cos2σ = 1.0;                      // σ = angular distance on the sphere from the equator to the midpoint of the line
    double sinα = 0.0, cosSqα = 1.0;         // α = azimuth of the geodesic at the equator
    double C = 0.0;

    double λʹ = 0.0;
    int iterations = 0.0;
    do {
        sinλ = std::sin(λ);
        cosλ = std::cos(λ);
        sinSqσ = (cosU2*sinλ) * (cosU2*sinλ) + (cosU1*sinU2-sinU1*cosU2*cosλ) * (cosU1*sinU2-sinU1*cosU2*cosλ);
        if (std::abs(sinSqσ) < ε) break;  // co-incident/antipodal points (falls back on λ/σ = L)
        sinσ = std::sqrt(sinSqσ);
        cosσ = sinU1*sinU2 + cosU1*cosU2*cosλ;
        σ = std::atan2(sinσ, cosσ);
        sinα = cosU1 * cosU2 * sinλ / sinσ;
        cosSqα = 1 - sinα*sinα;
        cos2σ = (cosSqα != 0) ? (cosσ - 2*sinU1*sinU2/cosSqα) : 0; // on equatorial line cos²α = 0 (§6)
        C = f/16*cosSqα*(4+f*(4-3*cosSqα));
        λʹ = λ;
        λ = L + (1-C) * f * sinα * (σ + C*sinσ*(cos2σ+C*cosσ*(-1+2*cos2σ*cos2σ)));
        const auto iterationCheck = antipodal ? std::abs(λ)-π : std::abs(λ);
        if (iterationCheck > π)
        {
            throw std::domain_error("λ > π");
        }
    } while (std::abs(λ-λʹ) > 1e-12 && ++iterations<1000);

    if (iterations >= 1000) {
        throw std::domain_error("VincentyInverse formula failed to converge");
    }

    const double uSq = cosSqα * (a*a - b*b) / (b*b);
    const double A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
    const double B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));
    const double Δσ = B*sinσ*(cos2σ+B/4*(cosσ*(-1+2*cos2σ*cos2σ)-
                                         B/6*cos2σ*(-3+4*sinσ*sinσ)*(-3+4*cos2σ*cos2σ)));

    const auto s = b*A*(σ-Δσ); // s = length of the geodesic

    // note special handling of exactly antipodal points where sin²σ = 0 (due to discontinuity
    // atan2(0, 0) = 0 but atan2(ε, 0) = π/2 / 90°) - in which case bearing is always meridional,
    // due north (or due south!)
    // α = azimuths of the geodesic; α2 the direction P₁ P₂ produced
    const auto α1 = std::abs(sinSqσ) < ε ? 0 : std::atan2(cosU2*sinλ,  cosU1*sinU2-sinU1*cosU2*cosλ);
    const auto α2 = std::abs(sinSqσ) < ε ? π : std::atan2(cosU1*sinλ, -sinU1*cosU2+cosU1*sinU2*cosλ);

    return {
            s,
            std::abs(s) < ε ? NAN : Dms::wrap360(toDegrees(α1)),
            std::abs(s) < ε ? NAN : Dms::wrap360(toDegrees(α2)),
            iterations,
    };
}

LatLonEllipsoidalVincenty LatLonEllipsoidalVincenty::destinationPoint(double distance, double initialBearing)
{
    return direct(distance, initialBearing).point;
}

VD LatLonEllipsoidalVincenty::direct(double distance, double initialBearing)
{
    if (std::isnan(distance)) { throw std::runtime_error("invalid distance"); }
    if (distance == 0) { return {*this, NAN, 0}; }
    if (std::isnan(initialBearing)) { throw std::runtime_error("invalid bearing"); }
    if (height() != 0) { throw std::range_error("point must be on the surface of the ellipsoid"); }

    const auto φ1 = toRadians(lat()), λ1 = toRadians(lon());
    const auto α1 = toRadians(initialBearing);
    const auto s = distance;

    // allow alternative ellipsoid to be specified
    const auto ellipsoid = m_datum.has_value() ? m_datum.value().ellipsoid : LatLonEllipsoidal::ellipsoids().WGS84;
    const auto [ a, b, f ] = ellipsoid;

    const auto sinα1 = std::sin(α1);
    const auto cosα1 = std::cos(α1);

    const auto tanU1 = (1-f) * std::tan(φ1), cosU1 = 1 / std::sqrt((1 + tanU1*tanU1)), sinU1 = tanU1 * cosU1;
    const auto σ1 = std::atan2(tanU1, cosα1); // σ1 = angular distance on the sphere from the equator to P1
    const auto sinα = cosU1 * sinα1;          // α = azimuth of the geodesic at the equator
    const auto cosSqα = 1 - sinα*sinα;
    const auto uSq = cosSqα * (a*a - b*b) / (b*b);
    const auto A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
    const auto B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));

    double σ = s / (b*A), sinσ = 0.0, cosσ = 0.0, Δσ = 0.0; // σ = angular distance P₁ P₂ on the sphere
    double cos2σ = 0.0; // σₘ = angular distance on the sphere from the equator to the midpoint of the line

    double σʹ = 0.0;
    int iterations = 0;
    do {
        cos2σ = std::cos(2*σ1 + σ);
        sinσ  = std::sin(σ);
        cosσ  = std::cos(σ);
        Δσ = B*sinσ*(cos2σ+B/4*(cosσ*(-1+2*cos2σ*cos2σ)-
                                 B/6*cos2σ*(-3+4*sinσ*sinσ)*(-3+4*cos2σ*cos2σ)));
        σʹ = σ;
        σ = s / (b*A) + Δσ;
    } while (std::abs(σ-σʹ) > 1e-12 && ++iterations<100);
    if (iterations >= 100)
        throw std::domain_error("Vincenty formula failed to converge"); // not possible?

    const auto x = sinU1*sinσ - cosU1*cosσ*cosα1;
    const auto φ2 = std::atan2(sinU1*cosσ + cosU1*sinσ*cosα1, (1-f)*std::sqrt(sinα*sinα + x*x));
    const auto λ = std::atan2(sinσ*sinα1, cosU1*cosσ - sinU1*sinσ*cosα1);
    const auto C = f/16*cosSqα*(4+f*(4-3*cosSqα));
    const auto L = λ - (1-C) * f * sinα * (σ + C*sinσ*(cos2σ+C*cosσ*(-1+2*cos2σ*cos2σ)));
    const auto λ2 = λ1 + L;

    const auto α2 = std::atan2(sinα, -x);

    const auto destinationPoint = LatLonEllipsoidalVincenty(toDegrees(φ2), toDegrees(λ2), 0, m_datum);

    return {
        destinationPoint,
        Dms::wrap360(toDegrees(α2)),
        iterations,
    };
}

double LatLonEllipsoidalVincenty::finalBearingTo(const LatLonEllipsoidal &point)
{
    try
    {
        const auto brng = inverse(point).finalBearing;
        return brng; // round to 0.001″ precision
    }
    catch (const std::runtime_error& e)
    {
        return NAN;
    }
    catch (const std::exception& e)
    {
        throw e;
    }
}

double LatLonEllipsoidalVincenty::initialBearingTo(const LatLonEllipsoidal &point)
{
    try
    {
        const auto brng = inverse(point).initialBearing;
        return brng;
    }
    catch (const std::runtime_error& e)
    {
        return NAN;
    }
    catch (const std::exception& e)
    {
        throw e;
    }
}

double LatLonEllipsoidalVincenty::finalBearingOn(double distance, double initialBearing)
{
    const auto brng = direct(distance, initialBearing).finalBearing;
    return brng; // round to 0.001″ precision
}

LatLonEllipsoidalVincenty LatLonEllipsoidalVincenty::intermediatePointTo(
        const LatLonEllipsoidal &point, double fraction)
{
    if (fraction == 0) return *this;
    const auto inverse = this->inverse(point);
    const auto dist = inverse.distance;
    const auto brng = inverse.initialBearing;
    return std::isnan(brng) ? *this : destinationPoint(dist * fraction, brng);
}
