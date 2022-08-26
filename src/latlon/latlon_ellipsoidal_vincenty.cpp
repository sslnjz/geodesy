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
    std::optional<std::string> epoch)
    : LatLonEllipsoidal(lat, lon, height, datum, reference, epoch)
{

}

LatLonEllipsoidalVincenty::LatLonEllipsoidalVincenty()
    : LatLonEllipsoidal()
{

}

double LatLonEllipsoidalVincenty::distanceTo(const LatLonEllipsoidal &point) const
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

VI LatLonEllipsoidalVincenty::inverse(const LatLonEllipsoidal &point) const
{
    if (height() != 0 || point.height() != 0)
    {
        throw std::range_error("point must be on the surface of the ellipsoid");
    }

    const LatLonEllipsoidal p1 = *this, p2 = point;
    const auto phi1 = toRadians(p1.lat()), lambda1 = toRadians(p1.lon());
    const auto phi2 = toRadians(p2.lat()), lambda2 = toRadians(p2.lon());

    // allow alternative ellipsoid to be specified
    const auto ellipsoid = m_datum.has_value() ? m_datum.value().ellipsoid : LatLonEllipsoidal::ellipsoids().WGS84;
    const auto [ a, b, f ] = ellipsoid;

    const auto L = lambda2 - lambda1; // L = difference in longitude, U = reduced latitude, defined by tan U = (1-f)·tanφ.
    const auto tanU1 = (1-f) * std::tan(phi1), cosU1 = 1 / std::sqrt((1 + tanU1*tanU1)), sinU1 = tanU1 * cosU1;
    const auto tanU2 = (1-f) * std::tan(phi2), cosU2 = 1 / std::sqrt((1 + tanU2*tanU2)), sinU2 = tanU2 * cosU2;

    const auto antipodal = std::abs(L) > pi/2 ||std::abs(phi2-phi1) > pi/2;

    double lambda = L, sinlambda = 0.0, coslambda = 0.0; // λ = difference in longitude on an auxiliary sphere
    double sigma = antipodal ? pi : 0, sinsigma = 0.0, cossigma = antipodal ? -1.0 : 1.0, sinSqsigma = 0.0; // σ = angular distance P₁ P₂ on the sphere
    double cos2sigma = 1.0;                      // σ = angular distance on the sphere from the equator to the midpoint of the line
    double sinalpha = 0.0, cosSqalpha = 1.0;         // α = azimuth of the geodesic at the equator
    double C = 0.0;

    double lambda_d = 0.0;
    int iterations = 0.0;
    do {
        sinlambda = std::sin(lambda);
        coslambda = std::cos(lambda);
        sinSqsigma = (cosU2*sinlambda) * (cosU2*sinlambda) + (cosU1*sinU2-sinU1*cosU2*coslambda) * (cosU1*sinU2-sinU1*cosU2*coslambda);
        if (std::abs(sinSqsigma) < epsilon) break;  // co-incident/antipodal points (falls back on λ/σ = L)
        sinsigma = std::sqrt(sinSqsigma);
        cossigma = sinU1*sinU2 + cosU1*cosU2*coslambda;
        sigma = std::atan2(sinsigma, cossigma);
        sinalpha = cosU1 * cosU2 * sinlambda / sinsigma;
        cosSqalpha = 1 - sinalpha*sinalpha;
        cos2sigma = (cosSqalpha != 0) ? (cossigma - 2*sinU1*sinU2/cosSqalpha) : 0; // on equatorial line cos²α = 0 (§6)
        C = f/16*cosSqalpha*(4+f*(4-3*cosSqalpha));
        lambda_d = lambda;
        lambda = L + (1-C) * f * sinalpha * (sigma + C*sinsigma*(cos2sigma+C*cossigma*(-1+2*cos2sigma*cos2sigma)));
        const auto iterationCheck = antipodal ? std::abs(lambda)-pi : std::abs(lambda);
        if (iterationCheck > pi)
        {
            throw std::domain_error("λ > π");
        }
    } while (std::abs(lambda-lambda_d) > 1e-12 && ++iterations<1000);

    if (iterations >= 1000) {
        throw std::domain_error("VincentyInverse formula failed to converge");
    }

    const double uSq = cosSqalpha * (a*a - b*b) / (b*b);
    const double A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
    const double B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));
    const double DELTAsigma = B*sinsigma*(cos2sigma+B/4*(cossigma*(-1+2*cos2sigma*cos2sigma)-
                                         B/6*cos2sigma*(-3+4*sinsigma*sinsigma)*(-3+4*cos2sigma*cos2sigma)));

    const auto s = b*A*(sigma-DELTAsigma); // s = length of the geodesic

    // note special handling of exactly antipodal points where sin²σ = 0 (due to discontinuity
    // atan2(0, 0) = 0 but atan2(epsilon, 0) = π/2 / 90°) - in which case bearing is always meridional,
    // due north (or due south!)
    // α = azimuths of the geodesic; α2 the direction P₁ P₂ produced
    const auto alpha1 = std::abs(sinSqsigma) < epsilon ? 0 : std::atan2(cosU2*sinlambda,  cosU1*sinU2-sinU1*cosU2*coslambda);
    const auto alpha2 = std::abs(sinSqsigma) < epsilon ? pi : std::atan2(cosU1*sinlambda, -sinU1*cosU2+cosU1*sinU2*coslambda);

    return {
            s,
            std::abs(s) < epsilon ? NAN : Dms::wrap360(toDegrees(alpha1)),
            std::abs(s) < epsilon ? NAN : Dms::wrap360(toDegrees(alpha2)),
            iterations,
    };
}

LatLonEllipsoidalVincenty LatLonEllipsoidalVincenty::destinationPoint(double distance, double initialBearing) const
{
    return direct(distance, initialBearing).point;
}

VD LatLonEllipsoidalVincenty::direct(double distance, double initialBearing) const
{
    if (std::isnan(distance)) { throw std::runtime_error("invalid distance"); }
    if (distance == 0) { return {*this, NAN, 0}; }
    if (std::isnan(initialBearing)) { throw std::runtime_error("invalid bearing"); }
    if (height() != 0) { throw std::range_error("point must be on the surface of the ellipsoid"); }

    const auto phi1 = toRadians(lat()), lambda1 = toRadians(lon());
    const auto alpha1 = toRadians(initialBearing);
    const auto s = distance;

    // allow alternative ellipsoid to be specified
    const auto ellipsoid = m_datum.has_value() ? m_datum.value().ellipsoid : LatLonEllipsoidal::ellipsoids().WGS84;
    const auto [ a, b, f ] = ellipsoid;

    const auto sinalpha1 = std::sin(alpha1);
    const auto cosalpha1 = std::cos(alpha1);

    const auto tanU1 = (1-f) * std::tan(phi1), cosU1 = 1 / std::sqrt((1 + tanU1*tanU1)), sinU1 = tanU1 * cosU1;
    const auto sigma1 = std::atan2(tanU1, cosalpha1); // σ1 = angular distance on the sphere from the equator to P1
    const auto sinalpha = cosU1 * sinalpha1;          // α = azimuth of the geodesic at the equator
    const auto cosSqalpha = 1 - sinalpha*sinalpha;
    const auto uSq = cosSqalpha * (a*a - b*b) / (b*b);
    const auto A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
    const auto B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));

    double sigma = s / (b*A), sinsigma = 0.0, cossigma = 0.0, DELTAsigma = 0.0; // σ = angular distance P₁ P₂ on the sphere
    double cos2sigma = 0.0; // σₘ = angular distance on the sphere from the equator to the midpoint of the line

    double sigma_d = 0.0;
    int iterations = 0;
    do {
        cos2sigma = std::cos(2*sigma1 + sigma);
        sinsigma  = std::sin(sigma);
        cossigma  = std::cos(sigma);
        DELTAsigma = B*sinsigma*(cos2sigma+B/4*(cossigma*(-1+2*cos2sigma*cos2sigma)-
                                 B/6*cos2sigma*(-3+4*sinsigma*sinsigma)*(-3+4*cos2sigma*cos2sigma)));
        sigma_d = sigma;
        sigma = s / (b*A) + DELTAsigma;
    } while (std::abs(sigma-sigma_d) > 1e-12 && ++iterations<100);
    if (iterations >= 100)
        throw std::domain_error("Vincenty formula failed to converge"); // not possible?

    const auto x = sinU1*sinsigma - cosU1*cossigma*cosalpha1;
    const auto phi2 = std::atan2(sinU1*cossigma + cosU1*sinsigma*cosalpha1, (1-f)*std::sqrt(sinalpha*sinalpha + x*x));
    const auto lambda = std::atan2(sinsigma*sinalpha1, cosU1*cossigma - sinU1*sinsigma*cosalpha1);
    const auto C = f/16*cosSqalpha*(4+f*(4-3*cosSqalpha));
    const auto L = lambda - (1-C) * f * sinalpha * (sigma + C*sinsigma*(cos2sigma+C*cossigma*(-1+2*cos2sigma*cos2sigma)));
    const auto lambda2 = lambda1 + L;

    const auto alpha2 = std::atan2(sinalpha, -x);

    const auto destinationPoint = LatLonEllipsoidalVincenty(toDegrees(phi2), toDegrees(lambda2), 0, m_datum);

    return {
        destinationPoint,
        Dms::wrap360(toDegrees(alpha2)),
        iterations,
    };
}

double LatLonEllipsoidalVincenty::finalBearingTo(const LatLonEllipsoidal &point) const
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

double LatLonEllipsoidalVincenty::initialBearingTo(const LatLonEllipsoidal &point) const
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

double LatLonEllipsoidalVincenty::finalBearingOn(double distance, double initialBearing) const
{
    const auto brng = direct(distance, initialBearing).finalBearing;
    return brng; // round to 0.001″ precision
}

LatLonEllipsoidalVincenty LatLonEllipsoidalVincenty::intermediatePointTo(
        const LatLonEllipsoidal &point, double fraction) const
{
    if (fraction == 0) return *this;
    const auto inverse = this->inverse(point);
    const auto dist = inverse.distance;
    const auto brng = inverse.initialBearing;
    return std::isnan(brng) ? *this : destinationPoint(dist * fraction, brng);
}
