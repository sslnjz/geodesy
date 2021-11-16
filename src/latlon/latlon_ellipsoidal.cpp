//
// Created by Bin on 2021/11/15.
//

#include "latlon_ellipsoidal.h"

#include <sstream>
#include <iomanip>

using namespace geodesy;

LatLonEllipsoidal::LatLonEllipsoidal(double lat, double lon, double height)
    : _lat(lat)
    , _lon(lon)
    , _height(height)
{
}

void LatLonEllipsoidal::setHeight(double height)
{
    if (std::isnan(_height)) throw std::invalid_argument("invalid height");
    _height = height;
}

std::wstring LatLonEllipsoidal::toString(int dp)
{
    return L"";
}

void LatLonEllipsoidal::setDatum(const Datum &datum)
{
    _datum = datum;
}

Cartesian LatLonEllipsoidal::toCartesian()
{
//            // x = (ν+h)⋅cosφ⋅cosλ, y = (ν+h)⋅cosφ⋅sinλ, z = (ν⋅(1-e²)+h)⋅sinφ
//            // where ν = a/√(1−e²⋅sinφ⋅sinφ), e² = (a²-b²)/a² or (better conditioned) 2⋅f-f²
//            const Ellipsoid ellipsoid = _datum
//                              ? _datum.ellipsoid
//                              : this.referenceFrame ? this.referenceFrame.ellipsoid : ellipsoids.WGS84;
//
//            const φ = this.lat.toRadians();
//            const λ = this.lon.toRadians();
//            const h = this.height;
//            const { a, f } = ellipsoid;
//
//            const sinφ = Math.sin(φ), cosφ = Math.cos(φ);
//            const sinλ = Math.sin(λ), cosλ = Math.cos(λ);
//
//            const eSq = 2*f - f*f;                      // 1st eccentricity squared ≡ (a²-b²)/a²
//            const ν = a / Math.sqrt(1 - eSq*sinφ*sinφ); // radius of curvature in prime vertical
//
//            const x = (ν+h) * cosφ * cosλ;
//            const y = (ν+h) * cosφ * sinλ;
//            const z = (ν*(1-eSq)+h) * sinφ;
//
//            return new Cartesian(x, y, z);
    return {0, 0, 0};
}

Cartesian::Cartesian()
    : vector3d()
{

}

Cartesian::Cartesian(double x, double y, double z)
    : vector3d(x, y, z)
{

}


