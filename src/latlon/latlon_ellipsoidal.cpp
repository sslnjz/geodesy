//
// Created by Bin on 2021/11/15.
//

#include "latlon_ellipsoidal.h"

#include <sstream>
#include <iomanip>

using namespace geodesy;

LatLonEllipsoidal::LatLonEllipsoidal(double lat, double lon, double height)
    : m_lat(lat)
    , m_lon(lon)
    , m_height(height)
    , m_epoch(0.0)
    , m_datum(nullptr)
    , m_referenceFrame(nullptr)
{
}

void LatLonEllipsoidal::setHeight(double height)
{
    m_height = height;
}

void LatLonEllipsoidal::setDatum(const Datum &datum)
{
    m_datum = new Datum{ datum.ellipsoid, datum.transforms };
}

Cartesian LatLonEllipsoidal::toCartesian() const
{
   // x = (ν+h)⋅cosφ⋅cosλ, y = (ν+h)⋅cosφ⋅sinλ, z = (ν⋅(1-e²)+h)⋅sinφ
   // where ν = a/√(1−e²⋅sinφ⋅sinφ), e² = (a²-b²)/a² or (better conditioned) 2⋅f-f²
   const Ellipsoid ellipsoid = m_datum
                     ? m_datum->ellipsoid
                     : m_referenceFrame ? m_referenceFrame->ellipsoid : ellipsoids().WGS84;

   const double φ = toRadians(m_lat);
   const double λ = toRadians(m_lon);
   const double h = m_height;
   const auto a = ellipsoid.a, f = ellipsoid.f;

   const double sinφ = std::sin(φ), cosφ = std::cos(φ);
   const double sinλ = std::sin(λ), cosλ = std::cos(λ);

   const double eSq = 2*f - f*f;                      // 1st eccentricity squared ≡ (a²-b²)/a²
   const double ν = a / std::sqrt(1 - eSq*sinφ*sinφ); // radius of curvature in prime vertical

   const double x = (ν+h) * cosφ * cosλ;
   const double y = (ν+h) * cosφ * sinλ;
   const double z = (ν*(1-eSq)+h) * sinφ;

   return { x, y, z };
}

std::wstring LatLonEllipsoidal::toString(Dms::eFormat format, int dph) const
{
   std::wstringstream hwss;
   hwss << (m_height >= 0 ? L" +" : L" ");
   hwss << std::fixed << std::setprecision(dph) << m_height << L"m";

   std::wstringstream llwss;
   if (format == Dms::N)
   {
      // signed numeric degrees
      llwss << std::fixed << std::setprecision(4);
      llwss << m_lat << ",";
      llwss << m_lon;
      llwss << hwss.str();
      return llwss.str();
   }

   llwss << Dms::toLatitude(m_lat, format) << ", ";
   llwss << Dms::toLatitude(m_lon, format);
   llwss << hwss.str();

   return llwss.str();
}

LatLonEllipsoidal::~LatLonEllipsoidal()
{
    if(nullptr != m_datum)
    {
        delete m_datum;
        m_datum = nullptr;
    }

    if(nullptr != m_referenceFrame)
    {
        delete m_referenceFrame;
        m_referenceFrame = nullptr;
    }
}

Cartesian::Cartesian()
    : vector3d()
{

}

Cartesian::Cartesian(double x, double y, double z)
    : vector3d(x, y, z)
{

}

LatLonEllipsoidal Cartesian::toLatLonEllipsoidal(Ellipsoid ellipsoid) const
{
   // note ellipsoid is available as a parameter for when toLatLon gets subclassed to
   // Ellipsoidal_Datum / Ellipsoidal_Referenceframe.
   const Cartesian cts = *this;
   const Ellipsoid els = ellipsoid;

   const double e2 = 2 * els.f - els.f * els.f; // 1st eccentricity squared ≡ (a²−b²)/a²
   const double ε2 = e2 / (1 - e2); // 2nd eccentricity squared ≡ (a²−b²)/b²
   const double p = std::sqrt(cts.x() * cts.x() + cts.y() * cts.y()); // distance from minor axis
   const double R = std::sqrt(p * p + cts.z() * cts.z()); // polar radius

   // parametric latitude (Bowring eqn.17, replacing tanβ = z·a / p·b)
   const double tanβ = (els.b * cts.z()) / (els.a * p) * (1 + ε2 * els.b / R);
   const double sinβ = tanβ / std::sqrt(1 + tanβ * tanβ);
   const double cosβ = sinβ / tanβ;

   // geodetic latitude (Bowring eqn.18: tanφ = z+ε²⋅b⋅sin³β / p−e²⋅cos³β)
   const double φ = std::isnan(cosβ)
                       ? 0
                       : std::atan2(cts.z() + ε2 * els.b * sinβ * sinβ * sinβ, p - e2 * els.a * cosβ * cosβ * cosβ);

   // longitude
   const double λ = std::atan2(cts.y(), cts.x());

   // height above ellipsoid (Bowring eqn.7)
   const double sinφ = std::sin(φ), cosφ = std::cos(φ);
   const double ν = els.a / std::sqrt(1 - e2 * sinφ * sinφ); // length of the normal terminated by the minor axis
   const double h = p * cosφ + cts.z() * sinφ - (els.a * els.a / ν);

   return {toDegrees(φ), toDegrees(λ), h};
}

