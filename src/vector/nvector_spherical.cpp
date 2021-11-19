#include "nvector_spherical.h"
#include "latlon_nvector_spherical.h"

using namespace geodesy;

NvectorSpherical::NvectorSpherical(double x, double y, double z)
   : vector3d( vector3d(x, y, z).unit())
{

}

vector3d NvectorSpherical::greatCircle(double bearing) const
{
   const double θ = toRadians(bearing);

   const vector3d N = vector3d(0, 0, 1); // n-vector representing north pole
   const vector3d e = N.cross(*this); // easting
   const vector3d n = this->cross(e); // northing
   const vector3d eʹ = e * (std::cos(θ) / e.length());
   const vector3d nʹ = n * (std::sin(θ) / n.length());
   const vector3d c = nʹ - (eʹ);

   return c;
}

LatLonNvectorSpherical NvectorSpherical::toLatLon() const
{
   // tanφ = z / √(x²+y²), tanλ = y / x (same as ellipsoidal calculation)
   const double φ = std::atan2(z(), std::sqrt(x() * x() + y() * y()));
   const double λ = std::atan2(y(), x());

   return { toDegrees(φ), toDegrees(λ) };
}
