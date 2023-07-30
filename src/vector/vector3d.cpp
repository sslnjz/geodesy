#include "vector3d.h"

using geodesy::vector3d;

/**
* Creates a 3-d vector.
*
* @param {number} x - X component of vector.
* @param {number} y - Y component of vector.
* @param {number} z - Z component of vector.
*
* @example
*   v = new vector3d(0.267, 0.535, 0.802);
*/

vector3d::vector3d() : vx(0), vy(0), vz(0) { }

vector3d::vector3d(double x, double y, double z) : vx(x), vy(y), vz(z)
{
    if (std::isnan(x) || std::isnan(y) || std::isnan(z)) {
        std::stringstream ss;
        ss << "invalid vector[" << x << ", " << y << ", " << z << "]";
        throw std::runtime_error(ss.str());
    }
}
