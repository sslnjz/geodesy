
#ifndef NVECTOR_SPHERICAL_H
#define NVECTOR_SPHERICAL_H

#include "vector3d.h"

namespace geodesy
{
   /**
    * An n-vector is a (unit) vector normal to the Earth's surface (a non-singular position
    * representation).
    *
    * For many applications, n-vectors are more convenient to work with than other position
    * representations such as latitude/longitude, UTM coordinates, etc.
    *
    * On a spherical model earth, an n-vector is equivalent to a (normalised) earth-centred earth-fixed
    * (ECEF) vector.
    *
    * @extends vector3d
    */

   class LatLonNvectorSpherical;
   class NvectorSpherical : public vector3d
   {
   public:
      /**
       * Creates a 3d n-vector normal to the Earth’s surface.
       *
       * @param {number} x - X component of n-vector (towards 0°N, 0°E).
       * @param {number} y - Y component of n-vector (towards 0°N, 90°E).
       * @param {number} z - Z component of n-vector (towards 90°N).
       *
       * @example
       *   const auto n = new Nvector(0.5000, 0.5000, 0.7071);
       */
      NvectorSpherical(double x, double y, double z);

      /**
       * Vector normal to great circle obtained by heading on given bearing from point given by
       * ‘this’ n-vector.
       *
       * Direction of vector is such that initial bearing vector b = c × n, where n is an n-vector
       * representing ‘this’ (start) point.
       *
       * @private
       * @param   {number}   bearing - Compass bearing in degrees.
       * @returns {Vector3d} Normalised vector representing great circle.
       *
       * @example
       *   const auto n1 = new LatLon(53.3206, -1.7297).toNvector();
       *   const auto gc = n1.greatCircle(96.0); // [-0.794,0.129,0.594]
       */
      [[nodiscard]] vector3d greatCircle(double bearing) const;


      /**
       * Converts ‘this’ n-vector to latitude/longitude point.
       *
       * @returns  {LatLon} Latitude/longitude point vector points to.
       *
       * @example
       *   const n = new Nvector(0.5000, 0.5000, 0.7071);
       *   const p = n.toLatLon(); // 45.0°N, 045.0°E
       */
      [[nodiscard]] LatLonNvectorSpherical toLatLon() const;
   };
}


#endif // NVECTOR_SPHERICAL_H
