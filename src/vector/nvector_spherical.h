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
