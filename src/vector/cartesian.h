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

#ifndef CARTESIAN_H
#define CARTESIAN_H

#include "vector3d.h"
#include "latlon_ellipsoidal.h"

namespace geodesy
{
   class Cartesian : public vector3d
   {
   public:
      Cartesian();
      Cartesian(double x, double y, double z);

      /**
       * Converts ‘this’ (geocentric) cartesian (x/y/z) coordinate to (geodetic) latitude/longitude
       * point on specified ellipsoid.
       *
       * Uses Bowring’s (1985) formulation for μm precision in concise form; ‘The accuracy of geodetic
       * latitude and height equations’, B R Bowring, Survey Review vol 28, 218, Oct 1985.
       *
       * @param   {LatLon.ellipsoids} [ellipsoid=WGS84] - Ellipsoid to use when converting point.
       * @returns {LatLon} Latitude/longitude point defined by cartesian coordinates, on given ellipsoid.
       * @throws  {TypeError} Invalid ellipsoid.
       *
       * @example
       *   const auto c = new Cartesian(4027893.924, 307041.993, 4919474.294);
       *   const auto p = c.toLatLon(); // 50.7978°N, 004.3592°E
       */
      [[nodiscard]] LatLonEllipsoidal toLatLon(Ellipsoid ellipsoid = s_ellipsoids.WGS84) const;
   };
}

#endif // CARTESIAN_H
