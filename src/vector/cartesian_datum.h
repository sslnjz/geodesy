/**********************************************************************************
*  MIT License                                                                    *
*                                                                                 *
*  Copyright (c) 2021 Binbin Song <ssln.jzs@gmail.com>                         *
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
#ifndef CARTESIAN_DATUM_H
#define CARTESIAN_DATUM_H

#include "cartesian.h"

namespace geodesy
{
   class LatLonEllipsoidalDatum;
   class CartesianDatum : public Cartesian
   {
   public:
      /**
        * Creates cartesian coordinate representing ECEF (earth-centric earth-fixed) point, on a given
        * datum. The datum will identify the primary meridian (for the x-coordinate), and is also
        * useful in transforming to/from geodetic (lat/lon) coordinates.
        *
        * @param  {number} x - X coordinate in metres (=> 0°N,0°E).
        * @param  {number} y - Y coordinate in metres (=> 0°N,90°E).
        * @param  {number} z - Z coordinate in metres (=> 90°N).
        * @param  {LatLon.datums} [datum] - Datum this coordinate is defined within.
        * @throws {TypeError} Unrecognised datum.
        *
        * @example
        *   import { Cartesian } from '/js/geodesy/latlon-ellipsoidal-datum.js';
        *   const coord = new Cartesian(3980581.210, -111.159, 4966824.522);
        */
      CartesianDatum(double x, double y, double z, const Datum& datum = g_datums.WGS84);

      /**
       * Datum this point is defined within.
       */
      Datum datum();
      void setDatum(const Datum& datum);

      /**
        * Converts ‘this’ (geocentric) cartesian (x/y/z) coordinate to (geodetic) latitude/longitude
        * point (based on the same datum, or WGS84 if unset).
        *
        * Shadow of Cartesian.toLatLon(), returning LatLon augmented with LatLonEllipsoidal_Datum
        * methods convertDatum, toCartesian, etc.
        *
        * @returns {LatLon} Latitude/longitude point defined by cartesian coordinates.
        * @throws  {TypeError} Unrecognised datum
        *
        * @example
        *   const c = new Cartesian(4027893.924, 307041.993, 4919474.294);
        *   const p = c.toLatLon(); // 50.7978°N, 004.3592°E
        */
       LatLonEllipsoidalDatum toLatLon(std::optional<Datum> deprecatedDatum = std::nullopt);


      /**
       * Converts ‘this’ cartesian coordinate to new datum using Helmert 7-parameter transformation.
       *
       * @param   {LatLon.datums} toDatum - Datum this coordinate is to be converted to.
       * @returns {Cartesian} This point converted to new datum.
       * @throws  {Error} Undefined datum.
       *
       * @example
       *   const c = new Cartesian(3980574.247, -102.127, 4966830.065, LatLon.datums.OSGB36);
       *   c.convertDatum(LatLon.datums.Irl1975); // [??,??,??]
       */
      CartesianDatum convertDatum(Datum toDatum);


      /**
       * Applies Helmert 7-parameter transformation to ‘this’ coordinate using transform parameters t.
       *
       * This is used in converting datums (geodetic->cartesian, apply transform, cartesian->geodetic).
       *
       * @private
       * @param   {number[]} t - Transformation to apply to this coordinate.
       * @returns {Cartesian} Transformed point.
       */
      CartesianDatum applyTransform(Transform t) const;

      /**
       * String representation of vector.
       *
       * @param   {number} [dp=3] - Number of decimal places to be used.
       * @returns {string} Vector represented as [x,y,z].
       */
      [[nodiscard]] std::string toString(int dp = 0) const override;

   private:
      Datum m_datum;
   };
}

#endif //CARTESIAN_DATUM_H
