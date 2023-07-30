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
#ifndef CARTESIAN_REFERENCEFRAME_H
#define CARTESIAN_REFERENCEFRAME_H

#include "geodesy/cartesian.h"

namespace geodesy {
   /**
    * Augments Cartesian with reference frame and observation epoch the cooordinate is based on, and
    * methods to convert between reference frames (using Helmert 14-parameter transforms) and to
    * convert cartesian to geodetic latitude/longitude point.
    *
    * @extends Cartesian
    */
   class LatLonEllipsoidalReferenceFrame;
   class CartesianReferenceFrame : public Cartesian
   {
   public:
      CartesianReferenceFrame();
      /**
       * Creates cartesian coordinate representing ECEF (earth-centric earth-fixed) point, on a given
       * reference frame. The reference frame will identify the primary meridian (for the x-coordinate),
       * and is also useful in transforming to/from geodetic (lat/lon) coordinates.
       *
       * @param  {number} x - X coordinate in metres (=> 0°N,0°E).
       * @param  {number} y - Y coordinate in metres (=> 0°N,90°E).
       * @param  {number} z - Z coordinate in metres (=> 90°N).
       * @param  {LatLon.referenceFrames} [referenceFrame] - Reference frame this coordinate is defined within.
       * @param  {number} [epoch=referenceFrame.epoch] - date of observation of coordinate (decimal year).
       * @throws {TypeError} Unrecognised reference frame, Invalid epoch.
       *
       * @example
       *   import { Cartesian } from '/js/geodesy/latlon-ellipsoidal-referenceframe.js';
       *   const coord = new Cartesian(3980581.210, -111.159, 4966824.522);
       */
      CartesianReferenceFrame(double x, double y, double z,
                              std::optional<ReferenceFrame> referenceFrame = std::nullopt,
                              std::optional<std::string> epoch = std::nullopt);

      /**
       * Reference frame this point is defined within.
       */
      [[nodiscard]] std::optional<ReferenceFrame> referenceFrame() const;
      void setReferenceFrame(const ReferenceFrame& referenceFrame);

      /**
       * Point’s observed epoch.
       */
      [[nodiscard]] std::optional<std::string> epoch();
      void setEpoch(std::string epoch);

      /**
       * Converts ‘this’ (geocentric) cartesian (x/y/z) coordinate to (geodetic) latitude/longitude
       * point (based on the same reference frame).
       *
       * Shadow of Cartesian.toLatLon(), returning LatLon augmented with LatLonEllipsoidal_ReferenceFrame
       * methods convertReferenceFrame, toCartesian, etc.
       *
       * @returns {LatLon} Latitude/longitude point defined by cartesian coordinates, in given reference frame.
       * @throws  {Error} No reference frame defined.
       *
       * @example
       *   const c = new Cartesian(4027893.924, 307041.993, 4919474.294, LatLon.referenceFrames.ITRF2000);
       *   const p = c.toLatLon(); // 50.7978°N, 004.3592°E
       */
      [[nodiscard]] LatLonEllipsoidalReferenceFrame toLatLon() const;

      /**
       * Converts ‘this’ cartesian coordinate to new reference frame using Helmert 14-parameter
       * transformation. The observation epoch is unchanged.
       *
       * Note that different conversions have different tolerences; refer to the literature if
       * tolerances are significant.
       *
       * @param   {LatLon.referenceFrames} toReferenceFrame - Reference frame this coordinate is to be converted to.
       * @returns {Cartesian} This point converted to new reference frame.
       * @throws  {Error} Undefined reference frame.
       *
       * @example
       *   const c = new Cartesian(3980574.247, -102.127, 4966830.065, LatLon.referenceFrames.ITRF2000);
       *   c.convertReferenceFrame(LatLon.referenceFrames.ETRF2000); // [3980574.395,-102.214,4966829.941](ETRF2000@1997.0)
       */
      [[nodiscard]] CartesianReferenceFrame convertReferenceFrame(const ReferenceFrame& to) const;

      /**
       * Applies Helmert 14-parameter transformation to ‘this’ coordinate using supplied transform
       * parameters and annual rates of change, with the secular variation given by the difference
       * between the reference epoch t0 and the observation epoch tc.
       *
       * This is used in converting reference frames.
       *
       * See e.g. 3D Coordinate Transformations, Deakin, 1998.
       *
       * @private
       * @param   {number[]} params - Transform parameters tx, ty, tz, s, rx, ry, rz..
       * @param   {number[]} rates - Rate of change of transform parameters ṫx, ṫy, ṫz, ṡ, ṙx, ṙy, ṙz.
       * @param   {number}   δt - Period between reference and observed epochs, t − t₀.
       * @returns {Cartesian} Transformed point (without reference frame).
       */
      [[nodiscard]] CartesianReferenceFrame applyTransform(const double params[7], const double rates[7], double deltat) const;


      /**
       * Returns a string representation of ‘this’ cartesian point. TRF is shown if set, and
       * observation epoch if different from reference epoch.
       *
       * @param   {number} [dp=0] - Number of decimal places to use.
       * @returns {string} Comma-separated latitude/longitude.
       */
      [[nodiscard]] std::string toString(int dp) const override;

   private:
      std::optional<ReferenceFrame> _referenceFrame;
      std::optional<std::string> _epoch;
   };
}

#endif //CARTESIAN_REFERENCEFRAME_H
