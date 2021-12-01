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

#include "cartesian.h"

namespace geodesy {
   /**
    * Augments Cartesian with reference frame and observation epoch the cooordinate is based on, and
    * methods to convert between reference frames (using Helmert 14-parameter transforms) and to
    * convert cartesian to geodetic latitude/longitude point.
    *
    * @extends Cartesian
    */
   class LatLonEllipsoidalReferenceFrame;
   class CartesianReferenceFrame : public Cartesian {
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
      std::optional<ReferenceFrame> referenceFrame() {
         return _referenceFrame;
      }

      void setReferenceFrame(const ReferenceFrame &referenceFrame) {
         if (!referenceFrame.epoch.has_value())
            throw std::invalid_argument("unrecognised reference frame");
         _referenceFrame = referenceFrame;
      }

      /**
       * Point’s observed epoch.
       */
      std::optional<std::string> epoch() {
         return _epoch ? *_epoch : _referenceFrame ? _referenceFrame->epoch : std::nullopt;
      }

      void setEpoch(std::string epoch) {
         if (_referenceFrame.has_value() &&
             _referenceFrame->epoch.has_value() &&
             _epoch != (*_referenceFrame).epoch.value()) {
            _epoch = epoch;
         }
      }

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
      LatLonEllipsoidalReferenceFrame toLatLon();


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
      CartesianReferenceFrame convertReferenceFrame(const ReferenceFrame& to);
//
//
//      /**
//       * Applies Helmert 14-parameter transformation to ‘this’ coordinate using supplied transform
//       * parameters and annual rates of change, with the secular variation given by the difference
//       * between the reference epoch t0 and the observation epoch tc.
//       *
//       * This is used in converting reference frames.
//       *
//       * See e.g. 3D Coordinate Transformations, Deakin, 1998.
//       *
//       * @private
//       * @param   {number[]} params - Transform parameters tx, ty, tz, s, rx, ry, rz..
//       * @param   {number[]} rates - Rate of change of transform parameters ṫx, ṫy, ṫz, ṡ, ṙx, ṙy, ṙz.
//       * @param   {number}   δt - Period between reference and observed epochs, t − t₀.
//       * @returns {Cartesian} Transformed point (without reference frame).
//       */
//      applyTransform(params, rates, δt
//      )   {
//         // this point
//         const x1 = this.x, y1 = this.y, z1 = this.z;
//
//         // base parameters
//         const tx = params[0] / 1000;                    // x-shift: normalise millimetres to metres
//         const ty = params[1] / 1000;                    // y-shift: normalise millimetres to metres
//         const tz = params[2] / 1000;                    // z-shift: normalise millimetres to metres
//         const s = params[3] / 1e9;                     // scale: normalise parts-per-billion
//         const rx = (params[4] / 3600 / 1000).toRadians(); // x-rotation: normalise milliarcseconds to radians
//         const ry = (params[5] / 3600 / 1000).toRadians(); // y-rotation: normalise milliarcseconds to radians
//         const rz = (params[6] / 3600 / 1000).toRadians(); // z-rotation: normalise milliarcseconds to radians
//
//         // rate parameters
//         const ṫx = rates[0] / 1000;                     // x-shift: normalise millimetres to metres
//         const ṫy = rates[1] / 1000;                     // y-shift: normalise millimetres to metres
//         const ṫz = rates[2] / 1000;                     // z-shift: normalise millimetres to metres
//         const ṡ = rates[3] / 1e9;                      // scale: normalise parts-per-billion
//         const ṙx = (rates[4] / 3600 / 1000).toRadians();  // x-rotation: normalise milliarcseconds to radians
//         const ṙy = (rates[5] / 3600 / 1000).toRadians();  // y-rotation: normalise milliarcseconds to radians
//         const ṙz = (rates[6] / 3600 / 1000).toRadians();  // z-rotation: normalise milliarcseconds to radians
//
//         // combined (normalised) parameters
//         const T = {x: tx + ṫx * δt, y: ty + ṫy * δt, z: tz + ṫz * δt};
//         const R = {x: rx + ṙx * δt, y: ry + ṙy * δt, z: rz + ṙz * δt};
//         const S = 1 + s + ṡ * δt;
//
//         // apply transform (shift, scale, rotate)
//         const x2 = T.x + x1 * S - y1 * R.z + z1 * R.y;
//         const y2 = T.y + x1 * R.z + y1 * S - z1 * R.x;
//         const z2 = T.z - x1 * R.y + y1 * R.x + z1 * S;
//
//         return new Cartesian_ReferenceFrame(x2, y2, z2);
//      }
//
//
//      /**
//       * Returns a string representation of ‘this’ cartesian point. TRF is shown if set, and
//       * observation epoch if different from reference epoch.
//       *
//       * @param   {number} [dp=0] - Number of decimal places to use.
//       * @returns {string} Comma-separated latitude/longitude.
//       */
//      toString(dp = 0
//      ) {
//         const { x, y, z } = this;
//         const epochFmt = {useGrouping: false, minimumFractionDigits: 1, maximumFractionDigits: 20};
//         const epoch = this.referenceFrame && this.epoch != this.referenceFrame.epoch ? this.epoch.toLocaleString('en',
//                                                                                                                  epochFmt)
//                                                                                      : '';
//         const trf = this.referenceFrame ? `(${this.referenceFrame.name}
//         ${epoch ? '@' + epoch : ''})` : '';
//         return `[${x.toFixed(dp)}, ${y.toFixed(dp)}, ${z.toFixed(dp)}]
//         ${trf}`;
//      }

   private:
      std::optional<ReferenceFrame> _referenceFrame;
      std::optional<std::string> _epoch;
   };
}

#endif //CARTESIAN_REFERENCEFRAME_H
