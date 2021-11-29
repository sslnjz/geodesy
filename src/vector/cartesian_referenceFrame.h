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

namespace geodesy
{
    /**
     * Augments Cartesian with reference frame and observation epoch the cooordinate is based on, and
     * methods to convert between reference frames (using Helmert 14-parameter transforms) and to
     * convert cartesian to geodetic latitude/longitude point.
     *
     * @extends Cartesian
     */
    class CartesianReferenceFrame : public Cartesian
    {
    public:
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
        std::optional<ReferenceFrame> referenceFrame()
        {
            return _referenceFrame;
        }
        void setReferenceFrame(const ReferenceFrame& referenceFrame)
        {
            if (!referenceFrame.epoch.has_value())
                throw std::invalid_argument("unrecognised reference frame");
            _referenceFrame = referenceFrame;
        }

        /**
         * Point’s observed epoch.
         */
        std::optional<std::string> epoch()
        {
            return _epoch ? *_epoch : _referenceFrame ? _referenceFrame->epoch : std::nullopt;
        }
        void setEpoch(std::string epoch)
        {
            if (_referenceFrame.has_value() &&
                _referenceFrame->epoch.has_value() &&
                _epoch != (*_referenceFrame).epoch.value())
            {
                _epoch = epoch;
            }
        }

    private:
        std::optional<ReferenceFrame> _referenceFrame;
        std::optional<std::string> _epoch;
    };
}

#endif //CARTESIAN_REFERENCEFRAME_H
