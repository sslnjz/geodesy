
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

#ifndef LATLON_NVECTOR_ELLIPSOIDAL_H
#define LATLON_NVECTOR_ELLIPSOIDAL_H

#include "latlon_ellipsoidal.h"
#include "cartesian.h"
#include "nvector_cartesian.h"

#include <optional>

/**
 * Tools for working with points on (ellipsoidal models of) the earth’s surface using a vector-based
 * approach using ‘n-vectors’ (rather than the more common spherical trigonometry).
 *
 * Based on Kenneth Gade’s ‘Non-singular Horizontal Position Representation’.
 *
 * Note that these formulations take x => 0°N,0°E, y => 0°N,90°E, z => 90°N (in order that n-vector
 * = cartesian vector at 0°N,0°E); Gade uses x => 90°N, y => 0°N,90°E, z => 0°N,0°E.
 *
 * @module latlon-nvector-ellipsoidal
 */

namespace geodesy
{
    class LatLonNvectorEllipsoidal : public LatLonEllipsoidal
    {
    public:
        LatLonNvectorEllipsoidal(double lat, double lon, double h, Datum datum);
        /**
         * Converts ‘this’ point from (geodetic) latitude/longitude coordinates to (geocentric) cartesian
         * (x/y/z) coordinates.
         *
         * @returns {Cartesian} Cartesian point equivalent to lat/lon point, with x, y, z in metres from
         *   earth centre.
         */
        Cartesian toCartesian();
    };
}



#endif //LATLON_NVECTOR_ELLIPSOIDAL_H
