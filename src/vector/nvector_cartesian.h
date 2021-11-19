
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

#ifndef CARTESIAN_NVECTOR_H
#define CARTESIAN_NVECTOR_H

#include "cartesian.h"

namespace geodesy
{
    class NvectorEllipsoidal;
    class NvectorCartesian : public Cartesian
    {
    public:
        NvectorCartesian();
        NvectorCartesian(double x, double y, double z);
        /**
         * Converts ‘this’ cartesian coordinate to an n-vector.
         *
         * qv Gade 2010 ‘A Non-singular Horizontal Position Representation’ eqn 23
         *
         * @param   {LatLon.datums} [datum=WGS84] - Datum to use for conversion.
         * @returns {Nvector} N-vector equivalent to this cartesian coordinate.
         *
         * @example
         *   const c = new Cartesian(3980581, 97, 4966825);
         *   const n = c.toNvector(); // { x: 0.6228, y: 0.0000, z: 0.7824, h: 0.0000 }
         */
        NvectorEllipsoidal toNvector(Datum datum = LatLonEllipsoidal::datums().WGS84);
    };
}

#endif //GEODESY_CARTESIAN_NVECTOR_H
