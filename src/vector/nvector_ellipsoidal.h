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

#ifndef NVECTOR_ELLIPSOIDAL_H
#define NVECTOR_ELLIPSOIDAL_H

#include <optional>

#include "geodesy/vector3d.h"
#include "geodesy/ellipsoids.h"
#include "geodesy/latlon_ellipsoidal.h"
#include "geodesy/latlon_nvector_ellipsoidal.h"

namespace geodesy
{
    /**
     * An n-vector is a position representation using a (unit) vector normal to the Earth ellipsoid.
     * Unlike latitude/longitude points, n-vectors have no singularities or discontinuities.
     *
     * For many applications, n-vectors are more convenient to work with than other position
     * representations such as latitude/longitude, earth-centred earth-fixed (ECEF) vectors, UTM
     * coordinates, etc.
     *
     * @extends Vector3d
     */
    class NvectorEllipsoidal : public vector3d
    {
    public:
        // note commonality with latlon-nvector-spherical

        /**
         * Creates a 3d n-vector normal to the Earth's surface.
         *
         * @param {number} x - X component of n-vector (towards 0°N, 0°E).
         * @param {number} y - Y component of n-vector (towards 0°N, 90°E).
         * @param {number} z - Z component of n-vector (towards 90°N).
         * @param {number} [h=0] - Height above ellipsoid surface in metres.
         * @param {LatLon.datums} [datum=WGS84] - Datum this n-vector is defined within.
         */
        NvectorEllipsoidal(double x, double y, double z, double h = 0.0,
                           Datum datum = LatLonEllipsoidal::datums().WGS84);


        /**
         * Converts ‘this’ n-vector to latitude/longitude point.
         *
         * @returns {LatLon} Latitude/longitude point equivalent to this n-vector.
         *
         * @example
         *   const p = new Nvector(0.500000, 0.500000, 0.707107).toLatLon(); // 45.0000°N, 045.0000°E
         */
        LatLonNvectorEllipsoidal toLatLon();


        /**
         * Converts ‘this’ n-vector to cartesian coordinate.
         *
         * qv Gade 2010 ‘A Non-singular Horizontal Position Representation’ eqn 22
         *
         * @returns {Cartesian} Cartesian coordinate equivalent to this n-vector.
         *
         * @example
         *   const c = new Nvector(0.500000, 0.500000, 0.707107).toCartesian(); // [3194419,3194419,4487349]
         *   const p = c.toLatLon();                                            // 45.0000°N, 045.0000°E
         */
        NvectorCartesian toCartesian();


        /**
         * Returns a string representation of ‘this’ (unit) n-vector. Height component is only shown if
         * dpHeight is specified.
         *
         * @param   {number} [dp=3] - Number of decimal places to display.
         * @param   {number} [dph=null] - Number of decimal places to use for height; default is no height display.
         * @returns {string} Comma-separated x, y, z, h values.
         *
         * @example
         *   new Nvector(0.5000, 0.5000, 0.7071).toString();        // [0.500,0.500,0.707]
         *   new Nvector(0.5000, 0.5000, 0.7071, 1).toString(6, 0); // [0.500002,0.500002,0.707103+1m]
         */
        std::string toString(int dp =3, std::optional<int> dph = std::nullopt);

    private:
        double m_h;
        Datum m_datum;
    };

}

#endif //GEODESY_NVECTOR_ELLIPSOIDAL_H
