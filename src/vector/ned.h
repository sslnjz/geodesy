
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

#ifndef NED_H
#define NED_H

#include <string>

namespace geodesy
{
    /**
     * North-east-down (NED), also known as local tangent plane (LTP), is a vector in the local
     * coordinate frame of a body.
     */
    class Ned
    {
    public:
        /**
         * Creates North-East-Down vector.
         *
         * @param {number} north - North component in metres.
         * @param {number} east - East component in metres.
         * @param {number} down - Down component (normal to the surface of the ellipsoid) in metres.
         *
         * @example
         *   const auto delta = new Ned(110569, 111297, 1936); // [N:110569,E:111297,D:1936]
         */
        Ned(double n, double e, double d);

        [[nodiscard]] double north() const { return m_north; }
        [[nodiscard]] double east() const  { return m_east; }
        [[nodiscard]] double down() const  { return m_down; }

        /**
         * Length of NED vector.
         *
         * @returns {number} Length of NED vector in metres.
         */
        double length() const;

        /**
         * Bearing of NED vector.
         *
         * @returns {number} Bearing of NED vector in degrees from north.
         */
        double bearing() const;

        /**
         * Elevation of NED vector.
         *
         * @returns {number} Elevation of NED vector in degrees from horizontal (ie tangent to ellipsoid surface).
         */
        double elevation() const;

        /**
         * Creates North-East-Down vector from distance, bearing, & elevation (in local coordinate system).
         *
         * @param   {number} dist - Length of NED vector in metres.
         * @param   {number} brng - Bearing (in degrees from north) of NED vector .
         * @param   {number} elev - Elevation (in degrees from local coordinate frame horizontal) of NED vector.
         * @returns {Ned} North-East-Down vector equivalent to distance, bearing, elevation.
         *
         * @example
         *   const delta = Ned.fromDistanceBearingElevation(116809.178, 222.493, -0.5416); // [N:-86127,E:-78901,D:1104]
         */
        [[maybe_unused]] static Ned fromDistanceBearingElevation(double dist, double brng, double elev);

        /**
         * Returns a string representation of ‘this’ NED vector.
         *
         * @param   {number} [dp=0] - Number of decimal places to display.
         * @returns {string} Comma-separated (labelled) n, e, d values.
         */
        std::string toString(int dp = 0) const;


    private:
        double m_north;
        double m_east;
        double m_down;
    };
}

#endif //NED_H
