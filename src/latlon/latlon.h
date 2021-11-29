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
#ifndef LATLON_H
#define LATLON_H

#include <string>

#include "dms.h"

namespace geodesy
{
    /**
     * Base LatLon class
     */
    class LatLon
    {
    public:
        LatLon(double lat, double lon);
        virtual ~LatLon() = default;

        [[nodiscard]] double lat() const;
        [[nodiscard]] double latitude() const;

        void setLat(double lat);
        void setLatitude(double lat);

        /**
         * Longitude in degrees east from international reference meridian (including aliases lon, lng,
         * longitude): can be set as numeric or hexadecimal (deg-min-sec); returned as numeric.
         */
        [[nodiscard]] double lon() const;
        [[nodiscard]] double lng() const;
        [[nodiscard]] double longitude() const;

        void setLon(double lon);
        void setLng(double lon);
        void setLongitude(double lon);

        /**
         * Returns a string representation of ‘this’ point, formatted as degrees, degrees+minutes, or
         * degrees+minutes+seconds.
         *
         * @param   {enum} [eFormat=D] - Format point as 'D', 'DM', 'DMS', or 'N' for signed numeric.
         * @returns {string} Comma-separated formatted latitude/longitude.
         *
         * @example
         *   const greenwich = new LatLon(51.47788, -0.00147);
         *   const d = greenwich.toString();                        // 51.4779°N, 000.0015°W
         *   const dms = greenwich.toString('dms', 2);              // 51°28′40.37″N, 000°00′05.29″W
         *   const [lat, lon] = greenwich.toString('n').split(','); // 51.4779, -0.0015
         */
        [[nodiscard]] std::string toString(Dms::eFormat e = Dms::D, std::optional<int> dp = std::nullopt) const;

        /**
         * Converts ‘this’ point to a GeoJSON object string.
         *
         * @returns {string} this point as a GeoJSON ‘Point’ string.
         *    { type: "Point", coordinates: [ lon, lat ] }
         */
        [[nodiscard]] std::string toGeoJSON() const;

        /**
       * Checks if another point is equal to ‘this’ point.
       *
       * @param   {LatLon} p1 - Point to be compared against p2.
       * @param   {LatLon} p2 - Point to be compared against p1.
       * @returns {bool}   True if points have identical latitude and longitude values.
       *
       */
        friend bool operator==(const LatLon& p1, const LatLon& p2);

        /**
         * Checks if another point is NOT equal to ‘this’ point.
         *
         * @param   {LatLon} p1 - Point to be compared against p2.
         * @param   {LatLon} p2 - Point to be compared against p1.
         * @returns {bool}   True if points have different latitude or longitude values.
         *
         */
        friend bool operator!=(const LatLon& p1, const LatLon& p2);

    protected:
        double m_lat; // Latitude in degrees north from equator
        double m_lon; // Longitude in degrees east from international reference meridian
    };
}

#endif //LATLON_H
