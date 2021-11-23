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
#ifndef LATLON_ELLIPSOIDAL_H
#define LATLON_ELLIPSOIDAL_H

#include <string>
#include <cmath>
#include <optional>

#include "ellipsoids.h"
#include "dms.h"

/**
 * A latitude/longitude point defines a geographic location on or above/below the earth’s surface,
 * measured in degrees from the equator & the International Reference Meridian and in metres above
 * the ellipsoid, and based on a given datum.
 *
 * As so much modern geodesy is based on WGS-84 (as used by GPS), this module includes WGS-84
 * ellipsoid parameters, and it has methods for converting geodetic (latitude/longitude) points to/from
 * geocentric cartesian points; the latlon-ellipsoidal-datum and latlon-ellipsoidal-referenceframe
 * modules provide transformation parameters for converting between historical datums and between
 * modern reference frames.
 *
 * This module is used for both trigonometric geodesy (eg latlon-ellipsoidal-vincenty) and n-vector
 * geodesy (eg latlon-nvector-ellipsoidal), and also for UTM/MGRS mapping.
 *
 * @module latlon-ellipsoidal
 */
namespace geodesy
{
    class Cartesian;
    /**
     * Latitude/longitude points on an ellipsoidal model earth, with ellipsoid parameters and methods
     * for converting points to/from cartesian (ECEF) coordinates.
     *
     * This is the core class, which will usually be used via LatLonEllipsoidal_Datum or
     * LatLonEllipsoidal_ReferenceFrame.
     */
    class LatLonEllipsoidal
    {
    public:
        LatLonEllipsoidal();
       /**
        * Creates a geodetic latitude/longitude point on a (WGS84) ellipsoidal model earth.
        *
        * @param  {number} lat - Latitude (in degrees).
        * @param  {number} lon - Longitude (in degrees).
        * @param  {number} [height=0] - Height above ellipsoid in metres.
        *
        * @example
        *   const auto p = new LatLonEllipsoidal(51.47788, -0.00147, 17);
        */
       LatLonEllipsoidal(double lat, double lon, double height = 0.0,
                         std::optional<Datum> datum = std::nullopt,
                         std::optional<ReferenceFrame> reference = std::nullopt,
                         std::optional<float> epoch = std::nullopt);
       virtual ~LatLonEllipsoidal();

       /**
        * Latitude in degrees north from equator (including aliases lat, latitude): can be set as
        * numeric or hexagesimal (deg-min-sec); returned as numeric.
        */
       [[nodiscard]] double lat() const;
       [[nodiscard]] double latitude() const;
       void setLat(double lat);
       void setLatitude(double lat);

       /**
        * Longitude in degrees east from international reference meridian (including aliases lon, lng,
        * longitude): can be set as numeric or hexagesimal (deg-min-sec); returned as numeric.
        */
       [[nodiscard]] double lon() const;
       [[nodiscard]] double lng() const;
       [[nodiscard]] double longitude() const;
       void setLon(const double lon);
       void setLng(const double lon);
       void setLongitude(const double lon);

       /**
        * Height in metres above ellipsoid.
        */
       [[nodiscard]] double height() const;
       void setHeight(const double height);

       /**
        * Datum.
        *
        * Note this is replicated within LatLonEllipsoidal in order that a LatLonEllipsoidal object can
        * be monkey-patched to look like a LatLonEllipsoidal_Datum, for VincentyInverse calculations on
        * different ellipsoids.
        *
        * @private
        */
       [[nodiscard]] Datum datum() const;
       void setDatum(const Datum& datum);


       /**
        * Ellipsoids with their parameters; this module only defines WGS84 parameters a = 6378137, b =
        * 6356752.314245, f = 1/298.257223563.
        *
        * @example
        *   const auto a = LatLon.ellipsoids.WGS84.a; // 6378137
        */
       static Ellipsoids ellipsoids();

       /**
        * Datums; this module only defines WGS84 datum, hence no datum transformations.
        *
        * @example
        *   const a = LatLon.datums.WGS84.ellipsoid.a; // 6377563.396
        */
       static Datums datums();

       /**
        * Converts ‘this’ point from (geodetic) latitude/longitude coordinates to (geocentric)
        * cartesian (x/y/z) coordinates.
        *
        * @returns {Cartesian} Cartesian point equivalent to lat/lon point, with x, y, z in metres from
        *   earth centre.
        */
       [[nodiscard]] Cartesian toCartesian() const;

       /**
        * Returns a string representation of ‘this’ point, formatted as degrees, degrees+minutes, or
        * degrees+minutes+seconds.
        *
        * @param   {string} [format=d] - Format point as 'd', 'dm', 'dms', or 'n' for signed numeric.
        * @param   {number} [dpHeight=null] - Number of decimal places to use for height; default is no height display.
        * @returns {string} Comma-separated formatted latitude/longitude.
        * @throws  {RangeError} Invalid format.
        *
        * @example
        *   const auto greenwich = new LatLon(51.47788, -0.00147, 46);
        *   const auto d = greenwich.toString();                        // 51.4779°N, 000.0015°W
        *   const auto dms = greenwich.toString('dms', 2);              // 51°28′40″N, 000°00′05″W
        *   const auto latlon = greenwich.toString('n').split(',');     // 51.4779, -0.0015
        *   const auto dmsh = greenwich.toString('dms', 0, 0);          // 51°28′40″N, 000°00′06″W +46m
        */
       [[nodiscard]] std::string toString(Dms::eFormat format = Dms::D, std::optional<int> dph = std::nullopt) const;

       /**
        * Checks if another point is equal to ‘this’ point.
        *
        * @param   {LatLon} point - Point to be compared against this point.
        * @returns {bool} True if points have identical latitude, longitude, height, and datum/referenceFrame.
        * @throws  {TypeError} Invalid point.
        *
        * @example
        *   const p1 = new LatLon(52.205, 0.119);
        *   const p2 = new LatLon(52.205, 0.119);
        *   const equal = p1.equals(p2); // true
        */
       [[nodiscard]] bool equals(const LatLonEllipsoidal& point) const;

       /**
        * Checks if another point is equal to ‘this’ point.
        *
        * @param   {LatLon} point - Point to be compared against this point.
        * @returns {bool} True if points have identical latitude, longitude, height, and datum/referenceFrame.
        *
        * @example
        *   const auto p1 = new LatLon(52.205, 0.119);
        *   const auto p2 = new LatLon(52.205, 0.119);
        *   const auto equal = p1.equals(p2); // true
        */
       inline bool operator==(const LatLonEllipsoidal& point) const;

    protected:
        std::optional<float>  m_epoch;
        std::optional<Datum>  m_datum;
        std::optional<ReferenceFrame> m_referenceFrame;

    private:
        double m_lat;
        double m_lon;
        double m_height;
    };

    inline bool LatLonEllipsoidal::operator==(const LatLonEllipsoidal& point) const
    {
       if (std::fabs(m_lat - point.m_lat) > std::numeric_limits<double>::epsilon()) return false;
       if (std::fabs(m_lon - point.m_lon) > std::numeric_limits<double>::epsilon()) return false;
       if (std::fabs(m_height - point.m_height) > std::numeric_limits<double>::epsilon()) return false;
       if (std::fabs(*m_epoch - *point.m_epoch) > std::numeric_limits<double>::epsilon()) return false;
       if (*m_datum != *point.m_datum) return false;
       if (*m_referenceFrame!= *point.m_referenceFrame) return false;

       return true;
    }
}

#endif //GEODESY_LATLON_ELLIPSOIDAL_H
