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
#include "latlon.h"

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
    class LatLonEllipsoidal : public LatLon
    {
    public:
        LatLonEllipsoidal();
       /**
        * Creates a geodetic latitude/longitude point on a (WGS84) ellipsoidal model earth.
        *
        * @param {number} lat - Latitude (in degrees).
        * @param {number} lon - Longitude (in degrees).
        * @param {Datum}  datum - optional datum
        * @param {number} [height=0] - Height above ellipsoid in metres.
        *
        * @example
        *   const auto p = new LatLonEllipsoidal(51.47788, -0.00147, 17);
        */
       LatLonEllipsoidal(double lat, double lon, double height = 0.0,
                         std::optional<Datum> datum = std::nullopt,
                         std::optional<ReferenceFrame> reference = std::nullopt,
                         std::optional<std::string> epoch = std::nullopt);
       ~LatLonEllipsoidal() override = default;

       /**
        * Height in metres above ellipsoid.
        */
       [[nodiscard]] double height() const;
       void setHeight(double height);

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
        * Parses a latitude/longitude point from a variety of formats.
        *
        * Latitude & longitude (in degrees) can be supplied as two separate parameters, as a single
        * comma-separated lat/lon string, or as a single object with { lat, lon } or GeoJSON properties.
        *
        * The latitude/longitude values may be numeric or strings; they may be signed decimal or
        * deg-min-sec (hexagesimal) suffixed by compass direction (NSEW); a variety of separators are
        * accepted. Examples -3.62, '3 37 12W', '3°37′12″W'.
        *
        * Thousands/decimal separators must be comma/dot; use Dms.fromLocale to convert locale-specific
        * thousands/decimal separators.
        *
        * @param   {number|string|Object} lat|latlon - Latitude (in degrees) or comma-separated lat/lon or lat/lon object.
        * @param   {number|string}        [lon]      - Longitude (in degrees).
        * @returns {LatLon} Latitude/longitude point.
        * @throws  {TypeError} Invalid point.
        *
        * @example
        *   const p1 = LatLon.parse(52.205, 0.119);                                    // numeric pair (≡ new LatLon)
        *   const p2 = LatLon.parse('52.205', '0.119');                                // numeric string pair (≡ new LatLon)
        *   const p3 = LatLon.parse('52.205, 0.119');                                  // single string numerics
        *   const p4 = LatLon.parse('52°12′18.0″N', '000°07′08.4″E');                  // DMS pair
        *   const p5 = LatLon.parse('52°12′18.0″N, 000°07′08.4″E');                    // single string DMS
        *   const p6 = LatLon.parse({ lat: 52.205, lon: 0.119 });                      // { lat, lon } object numeric
        *   const p7 = LatLon.parse({ lat: '52°12′18.0″N', lng: '000°07′08.4″E' });    // { lat, lng } object DMS
        *   const p8 = LatLon.parse({ type: 'Point', coordinates: [ 0.119, 52.205] }); // GeoJSON
        */
       static LatLonEllipsoidal parse(double lat, double lon, double height = 0);
       static LatLonEllipsoidal parse(const std::string& dms, double height = 0);
       static LatLonEllipsoidal parse(const std::string& lat, const std::string& lon, double height);
       static LatLonEllipsoidal parse(const std::string& lat, const std::string& lon, std::string height);


       /**
        * Converts ‘this’ point from (geodetic) latitude/longitude coordinates to (geocentric)
        * cartesian (x/y/z) coordinates.
        *
        * @returns {Cartesian} Cartesian point equivalent to lat/lon point, with x, y, z in metres from
        *   earth centre.
        */
       [[nodiscard]] Cartesian toCartesian() const;

       /**
        * Checks if another point is equal to ‘this’ point.
        *
        * @param   {LatLonEllipsoidal} point - Point to be compared against this point.
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
        * Returns a string representation of ‘this’ point, formatted as degrees, degrees+minutes, or
        * degrees+minutes+seconds.
        *
        * @param   {string} [format=d] - Format point as 'd', 'dm', 'dms', or 'n' for signed numeric.
        * @param   {number} [dp=4|2|0] - Number of decimal places to use: default 4 for d, 2 for dm, 0 for dms.
        * @param   {number} [dpHeight=null] - Number of decimal places to use for height; default is no height display.
        * @returns {string} Comma-separated formatted latitude/longitude.
        * @throws  {RangeError} Invalid format.
        *
        * @example
        *   const greenwich = new LatLon(51.47788, -0.00147, 46);
        *   const d = greenwich.toString();                        // 51.4779°N, 000.0015°W
        *   const dms = greenwich.toString('dms', 2);              // 51°28′40″N, 000°00′05″W
        *   const [lat, lon] = greenwich.toString('n').split(','); // 51.4779, -0.0015
        *   const dmsh = greenwich.toString('dms', 0, 0);          // 51°28′40″N, 000°00′06″W +46m
        */
       [[nodiscard]] std::string toString(Dms::eFormat format = Dms::D, std::optional<int> dp = std::nullopt,
                                    std::optional<int> dph = std::nullopt) const;

       /**
        * Checks if another point is equal to ‘this’ point.
        *
        * @param   {LatLonEllipsoidal} point - Point to be compared against this point.
        * @returns {bool} True if points have identical latitude, longitude, height, and datum/referenceFrame.
        *
        * @example
        *   const auto p1 = new LatLon(52.205, 0.119);
        *   const auto p2 = new LatLon(52.205, 0.119);
        *   const auto equal = p1.equals(p2); // true
        */
       bool operator==(const LatLonEllipsoidal& point) const;

    protected:
        std::optional<std::string>     m_epoch;
        std::optional<Datum>           m_datum;
        std::optional<ReferenceFrame>  m_referenceFrame;

     private:
        double m_height;
    };

   
}

#endif //GEODESY_LATLON_ELLIPSOIDAL_H
