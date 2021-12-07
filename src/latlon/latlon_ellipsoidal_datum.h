/**********************************************************************************
*  MIT License                                                                    *
*                                                                                 *
*  Copyright (c) 2021 Binbin Song <ssln.jzs@gmail.com>                         *
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
#ifndef LATLON_ELLIPSOIDAL_DATUM_H
#define LATLON_ELLIPSOIDAL_DATUM_H

#include "latlon_ellipsoidal.h"

namespace geodesy
{
   class CartesianDatum;
   class LatLonEllipsoidalDatum : public LatLonEllipsoidal
   {
   public:
      /**
       * Creates a geodetic latitude/longitude point on an ellipsoidal model earth using given datum.
       *
       * @param {number} lat - Latitude (in degrees).
       * @param {number} lon - Longitude (in degrees).
       * @param {number} [height=0] - Height above ellipsoid in metres.
       * @param {LatLon.datums} datum - Datum this point is defined within.
       *
       * @example
       *   import LatLon from '/js/geodesy/latlon-ellipsoidal-datum.js';
       *   const p = new LatLon(53.3444, -6.2577, 17, LatLon.datums.Irl1975);
       */
      LatLonEllipsoidalDatum( double lat, double lon, double height = 0, std::optional<Datum> datum= g_datums.WGS84);

      /**
       * Datum this point is defined within.
       */
      std::optional<Datum> datum()
      {
         return m_datum;
      }

      /**
       * Ellipsoids with their parameters; semi-major axis (a), semi-minor axis (b), and flattening (f).
       *
       * Flattening f = (a−b)/a; at least one of these parameters is derived from defining constants.
       *
       * @example
       *   const a = LatLon.ellipsoids.Airy1830.a; // 6377563.396
       */
      static Ellipsoids ellipsoids()
      {
         return g_ellipsoids;
      }


      /**
       * Datums; with associated ellipsoid, and Helmert transform parameters to convert from WGS-84
       * into given datum.
       *
       * Note that precision of various datums will vary, and WGS-84 (original) is not defined to be
       * accurate to better than ±1 metre. No transformation should be assumed to be accurate to
       * better than a metre, for many datums somewhat less.
       *
       * This is a small sample of commoner datums from a large set of historical datums. I will add
       * new datums on request.
       *
       * @example
       *   const a = LatLon.datums.OSGB36.ellipsoid.a;                    // 6377563.396
       *   const tx = LatLon.datums.OSGB36.transform;                     // [ tx, ty, tz, s, rx, ry, rz ]
       *   const availableDatums = Object.keys(LatLon.datums).join(', '); // ED50, Irl1975, NAD27, ...
       */
      static Datums datums()
      {
         return g_datums;
      }

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
       * @param   {number|string|Object} lat|latlon - Geodetic Latitude (in degrees) or comma-separated lat/lon or lat/lon object.
       * @param   {number}               [lon] - Longitude in degrees.
       * @param   {number}               [height=0] - Height above ellipsoid in metres.
       * @param   {LatLon.datums}        [datum=WGS84] - Datum this point is defined within.
       * @returns {LatLon} Latitude/longitude point on ellipsoidal model earth using given datum.
       * @throws  {TypeError} Unrecognised datum.
       *
       * @example
       *   const p = LatLon.parse('51.47736, 0.0000', 0, LatLon.datums.OSGB36);
       */
      static LatLonEllipsoidalDatum parse(double lat, double lon, double height = 0, Datum datum = g_datums.WGS84);
      static LatLonEllipsoidalDatum parse(const std::string& dms, double height = 0, Datum datum = g_datums.WGS84);
      static LatLonEllipsoidalDatum parse(const std::string& lat, const std::string& lon, double height, Datum datum = g_datums.WGS84);
      static LatLonEllipsoidalDatum parse(const std::string& lat, const std::string& lon, std::string height, Datum datum = g_datums.WGS84);

      /**
       * Converts ‘this’ lat/lon coordinate to new coordinate system.
       *
       * @param   {LatLon.datums} toDatum - Datum this coordinate is to be converted to.
       * @returns {LatLon} This point converted to new datum.
       * @throws  {TypeError} Unrecognised datum.
       *
       * @example
       *   const pWGS84 = new LatLon(51.47788, -0.00147, 0, LatLon.datums.WGS84);
       *   const pOSGB = pWGS84.convertDatum(LatLon.datums.OSGB36); // 51.4773°N, 000.0001°E
       */
      LatLonEllipsoidalDatum convertDatum(Datum toDatum) const;


      /**
       * Converts ‘this’ point from (geodetic) latitude/longitude coordinates to (geocentric) cartesian
       * (x/y/z) coordinates, based on the same datum.
       *
       * Shadow of LatLonEllipsoidal.toCartesian(), returning Cartesian augmented with
       * LatLonEllipsoidal_Datum methods/properties.
       *
       * @returns {Cartesian} Cartesian point equivalent to lat/lon point, with x, y, z in metres from
       *   earth centre, augmented with reference frame conversion methods and properties.
       */
      CartesianDatum toCartesian() const;

   };
}

#endif //LATLON_ELLIPSOIDAL_DATUM_H
