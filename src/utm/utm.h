/**********************************************************************************
*  MIT License                                                                    *
*                                                                                 *
*  Copyright (c) 2021 Binbin Song <ssln.jzs@gmail.com>                            *
*                                                                                 *
*  Geodesy tools for conversions between (historical) datums                      *
*  (c) Chris Veness 2005-2019                                                     *
*  www.movable-type.co.uk/scripts/latlong-utm-mgrs.html                           *
*  www.movable-type.co.uk/scripts/geodesy-library.html#utm                        *
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

#ifndef UTM_H
#define UTM_H

/**
 * The Universal Transverse Mercator (UTM) system is a 2-dimensional Cartesian coordinate system
 * providing locations on the surface of the Earth.
 *
 * UTM is a set of 60 transverse Mercator projections, normally based on the WGS-84 ellipsoid.
 * Within each zone, coordinates are represented as eastings and northings, measures in metres; e.g.
 * ‘31 N 448251 5411932’.
 *
 * This method based on Karney 2011 ‘Transverse Mercator with an accuracy of a few nanometers’,
 * building on Krüger 1912 ‘Konforme Abbildung des Erdellipsoids in der Ebene’.
 *
 */

#include <optional>
#include <string>

#include "ellipsoids.h"
#include "latlon_ellipsoidal.h"

namespace geodesy
{
   /**
    * UTM coordinates, with functions to parse them and convert them to LatLon points.
    */
   class LatLonUtm;
   class Utm
   {
   public:
      /**
       * Creates a Utm coordinate object comprising zone, hemisphere, easting, northing on a given
       * datum (normally WGS84).
       *
       * @param  {number}        zone - UTM 6° longitudinal zone (1..60 covering 180°W..180°E).
       * @param  {string}        hemisphere - N for northern hemisphere, S for southern hemisphere.
       * @param  {number}        easting - Easting in metres from false easting (-500km from central meridian).
       * @param  {number}        northing - Northing in metres from equator (N) or from false northing -10,000km (S).
       * @param  {LatLon.datums} [datum=WGS84] - Datum UTM coordinate is based on.
       * @param  {number}        [convergence=null] - Meridian convergence (bearing of grid north
       *                         clockwise from true north), in degrees.
       * @param  {number}        [scale=null] - Grid scale factor.
       * @param {boolean=true}  verifyEN - Check easting/northing is within 'normal' values (may be
       *                         suppressed for extended coherent coordinates or alternative datums
       *                         e.g. ED50 (epsg.io/23029).
       * @throws {TypeError} Invalid UTM coordinate.
       *
       * @example
       *   import Utm from '/js/geodesy/utm.js';
       *   const utmCoord = new Utm(31, 'N', 448251, 5411932);
       */

      enum class Hemisphere { N, S };


      Utm(int zone, Hemisphere h, double easting, double northing,
          std::optional<Datum> datum = LatLonEllipsoidal::datums().WGS84,
          std::optional<double> convergence = std::nullopt, std::optional<double> scale = std::nullopt,
          bool verifyEN = true);


      /**
       * Parses string representation of UTM coordinate.
       *
       * A UTM coordinate comprises (space-separated)
       *  - zone
       *  - hemisphere
       *  - easting
       *  - northing.
       *
       * @param   {string} utmCoord - UTM coordinate (WGS 84).
       * @param   {Datum}  [datum=WGS84] - Datum coordinate is defined in (default WGS 84).
       * @returns {Utm} Parsed UTM coordinate.
       * @throws  {TypeError} Invalid UTM coordinate.
       *
       * @example
       *   const utmCoord = Utm.parse('31 N 448251 5411932');
       *   // utmCoord: {zone: 31, hemisphere: 'N', easting: 448251, northing: 5411932 }
       */
      static Utm parse(const std::string& utmCoord, std::optional<Datum> datum = LatLonEllipsoidal::datums().WGS84);


      /**
       * Converts UTM zone/easting/northing coordinate to latitude/longitude.
       *
       * Implements Karney’s method, using Krüger series to order n⁶, giving results accurate to 5nm
       * for distances up to 3900km from the central meridian.
       *
       * @param   {Utm} utmCoord - UTM coordinate to be converted to latitude/longitude.
       * @returns {LatLon} Latitude/longitude of supplied grid reference.
       *
       * @example
       *   const grid = new Utm(31, 'N', 448251.795, 5411932.678);
       *   const latlong = grid.toLatLon(); // 48°51′29.52″N, 002°17′40.20″E
       */
      LatLonUtm toLatLon() const;


      /**
       * Returns a string representation of a UTM coordinate.
       *
       * To distinguish from MGRS grid zone designators, a space is left between the zone and the
       * hemisphere.
       *
       * Note that UTM coordinates get rounded, not truncated (unlike MGRS grid references).
       *
       * @param   {number} [digits=0] - Number of digits to appear after the decimal point (3 ≡ mm).
       * @returns {string} A string representation of the coordinate.
       *
       * @example
       *   const utm = new Utm('31', 'N', 448251, 5411932).toString(4);  // 31 N 448251.0000 5411932.0000
       */
      [[nodiscard]] std::string toString(int dp = 0) const;


      [[nodiscard]] double northing() const { return m_northing; }
      [[nodiscard]] double easting() const { return m_easting; }

   //protected:
      int m_zone;
      Hemisphere m_hemisphere;
      double m_easting;
      double m_northing;
      std::optional<Datum> m_datum;
      std::optional<double> m_convergence;
      std::optional<double> m_scale;
   };
}

#endif // UTM_H
