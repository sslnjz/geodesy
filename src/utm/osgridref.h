
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
#ifndef OSGRIDREF_H
#define OSGRIDREF_H

#include "latlon_ellipsoidal_datum.h"

#include <string>

namespace geodesy
{
   /**
    * British National Grid easting/northing in metres from the OS false origin.
    */
   class LatLonOsGridRef;
   class OsGridRef
   {
   public:
      /**
       * Creates an OsGridRef object.
       *
       * @param {number} easting - Easting in metres from OS Grid false origin.
       * @param {number} northing - Northing in metres from OS Grid false origin.
       *
       * @example
       *   import OsGridRef from '/js/geodesy/osgridref.js';
       *   const gridref = new OsGridRef(651409, 313177);
       */
      OsGridRef(double easting, double northing);

      [[nodiscard]] double easting() const;
      [[nodiscard]] double northing() const;

      /**
       * Converts ‘this’ Ordnance Survey Grid Reference easting/northing coordinate to latitude/longitude
       * (SW corner of grid square).
       *
       * While OS Grid References are based on OSGB-36, the Ordnance Survey have deprecated the use of
       * OSGB-36 for latitude/longitude coordinates (in favour of WGS-84), hence this function returns
       * WGS-84 by default, with OSGB-36 as an option. See www.ordnancesurvey.co.uk/blog/2014/12/2.
       *
       * Note formulation implemented here due to Thomas, Redfearn, etc is as published by OS, but is
       * inferior to Krüger as used by e.g. Karney 2011.
       *
       * @param   {LatLon.datum} [datum=WGS84] - Datum to convert grid reference into.
       * @returns {LatLon}       Latitude/longitude of supplied grid reference.
       *
       * @example
       *   const gridref = new OsGridRef(651409.903, 313177.270);
       *   const pWgs84 = gridref.toLatLon();                    // 52°39′28.723″N, 001°42′57.787″E
       *   // to obtain (historical) OSGB36 lat/lon point:
       *   const pOsgb = gridref.toLatLon(LatLon.datums.OSGB36); // 52°39′27.253″N, 001°43′04.518″E
       */
      LatLonOsGridRef toLatLon(Datum datum = LatLonEllipsoidalDatum::datums().WGS84) const;

      /**
     * Parses grid reference to OsGridRef object.
     *
     * Accepts standard grid references (eg 'SU 387 148'), with or without whitespace separators, from
     * two-digit references up to 10-digit references (1m × 1m square), or fully numeric comma-separated
     * references in metres (eg '438700,114800').
     *
     * @param   {string}    gridref - Standard format OS Grid Reference.
     * @returns {OsGridRef} Numeric version of grid reference in metres from false origin (SW corner of
     *   supplied grid square).
     * @throws  {Error}     Invalid grid reference.
     *
     * @example
     *   const grid = OsGridRef.parse('TG 51409 13177'); // grid: { easting: 651409, northing: 313177 }
     */
      static OsGridRef parse(const std::string& gridref);

      /**
       * Converts ‘this’ numeric grid reference to standard OS Grid Reference.
       *
       * @param   {number} [digits=10] - Precision of returned grid reference (10 digits = metres);
       *   digits=0 will return grid reference in numeric format.
       * @returns {string} This grid reference in standard format.
       *
       * @example
       *   const gridref = new OsGridRef(651409, 313177).toString(8); // 'TG 5140 1317'
       *   const gridref = new OsGridRef(651409, 313177).toString(0); // '651409,313177'
       */
      [[nodiscard]] std::string toString(int digits = 10) const;

   private:
      double m_easting;
      double m_northing;
   };
}


#endif // OSGRIDREF_H
