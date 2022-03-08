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

#ifndef MGRS_H
#define MGRS_H

#include <string>

#include "ellipsoids.h"
#include "latlon_ellipsoidal.h"

namespace geodesy
{
   /*
    * Latitude bands C..X 8° each, covering 80°S to 84°N
    */
   const inline std::string latBands = "CDEFGHJKLMNPQRSTUVWXX"; // X is repeated for 80-84°N

   /*
    * 100km grid square column ('e') letters repeat every third zone
    */
   const inline std::string e100kLetters[] = 
   {
      "ABCDEFGH",
      "JKLMNPQR",
      "STUVWXYZ"
   };

   /*
    * 100km grid square row ('n') letters repeat every other zone
    */
   const inline std::string n100kLetters[2] = {
      "ABCDEFGHJKLMNPQRSTUV",
      "FGHJKLMNPQRSTUVABCDE"
   };


   /**
    * Military Grid Reference System (MGRS/NATO) grid references provides geocoordinate references
    * covering the entire globe, based on UTM projections.
    *
    * MGRS references comprise a grid zone designator, a 100km square identification, and an easting
    * and northing (in metres); e.g. '31U DQ 48251 11932'.
    *
    * Depending on requirements, some parts of the reference may be omitted (implied), and
    * eastings/northings may be given to varying resolution.
    *
    * qv www.fgdc.gov/standards/projects/FGDC-standards-projects/usng/fgdc_std_011_2001_usng.pdf
    *
    * @module mgrs
    */
   class UtmMgrs;
   class Mgrs
   {
   public:
      /**
        * Creates an Mgrs grid reference object.
        *
        * @param  {number} zone - 6° longitudinal zone (1..60 covering 180°W..180°E).
        * @param  {string} band - 8° latitudinal band (C..X covering 80°S..84°N).
        * @param  {string} e100k - First letter (E) of 100km grid square.
        * @param  {string} n100k - Second letter (N) of 100km grid square.
        * @param  {number} easting - Easting in metres within 100km grid square.
        * @param  {number} northing - Northing in metres within 100km grid square.
        * @param  {LatLon.datums} [datum=WGS84] - Datum UTM coordinate is based on.
        * @throws {RangeError}  Invalid MGRS grid reference.
        *
        * @example
        *   import Mgrs from '/js/geodesy/mgrs.js';
        *   const mgrsRef = new Mgrs(31, 'U', 'D', 'Q', 48251, 11932); // 31U DQ 48251 11932
        */
      Mgrs(int zone, char band, char e100k, char n100k, double easting,
           double northing, Datum datum = LatLonEllipsoidal::datums().WGS84);


      /**
       * Converts MGRS grid reference to UTM coordinate.
       *
       * Grid references refer to squares rather than points (with the size of the square indicated
       * by the precision of the reference); this conversion will return the UTM coordinate of the SW
       * corner of the grid reference square.
       *
       * @returns {Utm} UTM coordinate of SW corner of this MGRS grid reference.
       *
       * @example
       *   const mgrsRef = Mgrs.parse('31U DQ 48251 11932');
       *   const utmCoord = mgrsRef.toUtm(); // 31 N 448251 5411932
       */
      UtmMgrs toUtm();

      /**
     * Parses string representation of MGRS grid reference.
     *
     * An MGRS grid reference comprises (space-separated)
     *  - grid zone designator (GZD)
     *  - 100km grid square letter-pair
     *  - easting
     *  - northing.
     *
     * @param   {string} mgrsGridRef - String representation of MGRS grid reference.
     * @returns {Mgrs}   Mgrs grid reference object.
     * @throws  {Error}  Invalid MGRS grid reference.
     *
     * @example
     *   const mgrsRef = Mgrs.parse('31U DQ 48251 11932');
     *   const mgrsRef = Mgrs.parse('31UDQ4825111932');
     *   //  mgrsRef: { zone:31, band:'U', e100k:'D', n100k:'Q', easting:48251, northing:11932 }
     */
      static Mgrs parse(const std::string& mgrsGridRef);

      /**
        * Returns a string representation of an MGRS grid reference.
        *
        * To distinguish from civilian UTM coordinate representations, no space is included within the
        * zone/band grid zone designator.
        *
        * Components are separated by spaces: for a military-style unseparated string, use
        *   Mgrs.toString().replace(/ /g, '');
        *
        * Note that MGRS grid references get truncated, not rounded (unlike UTM coordinates); grid
        * references indicate a bounding square, rather than a point, with the size of the square
        * indicated by the precision - a precision of 10 indicates a 1-metre square, a precision of 4
        * indicates a 1,000-metre square (hence 31U DQ 48 11 indicates a 1km square with SW corner at
        * 31 N 448000 5411000, which would include the 1m square 31U DQ 48251 11932).
        *
        * @param   {number}     [digits=10] - Precision of returned grid reference (eg 4 = km, 10 = m).
        * @returns {string}     This grid reference in standard format.
        * @throws  {RangeError} Invalid precision.
        *
        * @example
        *   const mgrsStr = new Mgrs(31, 'U', 'D', 'Q', 48251, 11932).toString(); // 31U DQ 48251 11932
        */
      std::string toString(unsigned digits=10);

   private:
      int m_zone;
      char m_band;
      std::string m_e100k;
      std::string m_n100k;
      double m_easting;
      double m_northing;
      Datum m_datum;
   };
}


#endif // MGRS_H
