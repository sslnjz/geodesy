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

#ifndef MGRS_H
#define MGRS_H

#include "ellipsoids.h"
#include "latlon_ellipsoidal.h"

#include <string>

namespace geodesy
{
   inline constexpr const char* latBands = "CDEFGHJKLMNPQRSTUVWXX";
   inline constexpr const char* e100kLetters[] = { "ABCDEFGH", "JKLMNPQR", "STUVWXYZ" };
   inline constexpr const char* n100kLetters[] = { "ABCDEFGHJKLMNPQRSTUV", "FGHJKLMNPQRSTUVABCDE" };

   class UtmMgrs;

   /**
    * Military Grid Reference System value type, backed by a UTM grid zone, latitude band,
    * 100 km square letters, and metre offsets within that square.
    */
   class Mgrs
   {
   public:
      Mgrs(int zone, char band, char e100k, char n100k, double easting,
           double northing, Datum datum = LatLonEllipsoidal::datums().WGS84);

      /**
       * Converts this MGRS grid reference to the UTM coordinate at the south-west corner
       * of the represented grid square.
       */
      [[nodiscard]] UtmMgrs toUtm() const;

      /**
       * Parses spaced or military-style MGRS text. Easting and northing digit counts must
       * match so the precision of the square is unambiguous.
       */
      [[nodiscard]] static Mgrs parse(const std::string& mgrsGridRef);

      /**
       * Formats an MGRS reference using an even number of easting+northing digits in the
       * range 0..10. MGRS digits are truncated because the text names a containing square.
       */
      [[nodiscard]] std::string toString(unsigned int digits = 10) const;

      [[nodiscard]] int zone() const { return m_zone; }
      [[nodiscard]] char band() const { return m_band; }
      [[nodiscard]] char e100k() const { return m_e100k; }
      [[nodiscard]] char n100k() const { return m_n100k; }
      [[nodiscard]] double easting() const { return m_easting; }
      [[nodiscard]] double northing() const { return m_northing; }
      [[nodiscard]] Datum datum() const { return m_datum; }

   private:
      int m_zone;
      char m_band;
      char m_e100k;
      char m_n100k;
      double m_easting;
      double m_northing;
      Datum m_datum;
   };
}

#endif // MGRS_H
