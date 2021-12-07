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

#ifndef UTM_MGRS_H
#define UTM_MGRS_H

#include "utm.h"

namespace geodesy
{
   /**
    * Extends Utm with method to convert UTM coordinate to MGRS reference.
    *
    * @extends Utm
    */
   class Mgrs;
   class UtmMgrs : public Utm
   {
   public:
      UtmMgrs(int zone, Hemisphere h, double easting, double northing, std::optional<Datum> datum = LatLonEllipsoidal::datums().WGS84);

      /**
       * Converts UTM coordinate to MGRS reference.
       *
       * @returns {Mgrs}
       * @throws  {TypeError} Invalid UTM coordinate.
       *
       * @example
       *   const utmCoord = new Utm(31, 'N', 448251, 5411932);
       *   const mgrsRef = utmCoord.toMgrs(); // 31U DQ 48251 11932
       */
      Mgrs toMgrs() const;
   };


}

#endif // UTM_MGRS_H
