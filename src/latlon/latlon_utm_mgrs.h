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

#ifndef LATLON_UTM_MGRS_H
#define LATLON_UTM_MGRS_H

#include "latlon_utm.h"
#include "utm_mgrs.h"

namespace geodesy
{
   /**
    * Latitude/longitude point with an MGRS-capable UTM conversion bridge.
    */
   class LatLonUtmMgrs : public LatLonUtm
   {
   public:
      LatLonUtmMgrs(double lat, double lon, double height = 0.0,
         std::optional<Datum> datum = std::nullopt,
         std::optional<ReferenceFrame> reference = std::nullopt,
         std::optional<std::string> epoch = std::nullopt);

      ~LatLonUtmMgrs() override = default;

      /**
       * Converts latitude/longitude in degrees to a UTM value that can be formatted as MGRS.
       */
      [[nodiscard]] UtmMgrs toUtm(std::optional<int> zoneOverride = std::nullopt) const;
   };
}

#endif // LATLON_UTM_MGRS_H
