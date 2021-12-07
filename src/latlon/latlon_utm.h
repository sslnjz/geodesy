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

#ifndef LATLON_UTM_H
#define LATLON_UTM_H

#include "latlon_ellipsoidal.h"

namespace geodesy
{
   /**
    * Extends LatLon with method to convert LatLon points to UTM coordinates.
    *
    * @extends LatLon
    */
   class Utm;
   class LatLonUtm : public LatLonEllipsoidal
   {
   public:
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
      LatLonUtm(double lat, double lon, double height = 0.0,
         std::optional<Datum> datum = std::nullopt,
         std::optional<ReferenceFrame> reference = std::nullopt,
         std::optional<std::string> epoch = std::nullopt);
      ~LatLonUtm() override = default;


      /**
       * Converts latitude/longitude to UTM coordinate.
       *
       * Implements Karney’s method, using Krüger series to order n⁶, giving results accurate to 5nm
       * for distances up to 3900km from the central meridian.
       *
       * @param   {number} [zoneOverride] - Use specified zone rather than zone within which point lies;
       *          note overriding the UTM zone has the potential to result in negative eastings, and
       *          perverse results within Norway/Svalbard exceptions.
       * @returns {Utm} UTM coordinate.
       * @throws  {TypeError} Latitude outside UTM limits.
       *
       * @example
       *   const latlong = new LatLon(48.8582, 2.2945);
       *   const utmCoord = latlong.toUtm(); // 31 N 448252 5411933
       */
      Utm toUtm(std::optional<int> zoneOverride = std::nullopt) const;

      /**
       * @brief Set Meridian convergence.
       * @param {number} convergence (bearing of grid north clockwise from true north), in degrees
      */
      void setConvergence(double convergence);
      [[nodiscard]] double getConvergence() const { return m_convergence; }
      /**
       * @brief Grid scale factor
       * @param scale set scale factor
      */
      void setScale(double scale);
      [[nodiscard]] double getScale() const { return m_scale; }

   private:
      double m_convergence = 0.0;
      double m_scale = 0.0;
   };

}

#endif // LATLON_UTM_H
