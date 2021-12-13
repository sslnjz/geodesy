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
#ifndef LATLON_OSGRIDREF_H
#define LATLON_OSGRIDREF_H

#include "latlon_ellipsoidal_datum.h"
#include "ellipsoids.h"

namespace geodesy
{
   static const struct _nationalGrid
   {
      inline static const struct _originGrid  // true origin of grid 49°N,2°W on OSGB36 datum
      {
         double lat;
         double lon;
      } &trueOrigin = *new _originGrid { 49, -2.0 };

      inline static const struct _falseOrign  // easting & northing of false origin, metres from true origin
      {
         double easting;
         double northing;
      } &falseOrigin = *new _falseOrign {-400e3, 100e3};

      double scaleFactor = 0.9996012717;      // scale factor on central meridian
      Ellipsoid ellipsoid = LatLonEllipsoidalDatum::ellipsoids().Airy1830;

   } &nationalGrid = *new _nationalGrid;
   // note Irish National Grid uses t/o 53°30′N, 8°W, f/o 200kmW, 250kmS, scale factor 1.000035, on Airy 1830 Modified ellipsoid

   class OsGridRef;
   class LatLonOsGridRef : public LatLonEllipsoidalDatum
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
      LatLonOsGridRef(double lat, double lon, double height = 0.0,
         std::optional<Datum> datum = std::nullopt);
      ~LatLonOsGridRef() override = default;

      /**
       * Converts latitude/longitude to Ordnance Survey grid reference easting/northing coordinate.
       *
       * @returns {OsGridRef} OS Grid Reference easting/northing.
       *
       * @example
       *   const grid = new LatLon(52.65798, 1.71605).toOsGrid(); // TG 51409 13177
       *   // for conversion of (historical) OSGB36 latitude/longitude point:
       *   const grid = new LatLon(52.65798, 1.71605).toOsGrid(LatLon.datums.OSGB36);
       */
      OsGridRef toOsGrid();

      /**
       * Override LatLonEllipsoidal.convertDatum() with version which returns LatLon_OsGridRef.
       */
      LatLonOsGridRef convertDatum(const Datum& toDatum) const;
   };
}

#endif // LATLON_OSGRIDREF_H
