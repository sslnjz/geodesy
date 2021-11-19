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
#ifndef LATLON_SPHERICAL_H
#define LATLON_SPHERICAL_H

#include <stdexcept>
#include <string>
#include <type_traits>
#include <cmath>

#include "dms.h"

namespace geodesy
{
    // note greek letters (e.g. φ, λ, θ) are used for angles in radians to distinguish from angles in
    // degrees (e.g. lat, lon, brng)
    /* LatLonSpherical - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
    /**
     * Latitude/longitude points on a spherical model earth, and methods for calculating distances,
     * bearings, destinations, etc on (orthodromic) great-circle paths and (loxodromic) rhumb lines.
     */
   class LatLonSpherical
   {
   public:
      LatLonSpherical(double lat, double lon);

      [[nodiscard]] double lat() const;
      [[nodiscard]] double latitude() const;

      void setLat(double lat);
      void setLatitude(double lat);

      /**
       * Longitude in degrees east from international reference meridian (including aliases lon, lng,
       * longitude): can be set as numeric or hexagesimal (deg-min-sec); returned as numeric.
       */
      [[nodiscard]] double lon() const;
      [[nodiscard]] double lng() const;
      [[nodiscard]] double longitude() const;

      void setLon(double lon);
      void setLng(double lon);
      void setLongitude(double lon);


      /**
       * Conversion factor metres to kilometres.
       * Conversion factors; 1000 * LatLon.metresToKm gives 1.
       */
      static inline double getMetresToKm();
      /** Conversion factors; 1000 * LatLon.metresToMiles gives 0.621371192237334. */
      static inline double getMetresToMiles();
      /** Conversion factors; 1000 * LatLon.metresToMiles gives 0.5399568034557236. */
      static inline double getMetresToNauticalMiles();

      /**
       * Returns the distance along the surface of the earth from ‘this’ point to destination point.
       *
       * Uses haversine formula: a = sin²(Δφ/2) + cosφ1·cosφ2 · sin²(Δλ/2); d = 2 · atan2(√a, √(a-1)).
       *
       * @param   {LatLon} point - Latitude/longitude of destination point.
       * @param   {number} [radius=6371e3] - Radius of earth (defaults to mean radius in meters).
       * @returns {number} Distance between this point and destination point, in same units as radius.
       * @throws  {TypeError} Invalid radius.
       *
       * @example
       *   const double p1 = new LatLon(52.205, 0.119);
       *   const double p2 = new LatLon(48.857, 2.351);
       *   const double d = p1.distanceTo(p2);       // 404.3×10³ m
       *   const double m = p1.distanceTo(p2, 3959); // 251.2 miles
       */
      [[nodiscard]] double distanceTo(const LatLonSpherical& point, double radius = 6371e3) const;


      /**
       * Returns the initial bearing from ‘this’ point to destination point.
       *
       * @param   {LatLon} point - Latitude/longitude of destination point.
       * @returns {number} Initial bearing in degrees from north (0°..360°).
       *
       * @example
       *   const double p1 = new LatLonSpherical(52.205, 0.119);
       *   const double p2 = new LatLonSpherical(48.857, 2.351);
       *   const double b1 = p1.initialBearingTo(p2); // 156.2°
       */
      [[nodiscard]] double initialBearingTo(const LatLonSpherical& point) const;


      /**
       * Returns final bearing arriving at destination point from ‘this’ point; the final bearing will
       * differ from the initial bearing by varying degrees according to distance and latitude.
       *
       * @param   {LatLon} point - Latitude/longitude of destination point.
       * @returns {number} Final bearing in degrees from north (0°..360°).
       *
       * @example
       *   const double p1 = new LatLonSpherical(52.205, 0.119);
       *   const double p2 = new LatLonSpherical(48.857, 2.351);
       *   const double b2 = p1.finalBearingTo(p2); // 157.9°
       */
      [[nodiscard]] double finalBearingTo(const LatLonSpherical& point) const;

      /**
       * Returns the midpoint between ‘this’ point and destination point.
       *
       * @param   {LatLonSpherical} point - Latitude/longitude of destination point.
       * @returns {LatLonSpherical} Midpoint between this point and destination point.
       *
       * @example
       *   const double p1 = new LatLonSpherical(52.205, 0.119);
       *   const double p2 = new LatLonSpherical(48.857, 2.351);
       *   const LatLonSpherical pMid = p1.midpointTo(p2); // 50.5363°N, 001.2746°E
       */
      [[nodiscard]] LatLonSpherical midpointTo(const LatLonSpherical& point) const;


      /**
       * Returns the point at given fraction between ‘this’ point and given point.
       *
       * @param   {LatLonSpherical} point - Latitude/longitude of destination point.
       * @param   {number} fraction - Fraction between the two points (0 = this point, 1 = specified point).
       * @returns {LatLonSpherical} Intermediate point between this point and destination point.
       *
       * @example
       *   const p1 = new LatLon(52.205, 0.119);
       *   const p2 = new LatLon(48.857, 2.351);
       *   const pInt = p1.intermediatePointTo(p2, 0.25); // 51.3721°N, 000.7073°E
       */
      [[nodiscard]] LatLonSpherical intermediatePointTo(const LatLonSpherical& point, double fraction) const;

      /**
       * Returns the destination point from ‘this’ point having travelled the given distance on the
       * given initial bearing (bearing normally varies around path followed).
       *
       * @param   {number} distance - Distance travelled, in same units as earth radius (default: metres).
       * @param   {number} bearing - Initial bearing in degrees from north.
       * @param   {number} [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
       * @returns {LatLonSpherical} Destination point.
       *
       * @example
       *   const LatLonSpherical p1 = new LatLonSpherical(51.47788, -0.00147);
       *   const LatLonSpherical p2 = p1.destinationPoint(7794, 300.7); // 51.5136°N, 000.0983°W
       */
      [[nodiscard]] LatLonSpherical destinationPoint(double distance, double bearing, double radius = 6371e3) const;

      /**
       * Returns the point of intersection of two paths defined by point and bearing.
       *
       * @param   {LatLon}      p1 - First point.
       * @param   {number}      brng1 - Initial bearing from first point.
       * @param   {LatLon}      p2 - Second point.
       * @param   {number}      brng2 - Initial bearing from second point.
       * @returns {LatLon}      Destination point (exception if no unique intersection defined).
       * @throws  {runtime_error} no unique intersection defined.
       *
       * @example
       *   const p1 = new LatLonSpherical(51.8853, 0.2545), brng1 = 108.547;
       *   const p2 = new LatLonSpherical(49.0034, 2.5735), brng2 =  32.435;
       *   const pInt = LatLon.intersection(p1, brng1, p2, brng2); // 50.9078°N, 004.5084°E
       */
      static LatLonSpherical intersection(const LatLonSpherical& p1, double brng1, const LatLonSpherical& p2,
                                          double brng2);

      /**
       * Returns (signed) distance from ‘this’ point to great circle defined by start-point and
       * end-point.
       *
       * @param   {LatLonSpherical} pathStart - Start point of great circle path.
       * @param   {LatLonSpherical} pathEnd - End point of great circle path.
       * @param   {number} [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
       * @returns {number} Distance to great circle (-ve if to left, +ve if to right of path).
       *
       * @example
       *   const LatLonSpherical pCurrent = LatLonSpherical(53.2611, -0.7972);
       *   const LatLonSpherical p1 = LatLonSpherical(53.3206, -1.7297);
       *   const LatLonSpherical p2 = LatLonSpherical(53.1887, 0.1334);
       *   const double d = pCurrent.crossTrackDistanceTo(p1, p2);  // -307.5 m
       */
      [[nodiscard]] double crossTrackDistanceTo(const LatLonSpherical& pathStart, const LatLonSpherical& pathEnd,
                                  double radius = 6371e3) const;


      /**
       * Returns how far ‘this’ point is along a path from from start-point, heading towards end-point.
       * That is, if a perpendicular is drawn from ‘this’ point to the (great circle) path, the
       * along-track distance is the distance from the start point to where the perpendicular crosses
       * the path.
       *
       * @param   {LatLonSpherical} pathStart - Start point of great circle path.
       * @param   {LatLonSpherical} pathEnd - End point of great circle path.
       * @param   {number} [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
       * @returns {number} Distance along great circle to point nearest ‘this’ point.
       *
       * @example
       *   const auto pCurrent = new LatLonSpherical(53.2611, -0.7972);
       *   const auto p1 = new LatLonSpherical(53.3206, -1.7297);
       *   const auto p2 = new LatLonSpherical(53.1887,  0.1334);
       *   const auto d = pCurrent.alongTrackDistanceTo(p1, p2);  // 62.331 km
       */
      [[nodiscard]] double alongTrackDistanceTo(const LatLonSpherical& pathStart, const LatLonSpherical& pathEnd,
                                  double radius = 6371e3) const;


      /**
       * Returns maximum latitude reached when travelling on a great circle on given bearing from
       * ‘this’ point (‘Clairaut’s formula’). Negate the result for the minimum latitude (in the
       * southern hemisphere).
       *
       * The maximum latitude is independent of longitude; it will be the same for all points on a
       * given latitude.
       *
       * @param   {number} bearing - Initial bearing.
       * @returns {number} Maximum latitude reached.
       */
      [[nodiscard]] double maxLatitude(double bearing) const;

      /**
       * Returns the pair of meridians at which a great circle defined by two points crosses the given
       * latitude. If the great circle doesn't reach the given latitude, null is returned.
       *
       * @param   {LatLon}      point1 - First point defining great circle.
       * @param   {LatLon}      point2 - Second point defining great circle.
       * @param   {number}      latitude - Latitude crossings are to be determined for.
       * @returns {Object|null} Object containing { lon1, lon2 } or null if given latitude not reached.
       */
      static std::pair<double, double> crossingParallels(const LatLonSpherical& point1, const LatLonSpherical& point2,
                                                         double latitude);


      /* Rhumb - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
      /**
       * Returns the distance travelling from ‘this’ point to destination point along a rhumb line.
       *
       * @param   {LatLon} point - Latitude/longitude of destination point.
       * @param   {number} [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
       * @returns {number} Distance in km between this point and destination point (same units as radius).
       *
       * @example
       *   const p1 = new LatLon(51.127, 1.338);
       *   const p2 = new LatLon(50.964, 1.853);
       *   const d = p1.distanceTo(p2); //  40.31 km
       */
      [[nodiscard]] double rhumbDistanceTo(const LatLonSpherical& point, double radius = 6371e3) const;

      /**
       * Returns the bearing from ‘this’ point to destination point along a rhumb line.
       *
       * @param   {LatLon}    point - Latitude/longitude of destination point.
       * @returns {number}    Bearing in degrees from north.
       *
       * @example
       *   const p1 = new LatLon(51.127, 1.338);
       *   const p2 = new LatLon(50.964, 1.853);
       *   const d = p1.rhumbBearingTo(p2); // 116.7°
       */
      [[nodiscard]] double rhumbBearingTo(const LatLonSpherical& point) const;


      /**
       * Returns the destination point having travelled along a rhumb line from ‘this’ point the given
       * distance on the given bearing.
       *
       * @param   {number} distance - Distance travelled, in same units as earth radius (default: metres).
       * @param   {number} bearing - Bearing in degrees from north.
       * @param   {number} [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
       * @returns {LatLon} Destination point.
       *
       * @example
       *   const p1 = new LatLon(51.127, 1.338);
       *   const p2 = p1.rhumbDestinationPoint(40300, 116.7); // 50.9642°N, 001.8530°E
       */
      [[nodiscard]] LatLonSpherical rhumbDestinationPoint(double distance, double bearing, double radius = 6371e3) const;

      /**
       * Returns the loxodromic midpoint (along a rhumb line) between ‘this’ point and second point.
       *
       * @param   {LatLon} point - Latitude/longitude of second point.
       * @returns {LatLon} Midpoint between this point and second point.
       *
       * @example
       *   const p1 = new LatLon(51.127, 1.338);
       *   const p2 = new LatLon(50.964, 1.853);
       *   const pMid = p1.rhumbMidpointTo(p2); // 51.0455°N, 001.5957°E
       */
      [[nodiscard]] LatLonSpherical rhumbMidpointTo(const LatLonSpherical& point) const;

      /* Area - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /**
     * Calculates the area of a spherical polygon where the sides of the polygon are great circle
     * arcs joining the vertices.
     *
     * @param   {LatLon[]} polygon - Array of points defining vertices of the polygon.
     * @param   {number}   [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
     * @returns {number}   The area of the polygon in the same units as radius.
     *
     * @example
     *   const polygon = [new LatLon(0,0), new LatLon(1,0), new LatLon(0,1)];
     *   const area = LatLon.areaOf(polygon); // 6.18e9 m²
     */
      static double areaOf(std::vector<LatLonSpherical>& polygon, double radius = 6371e3);

      /**
       * Returns a string representation of ‘this’ point, formatted as degrees, degrees+minutes, or
       * degrees+minutes+seconds.
       *
       * @param   {enum} [eFormat=D] - Format point as 'D', 'DM', 'DMS', or 'N' for signed numeric.
       * @returns {string} Comma-separated formatted latitude/longitude.
       *
       * @example
       *   const greenwich = new LatLon(51.47788, -0.00147);
       *   const d = greenwich.toString();                        // 51.4779°N, 000.0015°W
       *   const dms = greenwich.toString('dms', 2);              // 51°28′40.37″N, 000°00′05.29″W
       *   const [lat, lon] = greenwich.toString('n').split(','); // 51.4779, -0.0015
       */
      [[nodiscard]] std::string toString(Dms::eFormat e = Dms::D) const;

      /**
       * Converts ‘this’ point to a GeoJSON object string.
       *
       * @returns {string} this point as a GeoJSON ‘Point’ string.
       *    { type: "Point", coordinates: [ lon, lat ] }
       */
      [[nodiscard]] std::string toGeoJSON() const;


      /**
       * Checks if another point is equal to ‘this’ point.
       *
       * @param   {LatLon} p1 - Point to be compared against p2.
       * @param   {LatLon} p2 - Point to be compared against p1.
       * @returns {bool}   True if points have identical latitude and longitude values.
       *
       */
      friend inline bool operator==(const LatLonSpherical& p1, const LatLonSpherical& p2);

      /**
       * Checks if another point is NOT equal to ‘this’ point.
       *
       * @param   {LatLon} p1 - Point to be compared against p2.
       * @param   {LatLon} p2 - Point to be compared against p1.
       * @returns {bool}   True if points have different latitude or longitude values.
       *
       */
      friend inline bool operator!=(const LatLonSpherical& p1, const LatLonSpherical& p2);

   private:
      double m_lat; // Latitude in degrees north from equator
      double m_lon; // Longitude in degrees east from international reference meridian
   };
}


#endif // LATLON_SPHERICAL_H
