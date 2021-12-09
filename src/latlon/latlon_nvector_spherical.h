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
#ifndef LATLON_NVECTOR_SPHERICAL_H
#define LATLON_NVECTOR_SPHERICAL_H

#include "dms.h"
#include "latlon.h"
#include "nvector_spherical.h"

namespace geodesy
{
   class vector3d;
   class NvectorSpherical;
   class LatLonNvectorSpherical : public LatLon
   {
   public:
      LatLonNvectorSpherical();
      LatLonNvectorSpherical(double lat, double lon);
      LatLonNvectorSpherical(const std::string& lat, const std::string& lon);


      /**
        * Converts ‘this’ latitude/longitude point to an n-vector (normal to earth's surface).
        *
        * @returns {Nvector} Normalised n-vector representing lat/lon point.
        *
        * @example
        *   const p = new LatLon(45, 45);
        *   const v = p.toNvector();      // [0.5000,0.5000,0.7071]
        */
      [[nodiscard]] NvectorSpherical toNvector() const;

      /**
       * Vector normal to great circle obtained by heading on given bearing from ‘this’ point.
       *
       * Direction of vector is such that initial bearing vector b = c × n, where n is an n-vector
       * representing ‘this’ (start) point.
       *
       * @private
       * @param   {number}   bearing - Compass bearing in degrees.
       * @returns {Vector3d} Normalised vector representing great circle.
       *
       * @example
       *   const p1 = new LatLon(53.3206, -1.7297);
       *   const gc = p1.greatCircle(96.0);         // [-0.794,0.129,0.594]
       */
      [[nodiscard]] vector3d greatCircle(double bearing) const;

      /**
       * Returns the distance on the surface of the sphere from ‘this’ point to destination point.
       *
       * @param   {LatLon}    point - Latitude/longitude of destination point.
       * @param   {number}    [radius=6371e3] - Radius of earth (defaults to mean radius in metres).
       * @returns {number}    Distance between this point and destination point, in same units as radius.
       * @throws  {TypeError} Invalid point/radius.
       *
       * @example
       *   const p1 = new LatLon(52.205, 0.119);
       *   const p2 = new LatLon(48.857, 2.351);
       *   const d = p1.distanceTo(p2);          // 404.3 km
       */
      [[nodiscard]] double distanceTo(const LatLonNvectorSpherical& point, double radius = 6371e3) const;

      /**
       * Returns the initial bearing from ‘this’ point to destination point.
       *
       * @param   {LatLon}    point - Latitude/longitude of destination point.
       * @returns {number}    Initial bearing in degrees from north (0°..360°).
       * @throws  {TypeError} Invalid point.
       *
       * @example
       *   const p1 = new LatLon(52.205, 0.119);
       *   const p2 = new LatLon(48.857, 2.351);
       *   const b1 = p1.initialBearingTo(p2);   // 156.2°
       */
      [[nodiscard]] double initialBearingTo(const LatLonNvectorSpherical& point) const;

      /**
       * Returns final bearing arriving at destination point from ‘this’ point; the final bearing will
       * differ from the initial bearing by varying degrees according to distance and latitude.
       *
       * @param   {LatLon}    point - Latitude/longitude of destination point.
       * @returns {number}    Final bearing in degrees from north (0°..360°).
       * @throws  {TypeError} Invalid point.
       *
       * @example
       *   const p1 = new LatLon(52.205, 0.119);
       *   const p2 = new LatLon(48.857, 2.351);
       *   const b2 = p1.finalBearingTo(p2); // 157.9°
       */
      [[nodiscard]] double finalBearingTo(const LatLonNvectorSpherical& point) const;

      /**
       * Returns the midpoint between ‘this’ point and destination point.
       *
       * @param   {LatLon}    point - Latitude/longitude of destination point.
       * @returns {LatLon}    Midpoint between this point and destination point.
       * @throws  {TypeError} Invalid point.
       *
       * @example
       *   const p1 = new LatLon(52.205, 0.119);
       *   const p2 = new LatLon(48.857, 2.351);
       *   const pMid = p1.midpointTo(p2);       // 50.5363°N, 001.2746°E
       */
      [[nodiscard]] LatLonNvectorSpherical midpointTo(const LatLonNvectorSpherical& point) const;

      /**
       * Returns the point at given fraction between ‘this’ point and given point.
       *
       * @param   {LatLon}    point - Latitude/longitude of destination point.
       * @param   {number}    fraction - Fraction between the two points (0 = this point, 1 = specified point).
       * @returns {LatLon}    Intermediate point between this point and destination point.
       * @throws  {TypeError} Invalid point/fraction.
       *
       * @example
       *   const p1 = new LatLon(52.205, 0.119);
       *   const p2 = new LatLon(48.857, 2.351);
       *   const pInt = p1.intermediatePointTo(p2, 0.25); // 51.3721°N, 000.7072°E
       */
      [[nodiscard]] LatLonNvectorSpherical intermediatePointTo(const LatLonNvectorSpherical& point, double fraction) const;


      /**
       * Returns the latitude/longitude point projected from the point at given fraction on a straight
       * line between between ‘this’ point and given point.
       *
       * @param   {LatLon}    point - Latitude/longitude of destination point.
       * @param   {number}    fraction - Fraction between the two points (0 = this point, 1 = specified point).
       * @returns {LatLon}    Intermediate point between this point and destination point.
       *
       * @example
       *   const p1 = new LatLon(52.205, 0.119);
       *   const p2 = new LatLon(48.857, 2.351);
       *   const pInt = p1.intermediatePointTo(p2, 0.25); // 51.3723°N, 000.7072°E
       */
      [[nodiscard]] LatLonNvectorSpherical intermediatePointOnChordTo(const LatLonNvectorSpherical& point, double fraction) const;


      /**
       * Returns the destination point from ‘this’ point having travelled the given distance on the
       * given initial bearing (bearing normally varies around path followed).
       *
       * @param   {number} distance - Distance travelled, in same units as earth radius (default: metres).
       * @param   {number} bearing - Initial bearing in degrees from north.
       * @param   {number} [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
       * @returns {LatLon} Destination point.
       *
       * @example
       *   const p1 = new LatLon(51.47788, -0.00147);
       *   const p2 = p1.destinationPoint(7794, 300.7); // 51.5136°N, 000.0983°W
       */
      [[nodiscard]] LatLonNvectorSpherical destinationPoint(double distance, double bearing, double radius = 6371e3) const;


      /**
       * Returns the point of intersection of two paths each defined by point pairs or start point and bearing.
       *
       * @param   {LatLon}        path1start - Start point of first path.
       * @param   {LatLon} path1brngEnd - End point of first path.
       * @param   {LatLon}        path2start - Start point of second path.
       * @param   {LatLon} path2brngEnd - End point of second path.
       * @returns {LatLon}        Destination point (null if no unique intersection defined)
       * @throws  {TypeError}     Invalid parameter.
       *
       * @example
       *   const p1 = new LatLon(51.8853, 0.2545), brng1 = 108.55;
       *   const p2 = new LatLon(49.0034, 2.5735), brng2 =  32.44;
       *   const pInt = LatLon.intersection(p1, brng1, p2, brng2); // 50.9076°N, 004.5086°E
       */
      static LatLonNvectorSpherical intersection(const LatLonNvectorSpherical& path1start, const LatLonNvectorSpherical& path1brngEnd,
         const LatLonNvectorSpherical& path2start, const LatLonNvectorSpherical& path2brngEnd);


      /**
       * Returns the point of intersection of two paths each defined by point pairs or start point and bearing.
       *
       * @param   {LatLon}        path1start - Start point of first path.
       * @param   {LatLon} path1brngEnd - End point of first path.
       * @param   {LatLon}        path2start - Start point of second path.
       * @param   {LatLon} path2brngEnd - End point of second path.
       * @returns {LatLon}        Destination point (null if no unique intersection defined)
       * @throws  {TypeError}     Invalid parameter.
       *
       * @example
       *   const p1 = new LatLon(51.8853, 0.2545), brng1 = 108.55;
       *   const p2 = new LatLon(49.0034, 2.5735), brng2 =  32.44;
       *   const pInt = LatLon.intersection(p1, brng1, p2, brng2); // 50.9076°N, 004.5086°E
       */
      static LatLonNvectorSpherical intersection(const LatLonNvectorSpherical& path1start, double path1brngEnd,
         const LatLonNvectorSpherical& path2start, double path2brngEnd);


      /**
       * Returns the point of intersection of two paths each defined by point pairs or start point and bearing.
       *
       * @param   {LatLon}        path1start - Start point of first path.
       * @param   {LatLon} path1brngEnd - End point of first path.
       * @param   {LatLon}        path2start - Start point of second path.
       * @param   {LatLon} path2brngEnd - End point of second path.
       * @returns {LatLon}        Destination point (null if no unique intersection defined)
       * @throws  {TypeError}     Invalid parameter.
       *
       * @example
       *   const p1 = new LatLon(51.8853, 0.2545), brng1 = 108.55;
       *   const p2 = new LatLon(49.0034, 2.5735), brng2 =  32.44;
       *   const pInt = LatLon.intersection(p1, brng1, p2, brng2); // 50.9076°N, 004.5086°E
       */
      static LatLonNvectorSpherical intersection(const LatLonNvectorSpherical& path1start, const LatLonNvectorSpherical& path1brngEnd,
         const LatLonNvectorSpherical& path2start, double path2brngEnd);


      /**
       * Returns the point of intersection of two paths each defined by point pairs or start point and bearing.
       *
       * @param   {LatLon}        path1start - Start point of first path.
       * @param   {LatLon} path1brngEnd - End point of first path.
       * @param   {LatLon}        path2start - Start point of second path.
       * @param   {LatLon} path2brngEnd - End point of second path.
       * @returns {LatLon}        Destination point (null if no unique intersection defined)
       * @throws  {TypeError}     Invalid parameter.
       *
       * @example
       *   const p1 = new LatLon(51.8853, 0.2545), brng1 = 108.55;
       *   const p2 = new LatLon(49.0034, 2.5735), brng2 =  32.44;
       *   const pInt = LatLon.intersection(p1, brng1, p2, brng2); // 50.9076°N, 004.5086°E
       */
      static LatLonNvectorSpherical intersection(const LatLonNvectorSpherical& path1start, double path1brngEnd,
         const LatLonNvectorSpherical& path2start, const LatLonNvectorSpherical& path2brngEnd);

      /**
       * Returns (signed) distance from ‘this’ point to great circle defined by start-point and end-point/bearing.
       *
       * @param   {LatLon}        pathStart - Start point of great circle path.
       * @param   {LatLon}        pathBrngEnd - End point of great circle path or initial bearing from great circle start point.
       * @param   {number} [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
       * @returns {number}        Distance to great circle (-ve if to left, +ve if to right of path).
       * @throws  {TypeError}     Invalid parameter.
       *
       * @example
       *   const pCurrent = new LatLon(53.2611, -0.7972);
       *
       *   const p1 = new LatLon(53.3206, -1.7297), brng = 96.0;
       *   const d = pCurrent.crossTrackDistanceTo(p1, brng); // Number(d.toPrecision(4)): -305.7
       *
       *   const p1 = new LatLon(53.3206, -1.7297), p2 = new LatLon(53.1887, 0.1334);
       *   const d = pCurrent.crossTrackDistanceTo(p1, p2);   // Number(d.toPrecision(4)): -307.5
       */
      [[nodiscard]] double crossTrackDistanceTo(const LatLonNvectorSpherical& pathStart, const LatLonNvectorSpherical& pathBrngEnd,
                                                double radius = 6371e3) const;


      /**
       * Returns (signed) distance from ‘this’ point to great circle defined by start-point and end-point/bearing.
       *
       * @param   {LatLon}        pathStart - Start point of great circle path.
       * @param   {number}        pathBrngEnd - initial bearing from great circle start point.
       * @param   {number}       [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
       * @returns {number}        Distance to great circle (-ve if to left, +ve if to right of path).
       * @throws  {TypeError}     Invalid parameter.
       *
       * @example
       *   const pCurrent = new LatLon(53.2611, -0.7972);
       *
       *   const p1 = new LatLon(53.3206, -1.7297), brng = 96.0;
       *   const d = pCurrent.crossTrackDistanceTo(p1, brng); // Number(d.toPrecision(4)): -305.7
       *
       *   const p1 = new LatLon(53.3206, -1.7297), p2 = new LatLon(53.1887, 0.1334);
       *   const d = pCurrent.crossTrackDistanceTo(p1, p2);   // Number(d.toPrecision(4)): -307.5
       */
      [[nodiscard]] double crossTrackDistanceTo(const LatLonNvectorSpherical& pathStart, double pathBrngEnd,
         double radius = 6371e3) const;

      /**
       * Returns how far ‘this’ point is along a path from from start-point, heading on bearing or towards
       * end-point. That is, if a perpendicular is drawn from ‘this’ point to the (great circle) path, the
       * along-track distance is the distance from the start point to where the perpendicular crosses the
       * path.
       *
       * @param   {LatLon} pathStart - Start point of great circle path.
       * @param   {LatLon} pathBrngEnd - End point of great circle path.
       * @param   {number} [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
       * @returns {number} Distance along great circle to point nearest ‘this’ point.
       *
       * @example
       *   const pCurrent = new LatLon(53.2611, -0.7972);
       *   const p1 = new LatLon(53.3206, -1.7297);
       *   const p2 = new LatLon(53.1887,  0.1334);
       *   const d = pCurrent.alongTrackDistanceTo(p1, p2);  // 62.331 km
       */
      [[nodiscard]] double alongTrackDistanceTo(const LatLonNvectorSpherical& pathStart, const LatLonNvectorSpherical& pathBrngEnd,
                                  double radius = 6371e3) const;


      /**
       * Returns how far ‘this’ point is along a path from from start-point, heading on bearing or towards
       * end-point. That is, if a perpendicular is drawn from ‘this’ point to the (great circle) path, the
       * along-track distance is the distance from the start point to where the perpendicular crosses the
       * path.
       *
       * @param   {LatLon}        pathStart - Start point of great circle path.
       * @param   {number}        pathBrngEnd - initial bearing from great circle start point.
       * @param   {number}        [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
       * @returns {number}        Distance along great circle to point nearest ‘this’ point.
       *
       * @example
       *   const pCurrent = new LatLon(53.2611, -0.7972);
       *   const p1 = new LatLon(53.3206, -1.7297);
       *   const p2 = new LatLon(53.1887,  0.1334);
       *   const d = pCurrent.alongTrackDistanceTo(p1, p2);  // 62.331 km
       */
      [[nodiscard]] double alongTrackDistanceTo(const LatLonNvectorSpherical& pathStart, double pathBrngEnd,
         double radius = 6371e3) const;


      /**
       * Returns closest point on great circle segment between point1 & point2 to ‘this’ point.
       *
       * If this point is ‘within’ the extent of the segment, the point is on the segment between point1 &
       * point2; otherwise, it is the closer of the endpoints defining the segment.
       *
       * @param   {LatLon} point1 - Start point of great circle segment.
       * @param   {LatLon} point2 - End point of great circle segment.
       * @returns {LatLon} point on segment.
       *
       * @example
       *   const p1 = new LatLon(51.0, 1.0);
       *   const p2 = new LatLon(51.0, 2.0);
       *
       *   const p0 = new LatLon(51.0, 1.9);
       *   const p = p0.nearestPointOnSegment(p1, p2); // 51.0004°N, 001.9000°E
       *   const d = p.distanceTo(p);                  // 42.71 m
       *
       *   const p0 = new LatLon(51.0, 2.1);
       *   const p = p0.nearestPointOnSegment(p1, p2); // 51.0000°N, 002.0000°E
       */
      [[nodiscard]] LatLonNvectorSpherical nearestPointOnSegment(const LatLonNvectorSpherical& point1,
                                                   const LatLonNvectorSpherical& point2) const;

      /**
     * Locates a point given two known locations and bearings from those locations.
     *
     * @param   {LatLon} point1 - First reference point.
     * @param   {number} bearing1 - Bearing (in degrees from north) from first reference point.
     * @param   {LatLon} point2 - Second reference point.
     * @param   {number} bearing2 - Bearing (in degrees from north) from second reference point.
     * @returns {LatLon} Triangulated point.
     *
     * @example
     *   const p1 = new LatLon(50.7175,1.65139), p2 = new LatLon(50.9250,1.7094);
     *   const p = LatLon.triangulate(p1, 333.3508, p2, 310.1414); // 51.1297°N, 001.3214°E
     */
      static LatLonNvectorSpherical triangulate(const LatLonNvectorSpherical& point1, double bearing1,
                                                const LatLonNvectorSpherical& point2, double bearing2);


      /**
    * Locates a latitude/longitude point at given distances from three other points.
    *
    * @param   {LatLon} point1 - First reference point.
    * @param   {number} distance1 - Distance to first reference point (same units as radius).
    * @param   {LatLon} point2 - Second reference point.
    * @param   {number} distance2 - Distance to second reference point (same units as radius).
    * @param   {LatLon} point3 - Third reference point.
    * @param   {number} distance3 - Distance to third reference point (same units as radius).
    * @param   {number} [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
    * @returns {LatLon} Trilaterated point.
    *
    * @example
    *   LatLon.trilaterate(new LatLon(0, 0), 157e3, new LatLon(0, 1), 111e3, new LatLon(1, 0), 111e3); // 00.9985°N, 000.9986°E
    */
      static LatLonNvectorSpherical trilaterate(const LatLonNvectorSpherical& point1, double distance1,
                                                const LatLonNvectorSpherical& point2, double distance2,
                                                const LatLonNvectorSpherical& point3, double distance3,
                                                double radius = 6371e3);

       /**
      * Calculates the area of a spherical polygon where the sides of the polygon are great circle
      * arcs joining the vertices.
      *
      * Uses Girard’s theorem: A = [Σθᵢ − (n−2)·π]·R²
      *
      * @param   {LatLon[]} polygon - Array of points defining vertices of the polygon.
      * @param   {number}   [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
      * @returns {number}   The area of the polygon in the same units as radius.
      *
      * @example
      *   const polygon = [ new LatLon(0,0), new LatLon(1,0), new LatLon(0,1) ];
      *   const area = LatLon.areaOf(polygon); // 6.18e9 m²
      */
       static double areaOf(std::vector<LatLonNvectorSpherical>& polygon, double radius=6371e3);


      /**
      * Returns point representing geographic mean of supplied points.
      *
      * @param   {LatLon[]} points - Array of points to be averaged.
      * @returns {LatLon}   Point at the geographic mean of the supplied points.
      *
      * @example
      *   const p = LatLon.meanOf([ new LatLon(1, 1), new LatLon(4, 2), new LatLon(1, 3) ]); // 02.0001°N, 002.0000°E
      */
      static LatLonNvectorSpherical meanOf(const std::vector<LatLonNvectorSpherical>& points);

      /**
       * Returns whether this point is within the extent of a line segment joining point 1 & point 2.
       *
       * If this point is not on the great circle defined by point1 & point 2, returns whether it is
       * within the area bound by perpendiculars to the great circle at each point (in the same
       * hemisphere).
       *
       * @param   {LatLon}  point1 - First point defining segment.
       * @param   {LatLon}  point2 - Second point defining segment.
       * @returns {boolean} Whether this point is within extent of segment.
       *
       * @example
       *   const p1 = new LatLon(51, 1), p2 = new LatLon(52, 2);
       *   const within1 = new LatLon(52, 1).isWithinExtent(p1, p2); // true
       *   const within2 = new LatLon(51, 0).isWithinExtent(p1, p2); // false
       */
      [[nodiscard]] bool isWithinExtent(const LatLonNvectorSpherical& point1, const LatLonNvectorSpherical& point2) const;

       /**
        * Tests whether ‘this’ point is enclosed by the polygon defined by a set of points.
        *
        * @param   {LatLon[]} polygon - Ordered array of points defining vertices of polygon.
        * @returns {bool}     Whether this point is enclosed by polygon.
        *
        * @example
        *   const bounds = [ new LatLon(45,1), new LatLon(45,2), new LatLon(46,2), new LatLon(46,1) ];
        *   const p = new LatLon(45.1, 1.1);
        *   const inside = p.isEnclosedBy(bounds); // true
        */
      [[nodiscard]] bool isEnclosedBy(std::vector<LatLonNvectorSpherical>& polygon) const;
   };
}
#endif // LATLON_NVECTOR_SPHERICAL_H
