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

#ifndef LATLON_ELLIPSOIDAL_VINCENTY_H
#define LATLON_ELLIPSOIDAL_VINCENTY_H

/**
 * Distances & bearings between points, and destination points given start points & initial bearings,
 * calculated on an ellipsoidal earth model using ‘direct and inverse solutions of geodesics on the
 * ellipsoid’ devised by Thaddeus VincentyInverse.
 *
 * From: T VincentyInverse, "Direct and Inverse Solutions of Geodesics on the Ellipsoid with application of
 * nested equations", Survey Review, vol XXIII no 176, 1975. www.ngs.noaa.gov/PUBS_LIB/inverse.pdf.
 *
 * @module latlon-ellipsoidal-vincenty
 */

#include "latlon_ellipsoidal.h"
#include "algorithm.h"

namespace geodesy
{
    /**
     * Extends LatLonEllipsoidal with methods for calculating distances and bearings between points, and
     * destination points given distances and initial bearings, accurate to within 0.5mm distance,
     * 0.000015″ bearing.
     *
     * By default, these calculations are made on a WGS-84 ellipsoid. For geodesic calculations on other
     * ellipsoids, monkey-patch the LatLon point by setting the datum of ‘this’ point to make it appear
     * as a LatLonEllipsoidal_Datum or LatLonEllipsoidal_ReferenceFrame point: e.g.
     *
     *     import LatLon, { Dms } from '../latlon-ellipsoidal-vincenty.js';
     *     import { datums }      from '../latlon-ellipsoidal-datum.js';
     *     const le = new LatLon(50.065716, -5.713824);  // in OSGB-36
     *     const jog = new LatLon(58.644399, -3.068521); // in OSGB-36
     *     le.datum = datums.OSGB36;     // source point determines ellipsoid to use
     *     const d = le.distanceTo(jog); // = 969982.014; 27.848m more than on WGS-84 ellipsoid
     *
     * @extends LatLonEllipsoidal
     */
     struct VI;
     struct VD;
     class LatLonEllipsoidalVincenty : public LatLonEllipsoidal
     {
     public:
        LatLonEllipsoidalVincenty();
        LatLonEllipsoidalVincenty(double lat, double lon, double height = 0.0,
                                  std::optional<Datum> datum = std::nullopt,
                                  std::optional<ReferenceFrame> reference = std::nullopt,
                                  std::optional<std::string> epoch = std::nullopt);

        /**
         * Returns the distance between ‘this’ point and destination point along a geodesic on the
         * surface of the ellipsoid, using Vincenty inverse solution.
         *
         * @param   {LatLon} point - Latitude/longitude of destination point.
         * @returns {number} Distance in metres between points or NaN if failed to converge.
         *
         * @example
         *   const p1 = new LatLon(50.06632, -5.71475);
         *   const p2 = new LatLon(58.64402, -3.07009);
         *   const d = p1.distanceTo(p2); // 969,954.166 m
         */
        double distanceTo(const LatLonEllipsoidal& point) const;

        /**
         * Returns the initial bearing (forward azimuth) to travel along a geodesic from ‘this’ point to
         * the given point, using Vincenty inverse solution.
         *
         * @param   {LatLon} point - Latitude/longitude of destination point.
         * @returns {number} Initial bearing in degrees from north (0°..360°) or NaN if failed to converge.
         *
         * @example
         *   const p1 = new LatLon(50.06632, -5.71475);
         *   const p2 = new LatLon(58.64402, -3.07009);
         *   const b1 = p1.initialBearingTo(p2); // 9.1419°
         */
        double initialBearingTo(const LatLonEllipsoidal& point) const;

        /**
         * Returns the final bearing (reverse azimuth) having travelled along a geodesic from ‘this’
         * point to the given point, using Vincenty inverse solution.
         *
         * @param   {LatLon} point - Latitude/longitude of destination point.
         * @returns {number} Final bearing in degrees from north (0°..360°) or NaN if failed to converge.
         *
         * @example
         *   const p1 = new LatLon(50.06632, -5.71475);
         *   const p2 = new LatLon(58.64402, -3.07009);
         *   const b2 = p1.finalBearingTo(p2); // 11.2972°
         */
        double finalBearingTo(const LatLonEllipsoidal& point) const;

        /**
         * Returns the destination point having travelled the given distance along a geodesic given by
         * initial bearing from ‘this’ point, using Vincenty direct solution.
         *
         * @param   {number} distance - Distance travelled along the geodesic in metres.
         * @param   {number} initialBearing - Initial bearing in degrees from north.
         * @returns {LatLon} Destination point.
         *
         * @example
         *   const p1 = new LatLon(-37.95103, 144.42487);
         *   const p2 = p1.destinationPoint(54972.271, 306.86816); // 37.6528°S, 143.9265°E
         */
        LatLonEllipsoidalVincenty destinationPoint(double distance, double initialBearing) const;

         /**
          * Returns the final bearing (reverse azimuth) having travelled along a geodesic given by initial
          * bearing for a given distance from ‘this’ point, using Vincenty direct solution.
          * TODO: arg order? (this is consistent with destinationPoint, but perhaps less intuitive)
          *
          * @param   {number} distance - Distance travelled along the geodesic in metres.
          * @param   {LatLon} initialBearing - Initial bearing in degrees from north.
          * @returns {number} Final bearing in degrees from north (0°..360°).
          *
          * @example
          *   const p1 = new LatLon(-37.95103, 144.42487);
          *   const b2 = p1.finalBearingOn(54972.271, 306.86816); // 307.1736°
          */
         double finalBearingOn(double distance, double initialBearing) const;


         /**
          * Returns the point at given fraction between ‘this’ point and given point.
          *
          * @param   {LatLon} point - Latitude/longitude of destination point.
          * @param   {number} fraction - Fraction between the two points (0 = this point, 1 = specified point).
          * @returns {LatLon} Intermediate point between this point and destination point.
          *
          * @example
          *   const p1 = new LatLon(50.06632, -5.71475);
          *   const p2 = new LatLon(58.64402, -3.07009);
          *   const pInt = p1.intermediatePointTo(p2, 0.5); // 54.3639°N, 004.5304°W
          */
         LatLonEllipsoidalVincenty intermediatePointTo( const LatLonEllipsoidal& point, double fraction) const;

     private:
        /**
         * Vincenty direct calculation.
         *
         * Ellipsoid parameters are taken from datum of 'this' point. Height is ignored.
         *
         * @private
         * @param   {number} distance - Distance along bearing in metres.
         * @param   {number} initialBearing - Initial bearing in degrees from north.
         * @returns (Object} Object including point (destination point), finalBearing.
         * @throws  {RangeError} Point must be on surface of ellipsoid.
         * @throws  {EvalError}  Formula failed to converge.
         */
        VD direct(double distance, double initialBearing) const;

        /**
         * Vincenty inverse calculation.
         *
         * Ellipsoid parameters are taken from datum of 'this' point. Height is ignored.
         *
         * @private
         * @param   {LatLon} point - Latitude/longitude of destination point.
         * @returns {Object} Object including distance, initialBearing, finalBearing.
         * @throws  {TypeError}  Invalid point.
         * @throws  {RangeError} Points must be on surface of ellipsoid.
         * @throws  {EvalError}  Formula failed to converge.
         */
        VI inverse(const LatLonEllipsoidal& point) const;
    };
}

#endif // LATLON_ELLIPSOIDAL_VINCENTY_H
