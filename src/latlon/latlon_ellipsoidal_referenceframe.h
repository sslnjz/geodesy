/**********************************************************************************
*  MIT License                                                                    *
*                                                                                 *
*  Copyright (c) 2021 Binbin Song <ssln.jzs@gmail.com>                       *
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

#ifndef LATLON_ELLIPSOIDAL_REFERENCEFRAME_H
#define LATLON_ELLIPSOIDAL_REFERENCEFRAME_H

#include "ellipsoids.h"
#include "latlon_ellipsoidal_referenceframe_txparams.h"
#include "latlon_ellipsoidal.h"

/**
 * Modern geodetic reference frames: a latitude/longitude point defines a geographic location on or
 * above/below the earth’s surface, measured in degrees from the equator and the International
 * Reference Meridian and metres above the ellipsoid within a given terrestrial reference frame at a
 * given epoch.
 *
 * This module extends the core latlon-ellipsoidal module to include methods for converting between
 * different reference frames.
 *
 * This is scratching the surface of complexities involved in high precision geodesy, but may be of
 * interest and/or value to those with less demanding requirements.
 *
 * Note that ITRF solutions do not directly use an ellipsoid, but are specified by cartesian
 * coordinates; the GRS80 ellipsoid is recommended for transformations to geographical coordinates
 * (itrf.ensg.ign.fr).
 *
 * @module latlon-ellipsoidal-referenceframe
 */


/**
 * Sources:
 *
 * - Soler & Snay, “Transforming Positions and Velocities between the International Terrestrial Refer-
 *   ence Frame of 2000 and North American Datum of 1983”, Journal of Surveying Engineering May 2004;
 *   www.ngs.noaa.gov/CORS/Articles/SolerSnayASCE.pdf.
 *
 * - Dawson & Woods, “ITRF to GDA94 coordinate transformations”, Journal of Applied Geodesy 4 (2010);
 *   www.ga.gov.au/webtemp/image_cache/GA19050.pdf.
 */

namespace geodesy
{
    /**
     * Latitude/longitude points on an ellipsoidal model earth, with ellipsoid parameters and methods
     * for converting between reference frames and to geocentric (ECEF) cartesian coordinates.
     *
     * @extends LatLonEllipsoidal
     */
    class LatLonEllipsoidalReferenceFrame : public LatLonEllipsoidal
    {
    public:
        /**
         * Creates geodetic latitude/longitude point on an ellipsoidal model earth using using a
         * specified reference frame.
         *
         * Note that while the epoch defaults to the frame reference epoch, the accuracy of ITRF
         * realisations is meaningless without knowing the observation epoch.
         *
         * @param  {number} lat - Geodetic latitude in degrees.
         * @param  {number} lon - Geodetic longitude in degrees.
         * @param  {number} [height=0] - Height above ellipsoid in metres.
         * @param  {LatLon.referenceFrames} [referenceFrame=ITRF2014] - Reference frame this point is defined within.
         * @param  {number} [epoch=referenceFrame.epoch] - date of observation of coordinate (decimal year).
         *   defaults to reference epoch t₀ of reference frame.
         * @throws {TypeError} Unrecognised reference frame.
         *
         * @example
         *   import LatLon from '/js/geodesy/latlon-ellipsoidal-referenceframe.js';
         *   const p = new LatLon(51.47788, -0.00147, 0, LatLon.referenceFrames.ITRF2000);
         */
        LatLonEllipsoidalReferenceFrame(double lat, double lon, double height=0,
                                        std::optional<ReferenceFrame> referenceFrame = g_reference_frames.ITRF2014,
                                        std::optional<std::string> epoch = std::nullopt);

        /**
         * Reference frame this point is defined within.
         */
        [[nodiscard]] std::optional<ReferenceFrame> referenceFrame() const;


        /**
         * Point’s observed epoch.
         */
        std::optional<std::string> epoch();


        /**
         * Ellipsoid parameters; semi-major axis (a), semi-minor axis (b), and flattening (f).
         *
         * The only ellipsoids used in modern geodesy are WGS-84 and GRS-80 (while based on differing
         * defining parameters, the only effective difference is a 0.1mm variation in the minor axis b).
         *
         * @example
         *   const availableEllipsoids = Object.keys(LatLon.ellipsoids).join(); // WGS84,GRS80
         *   const a = LatLon.ellipsoids.Airy1830.a;                            // 6377563.396
         */
        static Ellipsoids ellipsoids();


        /**
         * Reference frames, with their base ellipsoids and reference epochs.
         *
         * @example
         *   const availableReferenceFrames = Object.keys(LatLon.referenceFrames).join(); // ITRF2014,ITRF2008, ...
         */
        static ReferenceFrames referenceFrames();

        /**
         * 14-parameter Helmert transformation parameters between (dynamic) ITRS frames, and from ITRS
         * frames to (static) regional TRFs NAD83, ETRF2000, and GDA94.
         *
         * This is a limited set of transformations; e.g. ITRF frames prior to ITRF2000 are not included.
         * More transformations could be added on request.
         *
         * Many conversions are direct; for NAD83, successive ITRF transformations are chained back to
         * ITRF2000.
         */
        static std::map<std::string, Helmert> transformParameters();

        /**
         * Converts ‘this’ lat/lon coordinate to new coordinate system.
         *
         * @param   {LatLon.referenceFrames} toReferenceFrame - Reference frame this coordinate is to be converted to.
         * @returns {LatLon} This point converted to new reference frame.
         * @throws  {Error}  Undefined reference frame, Transformation not available.
         *
         * @example
         *   const pEtrf = new LatLon(51.47788000, -0.00147000, 0, LatLon.referenceFrames.ITRF2000);
         *   const pItrf = pEtrf.convertReferenceFrame(LatLon.referenceFrames.ETRF2000); // 51.47787826°N, 000.00147125°W
         */
//        convertReferenceFrame(toReferenceFrame) {
//                if (!toReferenceFrame || toReferenceFrame.epoch == undefined) throw new TypeError('unrecognised reference frame');
//
//                const oldCartesian = this.toCartesian();                                   // convert geodetic to cartesian
//                const newCartesian = oldCartesian.convertReferenceFrame(toReferenceFrame); // convert TRF
//                const newLatLon = newCartesian.toLatLon();                                 // convert cartesian back to to geodetic
//
//                return newLatLon;
//        }
    };

}

#endif //LATLON_ELLIPSOIDAL_REFERENCEFRAME_H
