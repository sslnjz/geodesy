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

namespace geodesy {
   /**
    * Latitude/longitude points on an ellipsoidal model earth, with ellipsoid parameters and methods
    * for converting between reference frames and to geocentric (ECEF) cartesian coordinates.
    *
    * @extends LatLonEllipsoidal
    */
   class CartesianReferenceFrame;

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
      LatLonEllipsoidalReferenceFrame(double lat, double lon, double height = 0,
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
     * Parses a latitude/longitude point from a variety of formats.
     *
     * Latitude & longitude (in degrees) can be supplied as two separate parameters, as a single
     * comma-separated lat/lon string, or as a single object with { lat, lon } or GeoJSON properties.
     *
     * The latitude/longitude values may be numeric or strings; they may be signed decimal or
     * deg-min-sec (hexagesimal) suffixed by compass direction (NSEW); a variety of separators are
     * accepted. Examples -3.62, '3 37 12W', '3°37′12″W'.
     *
     * Thousands/decimal separators must be comma/dot; use Dms.fromLocale to convert locale-specific
     * thousands/decimal separators.
     *
     * @param   {number|string|Object}   lat|latlon - Geodetic Latitude (in degrees) or comma-separated lat/lon or lat/lon object.
     * @param   {number}                 [lon] - Longitude in degrees.
     * @param   {number}                 height - Height above ellipsoid in metres.
     * @param   {LatLon.referenceFrames} referenceFrame - Reference frame this point is defined within.
     * @param   {number} [epoch=referenceFrame.epoch] - date of observation of coordinate (decimal year).
     * @returns {LatLon} Latitude/longitude point on ellipsoidal model earth using given reference frame.
     * @throws {TypeError} Unrecognised reference frame.
     *
     * @example
     *   const p1 = LatLon.parse(51.47788, -0.00147, 17, LatLon.referenceFrames.ETRF2000);          // numeric pair
     *   const p2 = LatLon.parse('51°28′40″N, 000°00′05″W', 17, LatLon.referenceFrames.ETRF2000);   // dms string + height
     *   const p3 = LatLon.parse({ lat: 52.205, lon: 0.119 }, 17, LatLon.referenceFrames.ETRF2000); // { lat, lon } object numeric
     */
      static LatLonEllipsoidalReferenceFrame parse(double lat, double lon, double height,
                                                   std::optional<ReferenceFrame> referenceFrame = std::nullopt,
                                                   std::optional<std::string> epoch = std::nullopt);
      static LatLonEllipsoidalReferenceFrame parse(const std::string& dms, double height,
                                                   std::optional<ReferenceFrame> referenceFrame = std::nullopt,
                                                   std::optional<std::string> epoch = std::nullopt);
      static LatLonEllipsoidalReferenceFrame parse(const std::string& lat, const std::string& lon, double height,
                                                   std::optional<ReferenceFrame> referenceFrame = std::nullopt,
                                                   std::optional<std::string> epoch = std::nullopt);
      static LatLonEllipsoidalReferenceFrame parse(const std::string& lat, const std::string& lon, std::string height,
                                                   std::optional<ReferenceFrame> referenceFrame = std::nullopt,
                                                   std::optional<std::string> epoch = std::nullopt);


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
      static std::vector<HelmertTransforms> transformParameters();

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
      LatLonEllipsoidalReferenceFrame convertReferenceFrame(const ReferenceFrame &to) const;

      /**
      * Converts ‘this’ point from (geodetic) latitude/longitude coordinates to (geocentric) cartesian
      * (x/y/z) coordinates, based on same reference frame.
      *
      * Shadow of LatLonEllipsoidal.toCartesian(), returning Cartesian augmented with
      * LatLonEllipsoidal_ReferenceFrame methods/properties.
      *
      * @returns {Cartesian} Cartesian point equivalent to lat/lon point, with x, y, z in metres from
      *   earth centre, augmented with reference frame conversion methods and properties.
      */
      CartesianReferenceFrame toCartesian() const;

      /**
     * Returns a string representation of ‘this’ point, formatted as degrees, degrees+minutes, or
     * degrees+minutes+seconds.
     *
     * @param   {string} [format=d] - Format point as 'd', 'dm', 'dms'.
     * @param   {number} [dp=4|2|0] - Number of decimal places to use: default 4 for d, 2 for dm, 0 for dms.
     * @param   {number} [dpHeight=null] - Number of decimal places to use for height; default (null) is no height display.
     * @param   {boolean} [referenceFrame=false] - Whether to show reference frame point is defined on.
     * @returns {string} Comma-separated formatted latitude/longitude.
     *
     * @example
     *   new LatLon(51.47788, -0.00147, 0, LatLon.referenceFrames.ITRF2014).toString();             // 51.4778°N, 000.0015°W
     *   new LatLon(51.47788, -0.00147, 0, LatLon.referenceFrames.ITRF2014).toString('dms');        // 51°28′40″N, 000°00′05″W
     *   new LatLon(51.47788, -0.00147, 42, LatLon.referenceFrames.ITRF2014).toString('dms', 0, 0); // 51°28′40″N, 000°00′05″W +42m
     */
      std::string toString( Dms::eFormat format=Dms::D, std::optional<int> dp = std::nullopt,
                            std::optional<int> dph = std::nullopt, bool referenceFrame=false)
      {
         const auto ll = LatLonEllipsoidal::toString(format, dp, dph);
         const auto epoch = m_referenceFrame && m_epoch != m_referenceFrame->epoch ? m_epoch : "";

         const auto trf = referenceFrame ? "(" + m_referenceFrame->name + (epoch ? ("@" + *epoch) : " ") : "";
         return ll + trf;
      }
   };

}

#endif //LATLON_ELLIPSOIDAL_REFERENCEFRAME_H
