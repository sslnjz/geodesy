//
// Created by Bin on 2021/11/15.
//

#ifndef LATLON_ELLIPSOIDAL_H
#define LATLON_ELLIPSOIDAL_H

#include <string>
#include <cmath>

#include "ellipsoids.h"
#include "vector3d.hpp"
#include "dms.h"

/**
 * A latitude/longitude point defines a geographic location on or above/below the earth’s surface,
 * measured in degrees from the equator & the International Reference Meridian and in metres above
 * the ellipsoid, and based on a given datum.
 *
 * As so much modern geodesy is based on WGS-84 (as used by GPS), this module includes WGS-84
 * ellipsoid parameters, and it has methods for converting geodetic (latitude/longitude) points to/from
 * geocentric cartesian points; the latlon-ellipsoidal-datum and latlon-ellipsoidal-referenceframe
 * modules provide transformation parameters for converting between historical datums and between
 * modern reference frames.
 *
 * This module is used for both trigonometric geodesy (eg latlon-ellipsoidal-vincenty) and n-vector
 * geodesy (eg latlon-nvector-ellipsoidal), and also for UTM/MGRS mapping.
 *
 * @module latlon-ellipsoidal
 */
namespace geodesy
{
    class Cartesian;
    /**
     * Latitude/longitude points on an ellipsoidal model earth, with ellipsoid parameters and methods
     * for converting points to/from cartesian (ECEF) coordinates.
     *
     * This is the core class, which will usually be used via LatLonEllipsoidal_Datum or
     * LatLonEllipsoidal_ReferenceFrame.
     */
    class LatLonEllipsoidal
    {
    public:
        /**
         * Creates a geodetic latitude/longitude point on a (WGS84) ellipsoidal model earth.
         *
         * @param  {number} lat - Latitude (in degrees).
         * @param  {number} lon - Longitude (in degrees).
         * @param  {number} [height=0] - Height above ellipsoid in metres.
         *
         * @example
         *   const auto p = new LatLonEllipsoidal(51.47788, -0.00147, 17);
         */
        LatLonEllipsoidal(double lat, double lon, double height = 0);

        /**
         * Latitude in degrees north from equator (including aliases lat, latitude): can be set as
         * numeric or hexagesimal (deg-min-sec); returned as numeric.
         */
        double lat()       { return _lat; }
        double latitude()  { return _lat; }
        void setLat(double lat) { _lat = Dms::wrap90(lat); }
        void setLatitude(double lat) { _lat = Dms::wrap90(lat); }

        /**
         * Longitude in degrees east from international reference meridian (including aliases lon, lng,
         * longitude): can be set as numeric or hexagesimal (deg-min-sec); returned as numeric.
         */
        double lon()       { return _lon; }
        double lng()       { return _lon; }
        double longitude() { return _lon; }
        void setLon(double lon) { _lon = Dms::wrap180(lon); }
        void setLng(double lon) { _lon = Dms::wrap180(lon); }
        void setLongitude(double lon) { _lon = Dms::wrap180(lon); }

        /**
         * Height in metres above ellipsoid.
         */
        double height() { return _height; }
        void setHeight(double height);

        /**
         * Datum.
         *
         * Note this is replicated within LatLonEllipsoidal in order that a LatLonEllipsoidal object can
         * be monkey-patched to look like a LatLonEllipsoidal_Datum, for Vincenty calculations on
         * different ellipsoids.
         *
         * @private
         */
        Datum datum() { return _datum; }
        void setDatum(const Datum& datum);


        /**
         * Ellipsoids with their parameters; this module only defines WGS84 parameters a = 6378137, b =
         * 6356752.314245, f = 1/298.257223563.
         *
         * @example
         *   const a = LatLon.ellipsoids.WGS84.a; // 6378137
         */
        static Ellipsoids ellipsoids()
        {
            return s_ellipsoids;
        }

        /**
         * Datums; this module only defines WGS84 datum, hence no datum transformations.
         *
         * @example
         *   const a = LatLon.datums.WGS84.ellipsoid.a; // 6377563.396
         */
        static Datums datums()
        {
            return s_datums;
        }

        /**
         * Converts ‘this’ point from (geodetic) latitude/longitude coordinates to (geocentric)
         * cartesian (x/y/z) coordinates.
         *
         * @returns {Cartesian} Cartesian point equivalent to lat/lon point, with x, y, z in metres from
         *   earth centre.
         */
        Cartesian toCartesian();


        /**
         * Returns a string representation of ‘this’ cartesian point.
         *
         * @param   {number} [dp=0] - Number of decimal places to use.
         * @returns {string} Comma-separated latitude/longitude.
         */
        std::wstring toString(int dp = 0);

    private:
        double _lat;
        double _lon;
        double _height;

        Datum _datum;
    };

    class Cartesian : public vector3d
    {
    public:
        Cartesian();
        Cartesian(double x, double y, double z);

        /**
         * Converts ‘this’ (geocentric) cartesian (x/y/z) coordinate to (geodetic) latitude/longitude
         * point on specified ellipsoid.
         *
         * Uses Bowring’s (1985) formulation for μm precision in concise form; ‘The accuracy of geodetic
         * latitude and height equations’, B R Bowring, Survey Review vol 28, 218, Oct 1985.
         *
         * @param   {LatLon.ellipsoids} [ellipsoid=WGS84] - Ellipsoid to use when converting point.
         * @returns {LatLon} Latitude/longitude point defined by cartesian coordinates, on given ellipsoid.
         * @throws  {TypeError} Invalid ellipsoid.
         *
         * @example
         *   const c = new Cartesian(4027893.924, 307041.993, 4919474.294);
         *   const p = c.toLatLon(); // 50.7978°N, 004.3592°E
         */
        LatLonEllipsoidal toLatLonEllipsoidal(Ellipsoid ellipsoid = s_ellipsoids.WGS84)
        {
            // note ellipsoid is available as a parameter for when toLatLon gets subclassed to
            // Ellipsoidal_Datum / Ellipsoidal_Referenceframe.
            const Cartesian cart = *this;
            const Ellipsoid ell = ellipsoid;

            const double e2 = 2 * ell.f - ell.f * ell.f;           // 1st eccentricity squared ≡ (a²−b²)/a²
            const double ε2 = e2 / (1-e2);         // 2nd eccentricity squared ≡ (a²−b²)/b²
            const double p = std::sqrt(cart.x()*cart.x() + cart.y()*cart.y()); // distance from minor axis
            const double R = std::sqrt(p*p + cart.z()*cart.z()); // polar radius

            // parametric latitude (Bowring eqn.17, replacing tanβ = z·a / p·b)
            const double tanβ = (ell.b * cart.z())/(ell.a * p) * (1+ ε2 * ell.b/R);
            const double sinβ = tanβ / std::sqrt(1+tanβ*tanβ);
            const double cosβ = sinβ / tanβ;

            // geodetic latitude (Bowring eqn.18: tanφ = z+ε²⋅b⋅sin³β / p−e²⋅cos³β)
            const double φ = std::isnan(cosβ) ? 0 : std::atan2(cart.z() + ε2*ell.b*sinβ*sinβ*sinβ, p - e2*ell.a*cosβ*cosβ*cosβ);

            // longitude
            const double λ = std::atan2(cart.y(), cart.x());

            // height above ellipsoid (Bowring eqn.7)
            const double sinφ = std::sin(φ), cosφ = std::cos(φ);
            const double ν = ell.a / std::sqrt(1-e2*sinφ*sinφ); // length of the normal terminated by the minor axis
            const double h = p*cosφ + cart.z()*sinφ - (ell.a*ell.a/ν);

            return {toDegrees(φ), toDegrees(λ), h};
        }


        /**
         * Returns a string representation of ‘this’ cartesian point.
         *
         * @param   {number} [dp=0] - Number of decimal places to use.
         * @returns {string} Comma-separated latitude/longitude.
         */
//        toString(dp=0) {
//            const x = this.x.toFixed(dp), y = this.y.toFixed(dp), z = this.z.toFixed(dp);
//            return `[${x},${y},${z}]`;
//        }


    };
}

#endif //GEODESY_LATLON_ELLIPSOIDAL_H
