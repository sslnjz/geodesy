//
// Created by Bin on 2021/11/15.
//

#ifndef LATLON_ELLIPSOIDAL_H
#define LATLON_ELLIPSOIDAL_H

#include <string>
#include <cmath>

#include "ellipsoids.h"
#include "dms.h"

namespace geodesy
{
    class LatLonEllipsoidal
    {
    public:
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
        double height()       { return _height; }
        void setHeight(double height) {
            if (std::isnan(_height)) throw std::invalid_argument("invalid height");
            _height = height;
        }


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
    };
}

#endif //GEODESY_LATLON_ELLIPSOIDAL_H
