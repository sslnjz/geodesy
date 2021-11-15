//
// Created by Bin on 2021/11/15.
//

#include "latlon_ellipsoidal.h"

#include <sstream>
#include <iomanip>

using namespace geodesy;

LatLonEllipsoidal::LatLonEllipsoidal(double lat, double lon, double height)
    : _lat(lat)
    , _lon(lon)
    , _height(height)
{
}

std::wstring LatLonEllipsoidal::toString(int dp)
{
    return L"";
}
