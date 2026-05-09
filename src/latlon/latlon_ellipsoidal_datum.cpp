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

#include "latlon_ellipsoidal_datum.h"
#include "cartesian_datum.h"

#include <cmath>
#include <stdexcept>

using namespace geodesy;

namespace
{
   bool isValidEllipsoid(const Ellipsoid& ellipsoid)
   {
      return std::isfinite(ellipsoid.a) && std::isfinite(ellipsoid.b) && std::isfinite(ellipsoid.f)
         && ellipsoid.a > 1.0 && ellipsoid.b > 1.0 && ellipsoid.b <= ellipsoid.a
         && ellipsoid.f > 0.0 && ellipsoid.f < 1.0;
   }

   void validateDatum(const Datum& datum)
   {
      if (!isValidEllipsoid(datum.ellipsoid))
      {
         throw std::invalid_argument("unrecognised datum");
      }
   }
}

LatLonEllipsoidalDatum::LatLonEllipsoidalDatum(double lat, double lon, double height, std::optional<Datum> datum)
   : LatLonEllipsoidal(lat, lon, height, datum ? datum : std::optional<Datum>{g_datums.WGS84})
{
   validateDatum(*m_datum);
}

Datum LatLonEllipsoidalDatum::datum() const
{
   return *m_datum;
}

Ellipsoids LatLonEllipsoidalDatum::ellipsoids()
{
   return g_ellipsoids;
}

Datums LatLonEllipsoidalDatum::datums()
{
   return g_datums;
}

LatLonEllipsoidalDatum LatLonEllipsoidalDatum::parse(double lat, double lon, double height, Datum datum)
{
   validateDatum(datum);
   const auto latlon = LatLonEllipsoidal::parse(lat, lon, height);
   return LatLonEllipsoidalDatum(latlon.lat(), latlon.lon(), latlon.height(), datum);
}

LatLonEllipsoidalDatum LatLonEllipsoidalDatum::parse(const std::string& lat,
                                                     const std::string& lon, double height, Datum datum)
{
   validateDatum(datum);
   const auto latlon = LatLonEllipsoidal::parse(lat, lon, height);
   return LatLonEllipsoidalDatum(latlon.lat(), latlon.lon(), latlon.height(), datum);
}

LatLonEllipsoidalDatum LatLonEllipsoidalDatum::parse(const std::string& dms, double height, Datum datum)
{
   validateDatum(datum);
   const auto latlon = LatLonEllipsoidal::parse(dms, height);
   return LatLonEllipsoidalDatum(latlon.lat(), latlon.lon(), latlon.height(), datum);
}

LatLonEllipsoidalDatum LatLonEllipsoidalDatum::parse(const std::string& lat, const std::string& lon,
                                                     std::string height, Datum datum)
{
   validateDatum(datum);
   const auto latlon = LatLonEllipsoidal::parse(lat, lon, height);
   return LatLonEllipsoidalDatum(latlon.lat(), latlon.lon(), latlon.height(), datum);
}

LatLonEllipsoidalDatum LatLonEllipsoidalDatum::convertDatum(Datum toDatum) const
{
   validateDatum(toDatum);
   validateDatum(*m_datum);

   // Convert through ECEF metres so the datum transform is applied once in Cartesian space.
   return toCartesian().convertDatum(toDatum).toLatLon();
}

CartesianDatum LatLonEllipsoidalDatum::toCartesian() const
{
   const auto cartesian = LatLonEllipsoidal::toCartesian();
   return CartesianDatum(cartesian.x(), cartesian.y(), cartesian.z(), datum());
}
