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

#include "cartesian_datum.h"
#include "ellipsoids.h"
#include "latlon_ellipsoidal_datum.h"

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

   Transform inverseOf(const Transform& transform)
   {
      return {
         -transform.tx,
         -transform.ty,
         -transform.tz,
         -transform.s,
         -transform.rx,
         -transform.ry,
         -transform.rz
      };
   }
}

CartesianDatum::CartesianDatum(double x, double y, double z, const Datum& datum)
   : Cartesian(x, y, z)
   , m_datum(datum)
{
   validateDatum(m_datum);
}

Datum CartesianDatum::datum() const
{
   return m_datum;
}

void CartesianDatum::setDatum(const Datum& datum)
{
   validateDatum(datum);
   m_datum = datum;
}

LatLonEllipsoidalDatum CartesianDatum::toLatLon(std::optional<Datum> deprecatedDatum)
{
   const Datum datum = deprecatedDatum ? *deprecatedDatum : m_datum;
   validateDatum(datum);

   if (deprecatedDatum)
   {
      m_datum = datum;
   }

   const auto latLon = Cartesian::toLatLon(datum.ellipsoid);
   return LatLonEllipsoidalDatum(latLon.lat(), latLon.lon(), latLon.height(), datum);
}

CartesianDatum CartesianDatum::convertDatum(Datum toDatum) const
{
   validateDatum(toDatum);
   validateDatum(m_datum);

   if (m_datum == toDatum)
   {
      return CartesianDatum(x(), y(), z(), toDatum);
   }

   const CartesianDatum* oldCartesian = this;
   CartesianDatum wgs84Cartesian(0.0, 0.0, 0.0, g_datums.WGS84);
   Transform transform = toDatum.transforms;

   if (m_datum == g_datums.WGS84)
   {
      // Reference transforms are defined from WGS84 into the target datum.
      oldCartesian = this;
      transform = toDatum.transforms;
   }
   else if (toDatum == g_datums.WGS84)
   {
      // Inverting all seven Helmert parameters converts a source datum back to WGS84.
      oldCartesian = this;
      transform = inverseOf(m_datum.transforms);
   }
   else
   {
      // Non-WGS84 datum pairs are chained through WGS84, matching the JavaScript reference flow.
      wgs84Cartesian = convertDatum(g_datums.WGS84);
      oldCartesian = &wgs84Cartesian;
      transform = toDatum.transforms;
   }

   auto newCartesian = oldCartesian->applyTransform(transform);
   newCartesian.setDatum(toDatum);

   return newCartesian;
}

CartesianDatum CartesianDatum::applyTransform(Transform t) const
{
   const auto x1 = x(), y1 = y(), z1 = z();

   // Helmert transform parameters are stored as metres, parts-per-million, and arc-seconds.
   // Rotations and scale are normalised once here at the Cartesian calculation boundary.
   const auto tx = t.tx;
   const auto ty = t.ty;
   const auto tz = t.tz;
   const auto s  = t.s / 1e6 + 1.0;
   const auto rx = toRadians(t.rx / 3600.0);
   const auto ry = toRadians(t.ry / 3600.0);
   const auto rz = toRadians(t.rz / 3600.0);

   const auto x2 = tx + x1 * s  - y1 * rz + z1 * ry;
   const auto y2 = ty + x1 * rz + y1 * s  - z1 * rx;
   const auto z2 = tz - x1 * ry + y1 * rx + z1 * s;

   return CartesianDatum(x2, y2, z2, g_datums.WGS84);
}

std::string CartesianDatum::toString(int dp) const
{
   return vector3d::toString(dp);
}
