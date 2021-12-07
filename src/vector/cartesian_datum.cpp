/**********************************************************************************
*  MIT License                                                                    *
*                                                                                 *
*  Copyright (c) 2021 Binbin Song <ssln.jzs@gmail.com>                         *
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

#include <iostream>

using namespace geodesy;

CartesianDatum::CartesianDatum(double x, double y, double z, const Datum& datum)
   : Cartesian(x, y, z)
   , m_datum(datum)
{

}

Datum CartesianDatum::datum()
{
   return m_datum;
}

void CartesianDatum::setDatum(const Datum &datum)
{
   m_datum = datum;
}

LatLonEllipsoidalDatum CartesianDatum::toLatLon(std::optional<Datum> deprecatedDatum)
{
   if (deprecatedDatum)
   {
      std::cout << "datum parameter to Cartesian_Datum.toLatLon is deprecated: set datum before calling toLatLon()";
      m_datum = *deprecatedDatum;
   }
   const auto datum = m_datum ? m_datum : g_datums.WGS84;

   const auto latLon = Cartesian::toLatLon(datum.ellipsoid); // TODO: what if datum is not geocentric?
   const auto point = LatLonEllipsoidalDatum(latLon.lat(), latLon.lon(), latLon.height(), m_datum);
   return point;
}

CartesianDatum CartesianDatum::convertDatum(Datum toDatum)
{
   // TODO: what if datum is not geocentric?
   if (!toDatum || !toDatum.ellipsoid)
      throw std::invalid_argument("unrecognised datum");
   if (!m_datum)
      throw std::runtime_error("cartesian coordinate has no datum");

   CartesianDatum* oldCartesian;
   Transform transform;

   if (m_datum == g_datums.WGS84)
   {
      // converting from WGS 84
      oldCartesian = this;
      transform = toDatum.transforms;
   }
   if (toDatum == g_datums.WGS84) {
      // converting to WGS 84; use inverse transform
      oldCartesian = this;
      transform = m_datum.transforms.inverse();
   }
   if (!transform)
   {
      // neither this.datum nor toDatum are WGS84: convert this to WGS84 first
      CartesianDatum cas = convertDatum(g_datums.WGS84);
      oldCartesian = &cas;
      transform = toDatum.transforms;
   }

   auto newCartesian = oldCartesian->applyTransform(transform);
   newCartesian.m_datum = toDatum;

   return newCartesian;
}

CartesianDatum CartesianDatum::applyTransform(Transform t) const
{
   // this point
   const auto x1 = x(), y1 = y(), z1 = z();

   // transform parameters
   const auto tx = t.tx;                    // x-shift in metres
   const auto ty = t.ty;                    // y-shift in metres
   const auto tz = t.tz;                    // z-shift in metres
   const auto s  = t.s/1e6 + 1;            // scale: normalise parts-per-million to (s+1)
   const auto rx = toRadians((t.rx/3600)); // x-rotation: normalise arcseconds to radians
   const auto ry = toRadians((t.ry/3600)); // y-rotation: normalise arcseconds to radians
   const auto rz = toRadians((t.rz/3600)); // z-rotation: normalise arcseconds to radians

   // apply transform
   const auto x2 = tx + x1*s  - y1*rz + z1*ry;
   const auto y2 = ty + x1*rz + y1*s  - z1*rx;
   const auto z2 = tz - x1*ry + y1*rx + z1*s;

   return CartesianDatum(x2, y2, z2);
}


