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
#include "latlon_utm.h"

#include "algorithm.h"
#include "utm.h"

#include <array>
#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace geodesy;

namespace
{
constexpr double falseEastingMetres = 500e3;
constexpr double falseNorthingMetres = 10000e3;
constexpr double centralMeridianScale = 0.9996;
constexpr int krugerSeriesOrder = 6;

[[nodiscard]] int automaticZone(double longitudeDegrees)
{
   // Longitude 180 degrees belongs to zone 60, not a synthetic zone 61.
   return std::min(60, static_cast<int>(std::floor((longitudeDegrees + 180.0) / 6.0)) + 1);
}

[[nodiscard]] char latitudeBand(double latitudeDegrees)
{
   constexpr const char* mgrsLatBands = "CDEFGHJKLMNPQRSTUVWXX";
   const auto index = static_cast<size_t>(std::floor(latitudeDegrees / 8.0 + 10.0));
   return mgrsLatBands[index];
}

[[nodiscard]] double centralMeridianRadians(int zone)
{
   return toRadians((zone - 1) * 6.0 - 180.0 + 3.0);
}

void applyNorwayAndSvalbardExceptions(double latitudeDegrees, double longitudeDegrees,
   int& zone, double& lambda0)
{
   const char latBand = latitudeBand(latitudeDegrees);
   constexpr double zoneWidthRadians = pi / 30.0;

   // Norway and Svalbard use widened/narrowed zones to avoid narrow grid slivers.
   if (zone == 31 && latBand == 'V' && longitudeDegrees >= 3.0) { zone++; lambda0 += zoneWidthRadians; }
   if (zone == 32 && latBand == 'X' && longitudeDegrees < 9.0) { zone--; lambda0 -= zoneWidthRadians; }
   if (zone == 32 && latBand == 'X' && longitudeDegrees >= 9.0) { zone++; lambda0 += zoneWidthRadians; }
   if (zone == 34 && latBand == 'X' && longitudeDegrees < 21.0) { zone--; lambda0 -= zoneWidthRadians; }
   if (zone == 34 && latBand == 'X' && longitudeDegrees >= 21.0) { zone++; lambda0 += zoneWidthRadians; }
   if (zone == 36 && latBand == 'X' && longitudeDegrees < 33.0) { zone--; lambda0 -= zoneWidthRadians; }
   if (zone == 36 && latBand == 'X' && longitudeDegrees >= 33.0) { zone++; lambda0 += zoneWidthRadians; }
}
}

LatLonUtm::LatLonUtm(double lat, double lon, double height, std::optional<Datum> datum,
   std::optional<ReferenceFrame> reference, std::optional<std::string> epoch)
   : LatLonEllipsoidal(lat, lon, height, datum, reference, epoch)
   , m_convergence(0.0)
   , m_scale(0.0)
{
}

Utm LatLonUtm::toUtm(std::optional<int> zoneOverride) const
{
   if (!(-80.0 <= lat() && lat() <= 84.0))
      throw std::invalid_argument("latitude outside UTM limits");

   if (zoneOverride.has_value() && !(*zoneOverride >= 1 && *zoneOverride <= 60))
      throw std::invalid_argument("invalid UTM zone override: expected integer in range [1, 60]");

   auto zone = zoneOverride ? *zoneOverride : automaticZone(lon());
   auto lambda0 = centralMeridianRadians(zone);

   if (!zoneOverride.has_value())
      applyNorwayAndSvalbardExceptions(lat(), lon(), zone, lambda0);

   const auto phi = toRadians(lat());
   const auto lambda = toRadians(lon()) - lambda0;

   const Datum datum = m_datum.value_or(datums().WGS84);
   const auto ellipsoid = datum.ellipsoid;
   const auto a = ellipsoid.a, f = ellipsoid.f;

   // Karney 2011 equations 7-14 project geodetic latitude/longitude to the conformal sphere.
   const auto e = std::sqrt(f * (2.0 - f));
   const auto n = f / (2.0 - f);
   const auto n2 = n * n, n3 = n * n2, n4 = n * n3, n5 = n * n4, n6 = n * n5;

   const auto coslambda = std::cos(lambda), sinlambda = std::sin(lambda), tanlambda = std::tan(lambda);
   const auto tau = std::tan(phi);
   const auto sigma = std::sinh(e * std::atanh(e * tau / std::sqrt(1.0 + tau * tau)));
   const auto taud = tau * std::sqrt(1.0 + sigma * sigma) - sigma * std::sqrt(1.0 + tau * tau);

   const auto ksid = std::atan2(taud, coslambda);
   const auto etad = std::asinh(sinlambda / std::sqrt(taud * taud + coslambda * coslambda));

   const auto A = a / (1.0 + n) * (1.0 + 1.0/4.0 * n2 + 1.0/64.0 * n4 + 1.0/256.0 * n6);

   const std::array<double, krugerSeriesOrder> alpha = {
      1.0/2.0 * n - 2.0/3.0 * n2 + 5.0/16.0 * n3 + 41.0/180.0 * n4 - 127.0/288.0 * n5 + 7891.0/37800.0 * n6,
                    13.0/48.0 * n2 - 3.0/5.0 * n3 + 557.0/1440.0 * n4 + 281.0/630.0 * n5 - 1983433.0/1935360.0 * n6,
                                 61.0/240.0 * n3 - 103.0/140.0 * n4 + 15061.0/26880.0 * n5 + 167603.0/181440.0 * n6,
                                             49561.0/161280.0 * n4 - 179.0/168.0 * n5 + 6601661.0/7257600.0 * n6,
                                                                  34729.0/80640.0 * n5 - 3418889.0/1995840.0 * n6,
                                                                                      212378941.0/319334400.0 * n6
   };

   double epsilon = ksid;
   double eta = etad;
   for (int j = 1; j <= krugerSeriesOrder; j++)
   {
      // The coefficient tables are stored zero-based, but the Kruger series is one-based.
      epsilon += alpha[j - 1] * std::sin(2.0 * j * ksid) * std::cosh(2.0 * j * etad);
      eta += alpha[j - 1] * std::cos(2.0 * j * ksid) * std::sinh(2.0 * j * etad);
   }

   auto x = centralMeridianScale * A * eta;
   auto y = centralMeridianScale * A * epsilon;

   // Karney 2011 equations 23-25 provide grid convergence and scale from the same projection terms.
   auto pd = 1.0;
   auto qd = 0.0;
   for (int j = 1; j <= krugerSeriesOrder; j++)
   {
      pd += 2.0 * j * alpha[j - 1] * std::cos(2.0 * j * ksid) * std::cosh(2.0 * j * etad);
      qd += 2.0 * j * alpha[j - 1] * std::sin(2.0 * j * ksid) * std::sinh(2.0 * j * etad);
   }

   const auto gammad = std::atan(taud / std::sqrt(1.0 + taud * taud) * tanlambda);
   const auto gammasd = std::atan2(qd, pd);
   const auto gamma = gammad + gammasd;

   const auto sinphi = std::sin(phi);
   const auto kappad = std::sqrt(1.0 - e * e * sinphi * sinphi) * std::sqrt(1.0 + tau * tau)
      / std::sqrt(taud * taud + coslambda * coslambda);
   const auto kappasd = A / a * std::sqrt(pd * pd + qd * qd);
   const auto kappa = centralMeridianScale * kappad * kappasd;

   x += falseEastingMetres;
   if (y < 0.0)
      y += falseNorthingMetres;

   const auto easting = std::stod(geodesy::toFixed(x, 9));
   const auto northing = std::stod(geodesy::toFixed(y, 9));
   const auto convergence = std::stod(geodesy::toFixed(geodesy::toDegrees(gamma), 9));
   const auto scale = std::stod(geodesy::toFixed(kappa, 12));
   const auto hemisphere = lat() >= 0.0 ? Utm::Hemisphere::N : Utm::Hemisphere::S;

   return Utm(zone, hemisphere, easting, northing, datum, convergence, scale, zoneOverride.has_value());
}

void LatLonUtm::setConvergence(double convergence)
{
   if (!std::isfinite(convergence))
      throw std::invalid_argument("invalid UTM convergence: expected finite degrees");

   m_convergence = convergence;
}

void LatLonUtm::setScale(double scale)
{
   if (!std::isfinite(scale) || scale <= 0.0)
      throw std::invalid_argument("invalid UTM scale: expected positive finite scale factor");

   m_scale = scale;
}
