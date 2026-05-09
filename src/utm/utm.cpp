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
#include "utm.h"

#include "algorithm.h"
#include "latlon_utm.h"
#include "strutil.h"

#include <array>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

using namespace geodesy;

namespace
{
constexpr double minEastingMetres = 0.0;
constexpr double maxEastingMetres = 1000e3;
constexpr double maxNorthingNorthernMetres = 9329006.0;
constexpr double minNorthingSouthernMetres = 1116914.0;
constexpr double falseEastingMetres = 500e3;
constexpr double falseNorthingMetres = 10000e3;
constexpr double centralMeridianScale = 0.9996;
constexpr double tauConvergenceTolerance = 1e-12;
constexpr int maxTauIterations = 10;
constexpr int krugerSeriesOrder = 6;

[[nodiscard]] bool isFinite(double value)
{
   return std::isfinite(value);
}

[[nodiscard]] std::string hemisphereToString(Utm::Hemisphere hemisphere)
{
   switch (hemisphere)
   {
   case Utm::Hemisphere::N:
      return "N";
   case Utm::Hemisphere::S:
      return "S";
   }

   throw std::invalid_argument("invalid UTM hemisphere");
}

[[nodiscard]] Utm::Hemisphere parseHemisphere(const std::string& hemisphere)
{
   if (hemisphere == "N" || hemisphere == "n")
      return Utm::Hemisphere::N;
   if (hemisphere == "S" || hemisphere == "s")
      return Utm::Hemisphere::S;

   throw std::invalid_argument("invalid UTM hemisphere: expected N or S");
}
}

Utm::Utm(int zone, Hemisphere h, double easting, double northing, std::optional<Datum> datum,
                  std::optional<double> convergence, std::optional<double> scale, bool verifyEN)
   : m_zone(zone), m_hemisphere(h), m_easting(easting), m_northing(northing)
   , m_datum(datum), m_convergence(convergence), m_scale(scale)
{
   if (!(1 <= zone && zone <= 60))
      throw std::invalid_argument("invalid UTM zone: expected integer in range [1, 60]");

   (void)hemisphereToString(h);

   if (!isFinite(easting))
      throw std::invalid_argument("invalid UTM easting: expected finite metres");
   if (!isFinite(northing))
      throw std::invalid_argument("invalid UTM northing: expected finite metres");
   if (convergence.has_value() && !isFinite(*convergence))
      throw std::invalid_argument("invalid UTM convergence: expected finite degrees");
   if (scale.has_value() && (!isFinite(*scale) || *scale <= 0.0))
      throw std::invalid_argument("invalid UTM scale: expected positive finite scale factor");

   if (verifyEN)
   {
      // UTM eastings and northings are metres from false origins; these rough validity
      // ranges match the reference library while allowing opt-out for extended coordinates.
      if (!(minEastingMetres <= easting && easting <= maxEastingMetres))
         throw std::invalid_argument("invalid UTM easting: out of range [0, 1000000]");

      switch (h)
      {
      case Hemisphere::N:
         if (!(0.0 <= northing && northing < maxNorthingNorthernMetres))
            throw std::invalid_argument("invalid UTM northing: out of range [0, 9329006)");
         break;
      case Hemisphere::S:
         if (!(minNorthingSouthernMetres <= northing && northing < falseNorthingMetres))
            throw std::invalid_argument("invalid UTM northing: out of range [1116914, 10000000)");
         break;
      }
   }

   if (!datum.has_value())
      throw std::invalid_argument("unrecognized datum: empty");
}

Utm Utm::parse(const std::string& utmCoord, std::optional<Datum> datum)
{
   // A valid UTM text coordinate is four whitespace-separated fields:
   // zone, hemisphere, easting in metres, and northing in metres.
   std::istringstream stream(strutil::strip(utmCoord));
   std::string zoneText;
   std::string hemisphereText;
   std::string eastingText;
   std::string northingText;
   std::string extra;

   if (!(stream >> zoneText >> hemisphereText >> eastingText >> northingText) || (stream >> extra))
      throw std::invalid_argument("invalid UTM coordinate: expected zone hemisphere easting northing");

   try
   {
      std::size_t parsed = 0;
      const int zone = std::stoi(zoneText, &parsed);
      if (parsed != zoneText.size())
         throw std::invalid_argument("invalid UTM zone");

      parsed = 0;
      const double easting = std::stod(eastingText, &parsed);
      if (parsed != eastingText.size())
         throw std::invalid_argument("invalid UTM easting");

      parsed = 0;
      const double northing = std::stod(northingText, &parsed);
      if (parsed != northingText.size())
         throw std::invalid_argument("invalid UTM northing");

      return Utm(zone, parseHemisphere(hemisphereText), easting, northing, datum);
   }
   catch (const std::invalid_argument&)
   {
      throw;
   }
   catch (const std::out_of_range&)
   {
      throw std::invalid_argument("invalid UTM coordinate: numeric field out of range");
   }
}

LatLonUtm Utm::toLatLon() const
{
   const auto z = m_zone;
   const auto h = m_hemisphere;
   const auto a = m_datum->ellipsoid.a, f = m_datum->ellipsoid.f;

   const auto x = m_easting - falseEastingMetres;
   const auto y = h == Hemisphere::S ? m_northing - falseNorthingMetres : m_northing;

   // Karney 2011 equations 15-22 invert the conformal-sphere projection back to geodetic latitude.
   const auto e = std::sqrt(f * (2.0 - f));
   const auto n = f / (2.0 - f);
   const auto n2 = n * n, n3 = n * n2, n4 = n * n3, n5 = n * n4, n6 = n * n5;

   const auto A = a / (1.0 + n) * (1.0 + 1.0/4.0 * n2 + 1.0/64.0 * n4 + 1.0/256.0 * n6);

   const auto eta = x / (centralMeridianScale * A);
   const auto epsilon = y / (centralMeridianScale * A);

   const std::array<double, krugerSeriesOrder> beta = {
      1.0/2.0 * n - 2.0/3.0 * n2 + 37.0/96.0 * n3 - 1.0/360.0 * n4 - 81.0/512.0 * n5 + 96199.0/604800.0 * n6,
      1.0/48.0 * n2 + 1.0/15.0 * n3 - 437.0/1440.0 * n4 + 46.0/105.0 * n5 - 1118711.0/3870720.0 * n6,
      17.0/480.0 * n3 - 37.0/840.0 * n4 - 209.0/4480.0 * n5 + 5569.0/90720.0 * n6,
      4397.0/161280.0 * n4 - 11.0/504.0 * n5 - 830251.0/7257600.0 * n6,
      4583.0/161280.0 * n5 - 108847.0/3991680.0 * n6,
      20648693.0/638668800.0 * n6
   };

   auto epsilond = epsilon;
   auto etad = eta;
   for (int j = 1; j <= krugerSeriesOrder; j++)
   {
      // The coefficient tables are stored zero-based, but the Kruger series is one-based.
      epsilond -= beta[j - 1] * std::sin(2.0 * j * epsilon) * std::cosh(2.0 * j * eta);
      etad -= beta[j - 1] * std::cos(2.0 * j * epsilon) * std::sinh(2.0 * j * eta);
   }

   const auto sinhetad = std::sinh(etad);
   const auto sinepsilond = std::sin(epsilond), cosepsilond = std::cos(epsilond);
   const auto taud = sinepsilond / std::sqrt(sinhetad * sinhetad + cosepsilond * cosepsilond);

   double deltataui = 0.0;
   auto taui = taud;
   int iteration = 0;
   do
   {
      const auto etai = std::sinh(e * std::atanh(e * taui / std::sqrt(1.0 + taui * taui)));
      const auto tauid = taui * std::sqrt(1.0 + etai * etai) - etai * std::sqrt(1.0 + taui * taui);
      deltataui = (taud - tauid) / std::sqrt(1.0 + tauid * tauid)
         * (1.0 + (1.0 - e * e) * taui * taui) / ((1.0 - e * e) * std::sqrt(1.0 + taui * taui));
      taui += deltataui;
      iteration++;
   } while (std::abs(deltataui) > tauConvergenceTolerance && iteration < maxTauIterations);

   if (std::abs(deltataui) > tauConvergenceTolerance)
      throw std::runtime_error("UTM inverse projection failed to converge");

   const auto tau = taui;
   const auto phi = std::atan(tau);
   auto lambda = std::atan2(sinhetad, cosepsilond);

   // Karney 2011 equations 26-28 derive grid convergence and scale for the recovered point.
   auto p = 1.0;
   auto q = 0.0;
   for (int j = 1; j <= krugerSeriesOrder; j++)
   {
      p -= 2.0 * j * beta[j - 1] * std::cos(2.0 * j * epsilon) * std::cosh(2.0 * j * eta);
      q += 2.0 * j * beta[j - 1] * std::sin(2.0 * j * epsilon) * std::sinh(2.0 * j * eta);
   }

   const auto gammad = std::atan(std::tan(epsilond) * std::tanh(etad));
   const auto gammasd = std::atan2(q, p);
   const auto gamma = gammad + gammasd;

   const auto sinphi = std::sin(phi);
   const auto ksid = std::sqrt(1.0 - e * e * sinphi * sinphi) * std::sqrt(1.0 + tau * tau)
      * std::sqrt(sinhetad * sinhetad + cosepsilond * cosepsilond);
   const auto ksidd = A / a / std::sqrt(p * p + q * q);
   const auto k = centralMeridianScale * ksid * ksidd;

   const auto lambda0 = toRadians((z - 1) * 6.0 - 180.0 + 3.0);
   lambda += lambda0;

   const auto lat = std::stod(geodesy::toFixed(geodesy::toDegrees(phi), 14));
   const auto lon = std::stod(geodesy::toFixed(geodesy::wrap180(geodesy::toDegrees(lambda)), 14));
   const auto convergence = std::stod(geodesy::toFixed(geodesy::toDegrees(gamma), 9));
   const auto scale = std::stod(geodesy::toFixed(k, 12));

   LatLonUtm latLong = LatLonUtm(lat, lon, 0.0, m_datum);
   latLong.setConvergence(convergence);
   latLong.setScale(scale);

   return latLong;
}

std::string Utm::toString(int dp) const
{
   if (dp < 0)
      throw std::invalid_argument("invalid UTM precision: expected non-negative decimal places");

   std::stringstream ss;
   ss << std::setw(2) << std::setfill('0') << m_zone;

   return ss.str() + " " + hemisphereToString(m_hemisphere) + " "
      + geodesy::toFixed(m_easting, dp) + " " + geodesy::toFixed(m_northing, dp);
}
