/**********************************************************************************
*  MIT License                                                                    *
*                                                                                 *
*  Copyright (c) 2021 Binbin Song <ssln.jzs@gmail.com>                            *
*                                                                                 *
*  Geodesy tools for conversions between (historical) datums                      *
*  (c) Chris Veness 2005-2019                                                     *
*  www.movable-type.co.uk/scripts/latlong-utm-mgrs.html                           *
*  www.movable-type.co.uk/scripts/geodesy-library.html#utm                        *
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

#include <iomanip>
#include <regex>
#include <sstream>

#include "strutil.h"
#include "algorithm.h"
#include "latlon_utm.h"

using namespace geodesy;

Utm::Utm(int zone, Hemisphere h, double easting, double northing, std::optional<Datum> datum,
                  std::optional<double> convergence, std::optional<double> scale, bool verifyEN)
   : m_zone(zone), m_hemisphere(h), m_easting(easting), m_northing(northing)
   , m_datum(datum), m_convergence(convergence), m_scale(scale)
{
   if (!(1 <= zone && zone <= 60))
      throw std::range_error("invalid UTM zone: out of range [1, 60]");

   if (verifyEN)
   {
      // range-check E/N values
      if (!(0 <= easting && easting <= 1000e3))
      {
         throw std::range_error("invalid UTM easting: out of range [0, 1000e3]");
      }

      switch (h)
      {
      case Hemisphere::N:
         {
            if (0 <= northing && northing < 9328094)
            {
               throw std::range_error("invalid UTM northing: out of range [0, 9328094]");
            }
         }
         break;
      case Hemisphere::S:
         if (1118414 <= northing && northing < 10000e3)
         {
            throw std::range_error("invalid UTM northing: out of range [1118414, 10000e3]");
         }
         break;
      }
   }

   if (!datum.has_value())
   {
      throw std::invalid_argument("unrecognized datum: empty");
   }
}

Utm Utm::parse(const std::string& utmCoord, std::optional<Datum> datum)
{
   // match separate elements (separated by whitespace)
   const std::string utm = strutil::strip(utmCoord);
   std::smatch sm;
   std::regex_match(utm, sm, std::regex("\\S+"));

   if (sm.size() != 4) 
      throw std::invalid_argument("invalid UTM coordinate");

   const std::string zone = sm[0].str(), hemisphere = sm[1].str(), easting = sm[2].str(), northing = sm[3].str();

   try
   {
      Hemisphere h = sm[1].str() == "N" || sm[1].str() == "n" ? Hemisphere::N : Hemisphere::S;
      return Utm(std::stoi(zone), h, std::stod(easting), std::stod(northing), datum); // 'new this' as may return subclassed types
   }
   catch (const std::exception& e)
   {
      throw std::invalid_argument("invalid UTM coordinate");
   }
}

LatLonUtm Utm::toLatLon() const
{
   const auto z = m_zone;
   const auto h = m_hemisphere;

   const auto falseEasting = 500e3, falseNorthing = 10000e3;

   const auto a = m_datum->ellipsoid.a, f = m_datum->ellipsoid.f; // WGS-84: a = 6378137, f = 1/298.257223563;

   const auto k0 = 0.9996; // UTM scale on the central meridian

   const auto x = m_easting - falseEasting;                            // make x ± relative to central meridian
   const auto y = h == Hemisphere::S ? m_northing - falseNorthing : m_northing; // make y ± relative to equator

   // ---- from Karney 2011 Eq 15-22, 36:

   const auto e = std::sqrt(f * (2 - f)); // eccentricity
   const auto n = f / (2 - f);        // 3rd flattening
   const auto n2 = n * n, n3 = n * n2, n4 = n * n3, n5 = n * n4, n6 = n * n5;

   const auto A = a / (1 + n) * (1 + 1 / 4 * n2 + 1 / 64 * n4 + 1 / 256 * n6); // 2πA is the circumference of a meridian

   const auto eta = x / (k0 * A);
   const auto epsilon = y / (k0 * A);

   const double beta[6] = {
      1.0/2.0 * n - 2.0/3.0 * n2 + 37.0/96.0 * n3 - 1.0/360.0 * n4 - 81.0/512.0 * n5 + 96199.0/604800.0 * n6,
      1.0/48.0 * n2 + 1.0/15.0 * n3 - 437.0/1440.0 * n4 + 46.0/105.0 * n5 - 1118711.0/3870720.0 * n6,
      17.0/480.0 * n3 - 37.0/840.0 * n4 - 209.0/4480.0 * n5 + 5569.0/90720.0 * n6,
      4397.0/161280.0 * n4 - 11.0/504.0 * n5 - 830251.0/7257600.0 * n6,
      4583.0/161280.0 * n5 - 108847.0/3991680.0 * n6,
      20648693.0/638668800.0 * n6 };

   auto epsilond = epsilon;
   for (auto j = 0; j < 6; j++) epsilond -= beta[j] * std::sin(2 * j * epsilon) * std::cosh(2 * j * eta);

   auto etad = eta;
   for (auto j = 0; j < 6; j++) etad -= beta[j] * std::cos(2 * j * epsilon) * std::sinh(2 * j * eta);

   const auto sinhetad = std::sinh(etad);
   const auto sinepsilond = std::sin(epsilond), cosepsilond = std::cos(epsilond);

   const auto taud = sinepsilond / std::sqrt(sinhetad * sinhetad + cosepsilond * cosepsilond);

   double deltataui;
   auto taui = taud;
   do {
      const auto etai = std::sinh(e * std::atanh(e * taui / std::sqrt(1 + taui * taui)));
      const auto tauid = taui * std::sqrt(1 + etai * etai) - etai * std::sqrt(1 + taui * taui);
      deltataui = (taud - tauid) / std::sqrt(1 + tauid * tauid)
         * (1 + (1 - e * e) * taui * taui) / ((1 - e * e) * std::sqrt(1 + taui * taui));
      taui += deltataui;
   } while (std::abs(deltataui) > 1e-12); // using IEEE 754 δτi -> 0 after 2-3 iterations
   // note relatively large convergence test as δτi toggles on ±1.12e-16 for eg 31 N 400000 5000000
   const auto tau = taui;

   const auto phi = std::atan(tau);

   auto lambda = std::atan2(sinhetad, cosepsilond);

   // ---- convergence: Karney 2011 Eq 26, 27

   auto p = 1.0;
   for (auto j = 0; j < 6; j++) p -= 2 * j * beta[j] * std::cos(2 * j * epsilon) * std::cosh(2 * j * eta);
   auto q = 0.0;
   for (auto j = 0; j < 6; j++) q += 2 * j * beta[j] * std::sin(2 * j * epsilon) * std::sinh(2 * j * eta);

   const auto gammad = std::atan(std::tan(epsilond) * std::tanh(etad));
   const auto gammasd = std::atan2(q, p);

   const auto gamma = gammad + gammasd;

   // ---- scale: Karney 2011 Eq 28

   const auto sinphi = std::sin(phi);
   const auto ksid = std::sqrt(1 - e * e * sinphi * sinphi) * std::sqrt(1 + tau * tau) * std::sqrt(sinhetad * sinhetad + cosepsilond * cosepsilond);
   const auto ksidd = A / a / std::sqrt(p * p + q * q);

   const auto k  = k0 * ksid * ksidd;

   // ------------

   const auto lambda0 = toRadians((z - 1) * 6 - 180 + 3); // longitude of central meridian
   lambda += lambda0; // move λ from zonal to global coordinates

   // round to reasonable precision
   const auto lat = std::stod(geodesy::toFixed(geodesy::toDegrees(phi), 14)); // nm precision (1nm = 10^-14°)
   const auto lon = std::stod(geodesy::toFixed(geodesy::toDegrees(lambda),14)); // (strictly lat rounding should be φ⋅cosφ!)
   const auto convergence = std::stod(geodesy::toFixed(geodesy::toDegrees(gamma),9));
   const auto scale = std::stod(geodesy::toFixed(k,12));

   LatLonUtm latLong = LatLonUtm(lat, lon, 0, m_datum);
   // ... and add the convergence and scale into the LatLon object ... wonderful JavaScript!
   latLong.setConvergence(convergence);
   latLong.setScale(scale);

   return latLong;
}

std::string Utm::toString(int dp) const
{
   std::stringstream ss;
   ss << std::setw(2) << std::left << std::setfill('0') << m_zone;

   return ss.str() + " " + (m_hemisphere == Hemisphere::N ? "N " : "S ")
      + geodesy::toFixed(m_easting, dp) + geodesy::toFixed(m_northing, dp);
}
