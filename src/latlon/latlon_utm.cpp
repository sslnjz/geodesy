#include "latlon_utm.h"

#include "utm.h"
#include "algorithm.h"

using namespace geodesy;

LatLonUtm::LatLonUtm(double lat, double lon, double height, std::optional<Datum> datum,
   std::optional<ReferenceFrame> reference, std::optional<std::string> epoch)
   : LatLonEllipsoidal(lat, lon, height, datum, reference, epoch)
   , m_convergence(0.0)
   , m_scale(0.0)
{
}

Utm LatLonUtm::toUtm(std::optional<int> zoneOverride) const
{
   if (!(-80 <= lat() && lat() <= 84)) 
      throw std::out_of_range("latitude outside UTM limits");

   const auto falseEasting = 500e3;

   auto zone = zoneOverride ? *zoneOverride : static_cast<int>(std::floor((lon() + 180) / 6) + 1); // longitudinal zone
   auto lambda0 = toRadians((zone - 1) * 6 - 180 + 3); // longitude of central meridian

   const auto index = static_cast<size_t>(std::floor(lat() / 8 + 10));
   const auto mgrsLatBands = "CDEFGHJKLMNPQRSTUVWXX";
   const char latBand = index < strlen(mgrsLatBands) ? mgrsLatBands[index] : 0;

   const auto radians6 = toRadians(6);
   // adjust zone & central meridian for Norway
   if (zone == 31 && latBand == 'V' && lon() >= 3) { zone++; lambda0 += toRadians((6)); }
   // adjust zone & central meridian for Svalbard
   if (zone == 32 && latBand == 'X' && lon() < 9) { zone--; lambda0 -= radians6; }
   if (zone == 32 && latBand == 'X' && lon() >= 9) { zone++; lambda0 += radians6; }
   if (zone == 34 && latBand == 'X' && lon() < 21) { zone--; lambda0 -= radians6; }
   if (zone == 34 && latBand == 'X' && lon() >= 21) { zone++; lambda0 += radians6; }
   if (zone == 36 && latBand == 'X' && lon() < 33) { zone--; lambda0 -= radians6; }
   if (zone == 36 && latBand == 'X' && lon() >= 33) { zone++; lambda0 += radians6; }

   const auto phi = toRadians(lat());      // latitude ± from equator
   const auto lambda = toRadians(lon() - lambda0); // longitude ± from central meridian

   // allow alternative ellipsoid to be specified
   const auto ellipsoid = m_datum ? m_datum->ellipsoid : ellipsoids().WGS84;
   const auto a = ellipsoid.a, f = ellipsoid.f; // WGS-84: a = 6378137, f = 1/298.257223563;
   const auto kappa0 = 0.9996; // UTM scale on the central meridian

   // ---- easting, northing: Karney 2011 Eq 7-14, 29, 35:
   const auto e = std::sqrt(f * (2 - f)); // eccentricity
   const auto n = f / (2 - f);        // 3rd flattening
   const auto n2 = n * n, n3 = n * n2, n4 = n * n3, n5 = n * n4, n6 = n * n5;

   const auto coslambda = std::cos(lambda), sinlambda = std::sin(lambda), tanlambda = std::tan(lambda);

   const auto tau = std::tan(phi); // τ ≡ tanφ, τʹ ≡ tanφʹ; prime (ʹ) indicates angles on the conformal sphere
   const auto sigma = std::sinh(e * std::atanh(e * tau / std::sqrt(1 + tau * tau)));

   const auto taud = tau * std::sqrt(1 + sigma * sigma) - sigma * std::sqrt(1 + tau * tau);

   const auto ksid = std::atan2(taud, coslambda);
   const auto etad = std::asinh(sinlambda / std::sqrt(taud * taud + coslambda * coslambda));

   const auto A = a / (1 + n) * (1 + 1.0/4.0 * n2 + 1.0/64.0 * n4 + 1.0/256.0 * n6); // 2πA is the circumference of a meridian

   const double alpha[] = {
      1.0/2.0 * n - 2.0/3.0 * n2 + 5.0/16.0 * n3 + 41.0/180.0 * n4 - 127.0/288.0 * n5 + 7891.0 / 37800.0 * n6,
                    13.0/48.0 * n2 - 3.0/5 * n3 + 557.0/1440.0 * n4 + 281.0/630.0 * n5 - 1983433.0/1935360.0 * n6,
                                 61.0/240.0 * n3 - 103.0/140.0 * n4 + 15061.0/26880.0 * n5 + 167603.0/ 181440.0 * n6,
                                             49561.0/161280.0 * n4 - 179.0 /168.0 * n5 + 6601661.0/7257600.0 * n6,
                                                                  34729.0/ 80640.0 * n5 - 3418889.0/1995840.0 * n6,
                                                                                          212378941.0/ 319334400.0 * n6 };

   double epsilon = ksid;
   for (auto j = 0; j < 6; j++) epsilon += alpha[j] * std::sin(2 * j * ksid) * std::cosh(2 * j * etad);

   double eta = etad;
   for (auto j = 0; j < 6; j++) eta += alpha[j] * std::cos(2 * j * ksid) * std::sinh(2 * j * etad);

   auto x = kappa0 * A * eta;
   auto y = kappa0 * A * epsilon;

   // ---- convergence: Karney 2011 Eq 23, 24
   auto pd = 1.0;
   for (auto j = 0; j < 6; j++) pd += 2.0 * j * alpha[j] * std::cos(2 * j * ksid) * std::cosh(2 * j * etad);
   auto qd = 0.0;
   for (auto j = 0; j < 6; j++) qd += 2.0 * j * alpha[j] * std::sin(2 * j * ksid) * std::sinh(2 * j * etad);

   const auto gammad = std::atan(taud / std::sqrt(1 + taud * taud) * tanlambda);
   const auto gammasd = std::atan2(qd, pd);
   const auto gamma = gammad + gammasd;

   // ---- scale: Karney 2011 Eq 25
   const auto sinφ = std::sin(phi);
   const auto kappad = std::sqrt(1 - e * e * sinφ * sinφ) * std::sqrt(1 + tau * tau) / std::sqrt(taud * taud + coslambda * coslambda);
   const auto kappasd = A / a * std::sqrt(pd * pd + qd * qd);
   const auto kappa = kappa0 * kappad * kappasd;

   // shift x/y to false origins
   x = x + falseEasting;             // make x relative to false easting
   if (const auto falseNorthing = 10000e3; y < 0) y = y + falseNorthing; // make y in southern hemisphere relative to false northing

   // round to reasonable precision
   x = std::stod(toFixed(x,9)); // nm precision
   y = std::stod(toFixed(y,9)); // nm precision
   const auto convergence = std::stod(toFixed(toDegrees(gamma), 9));
   const auto scale = std::stod(toFixed(kappa,12));
   const auto h = lat() >= 0 ? Utm::Hemisphere::N : Utm::Hemisphere::S; // hemisphere

   return Utm(zone, h, x, y, m_datum, convergence, scale, zoneOverride.has_value());
}

void LatLonUtm::setConvergence(double convergence)
{
   m_convergence = convergence;
}

void LatLonUtm::setScale(double scale)
{
   m_scale = scale;
}
