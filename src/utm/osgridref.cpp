#include "osgridref.h"

#include "latlon_osgridref.h"

using namespace geodesy;

OsGridRef::OsGridRef(double easting, double northing)
   : m_easting(easting)
   , m_northing(northing)
{
   if (std::isnan(easting) || m_easting < 0 || m_easting > 700e3) throw std::range_error("invalid easting’");
   if (std::isnan(northing) || m_northing < 0 || m_northing > 1300e3) throw std::range_error("invalid northing’");
}

LatLonOsGridRef OsGridRef::toLatLon(Datum datum) const
{
   const auto E = m_easting, N = m_northing;
   const auto [ a, b, f] = nationalGrid.ellipsoid;            // a = 6377563.396, b = 6356256.909
   const auto phi0 = toRadians(nationalGrid.trueOrigin.lat); // latitude of true origin, 49°N
   const auto lambda0 = toRadians(nationalGrid.trueOrigin.lon); // longitude of true origin, 2°W
   const auto E0 = -nationalGrid.falseOrigin.easting;       // easting of true origin, 400km
   const auto N0 = -nationalGrid.falseOrigin.northing;      // northing of true origin, -100km
   const auto F0 = nationalGrid.scaleFactor;                // 0.9996012717

   const auto e2 = 1 - (b * b) / (a * a);                         // eccentricity squared
   const auto n = (a - b) / (a + b), n2 = n * n, n3 = n * n * n;        // n, n², n³

   auto phi = phi0, M = 0.0;
   do {
      phi = (N - N0 - M) / (a * F0) + phi;

      const auto Ma = (1 + n + (5 / 4) * n2 + (5 / 4) * n3) * (phi - phi0);
      const auto Mb = (3 * n + 3 * n * n + (21 / 8) * n3) * std::sin(phi - phi0) * std::cos(phi + phi0);
      const auto Mc = ((15 / 8) * n2 + (15 / 8) * n3) * std::sin(2 * (phi - phi0)) * std::cos(2 * (phi + phi0));
      const auto Md = (35 / 24) * n3 * std::sin(3 * (phi - phi0)) * std::cos(3 * (phi + phi0));
      M = b * F0 * (Ma - Mb + Mc - Md);                // meridional arc

   } while (std::abs(N - N0 - M) >= 0.00001);  // ie until < 0.01mm

   const auto cosphi = std::cos(phi), sinphi = std::sin(phi);
   const auto nu = a * F0 / std::sqrt(1 - e2 * sinphi * sinphi);            // nu = transverse radius of curvature
   const auto rho = a * F0 * (1 - e2) / std::pow(1 - e2 * sinphi * sinphi, 1.5); // rho = meridional radius of curvature
   const auto eta2 = nu / rho - 1;                                    // eta = ?

   const auto tanphi = std::tan(phi);
   const auto tan2phi = tanphi * tanphi, tan4phi = tan2phi * tan2phi, tan6phi = tan4phi * tan2phi;
   const auto secphi = 1 / cosphi;
   const auto nu3 = nu * nu * nu, nu5 = nu3 * nu * nu, nu7 = nu5 * nu * nu;
   const auto VII = tanphi / (2 * rho * nu);
   const auto VIII = tanphi / (24 * rho * nu3) * (5 + 3 * tan2phi + eta2 - 9 * tan2phi * eta2);
   const auto IX = tanphi / (720 * rho * nu5) * (61 + 90 * tan2phi + 45 * tan4phi);
   const auto X = secphi / nu;
   const auto XI = secphi / (6 * nu3) * (nu / rho + 2 * tan2phi);
   const auto XII = secphi / (120 * nu5) * (5 + 28 * tan2phi + 24 * tan4phi);
   const auto XIIA = secphi / (5040 * nu7) * (61 + 662 * tan2phi + 1320 * tan4phi + 720 * tan6phi);

   const auto dE = (E - E0), dE2 = dE * dE, dE3 = dE2 * dE, dE4 = dE2 * dE2, dE5 = dE3 * dE2, dE6 = dE4 * dE2, dE7 = dE5 * dE2;
   phi = phi - VII * dE2 + VIII * dE4 - IX * dE6;
   const auto lambda = lambda0 + X * dE - XI * dE3 + XII * dE5 - XIIA * dE7;

   auto point = LatLonOsGridRef(toDegrees(phi), toDegrees(lambda), 0, LatLonEllipsoidal::datums().OSGB36);

   if (datum != LatLonEllipsoidal::datums().OSGB36) {
      // if point is required in datum other than OSGB36, convert it
      point = point.convertDatum(datum);
      // convertDatum() gives us a LatLon: convert to LatLon_OsGridRef which includes toOsGrid()
      point = LatLonOsGridRef(point.lat(), point.lon(), point.height(), point.datum());
   }

   return point;
}

OsGridRef OsGridRef::parse(const std::string& gridref)
{
   std::string ref = strutil::strip(gridref);

   // check for fully numeric comma-separated gridref format
   std::smatch sm1;
   if(std::regex_match(ref, sm1, std::regex("^(\\d+),\\s*(\\d+)$")))
   {
      if(sm1.size() == 2)
      {
         return { std::stod(sm1[0].str()), std::stod(sm1[1].str()) };
      }
   }

   // validate format
   std::smatch sm2;
   if (!std::regex_match(ref, sm2, std::regex("[HNST][ABCDEFGHJKLMNOPQRSTUVWXYZ]\\s*[0-9]+\\s*[0-9]+$", std::regex_constants::icase)))
   {
      throw std::runtime_error("invalid grid reference");
   }

   // get numeric values of letter references, mapping A->0, B->1, C->2, etc:
   std::transform(ref.begin(), ref.end(), ref.begin(), toupper);


   unsigned l1 = ref.at(0) - 'A'; // 500km square
   unsigned l2 = ref.at(1) - 'A'; // 100km square
   // shuffle down letters after 'I' since 'I' is not used in grid:
   if (l1 > 7) l1--;
   if (l2 > 7) l2--;

   // convert grid letters into 100km-square indexes from false origin (grid square SV):
   const auto e100km = ((l1 - 2) % 5) * 5 + (l2 % 5);
   const auto n100km = (19 - std::floor(l1 / 5) * 5) - std::floor(l2 / 5);

   // skip grid letters to get numeric (easting/northing) part of ref
   auto en = strutil::split_regex(strutil::strip(ref.substr(2)), "\\s+");
   // if e/n not whitespace separated, split half way
   if (en.size() == 1) 
      en = { en[0].substr(0, en[0].length() / 2), en[0].substr(en[0].length() / 2) };

   // validation
   if (en[0].length() != en[1].length()) 
      throw std::runtime_error("invalid grid reference’");

   // standardise to 10-digit refs (metres)
   en[0] = strutil::padRight(en[0], 5, '0');
   en[1] = strutil::padRight(en[1], 5, '0');

   const auto e = e100km + std::stod(en[0]);
   const auto n = n100km + std::stod(en[1]);

   return OsGridRef(e, n);
}

std::string OsGridRef::toString(int digits)
{
   std::vector precisions = { 0, 2, 4, 6, 8, 10, 12, 14, 16 };

   if(0 == std::count(precisions.begin(), precisions.end(), digits))
   {
      throw std::range_error("invalid precision"); // eslint-disable-line comma-spacing
   }

   auto e = m_easting,  n = m_northing;
   // use digits = 0 to return numeric format (in metres) - note northing may be >= 1e7
   if (digits == 0)
   {
      std::stringstream sePad, snPad;
      sePad << std::setw(6) << std::setfill('0') << static_cast<int>(e);
      sePad << std::setw(3) << std::setfill('0') << e - static_cast<int>(e);
      snPad << std::setw(6) << std::setfill('0') << static_cast<int>(n);
      snPad << std::setw(3) << std::setfill('0') << e - static_cast<int>(n);
      return sePad.str() + ", " + snPad.str();
   }

   // get the 100km-grid indices
   const int e100km = std::floor(e / 100000), n100km = std::floor(n / 100000);

   // translate those into numeric equivalents of the grid letters
   // translate those into numeric equivalents of the grid letters
   int l1 = (19 - n100km) - (19 - n100km) % 5 + (e100km + 10) / 5;
   int l2 = (19 - n100km) * 5 % 25 + e100km % 5;

   // compensate for skipped 'I' and build grid letter-pairs
   if (l1 > 7) ++l1;
   if (l2 > 7) ++l2;

   std::string letterPair;
   letterPair += l1 + 'A';
   letterPair += l2 + 'A';

   // strip 100km-grid indices from easting & northing, and reduce precision
   e = std::floor(std::fmod(e, 100000) / std::pow(10, 5 - digits / 2));
   n = std::floor(std::fmod(n, 100000) / std::pow(10, 5 - digits / 2));

   // pad eastings & northings with leading zeros
   std::string se = strutil::padLeft(std::to_string(e), digits / 2, '0');
   std::string sn = strutil::padLeft(std::to_string(n), digits / 2, '0');

   return letterPair + " " + se + " " + sn;
}
