#include "latlon_osgridref.h"
#include "algorithm.h"

#include "osgridref.h"

using namespace geodesy;

LatLonOsGridRef::LatLonOsGridRef(double lat, double lon, double height, std::optional<Datum> datum)
      : LatLonEllipsoidalDatum(lat, lon, height, datum)
{
}

OsGridRef LatLonOsGridRef::toOsGrid()
{
   // if necessary convert to OSGB36 first
   const auto point = this->datum() == datums().OSGB36
                         ? *this
                         : convertDatum(datums().OSGB36);

   const auto phi = toRadians(point.lat());
   const auto lambda = toRadians(point.lon());

   const auto [ a, b, f ] = nationalGrid.ellipsoid;            // a = 6377563.396, b = 6356256.909
   const auto phi0 = toRadians(nationalGrid.trueOrigin.lat); // latitude of true origin, 49°N
   const auto lambda0 = toRadians(nationalGrid.trueOrigin.lon); // longitude of true origin, 2°W
   const auto E0 = -nationalGrid.falseOrigin.easting;       // easting of true origin, 400km
   const auto N0 = -nationalGrid.falseOrigin.northing;      // northing of true origin, -100km
   const auto F0 = nationalGrid.scaleFactor;                // 0.9996012717

   const auto e2 = 1 - (b * b) / (a * a);                          // eccentricity squared
   const auto n = (a - b) / (a + b), n2 = n * n, n3 = n * n * n;         // n, n², n³

   const auto cosphi = std::cos(phi), sinphi = std::sin(phi);
   const auto nu = a * F0 / std::sqrt(1 - e2 * sinphi * sinphi);            // nu = transverse radius of curvature
   const auto rho = a * F0 * (1 - e2) / std::pow(1 - e2 * sinphi * sinphi, 1.5); // rho = meridional radius of curvature
   const auto eta2 = nu / rho - 1;                                    // eta = ?

   const auto Ma = (1 + n + (5 / 4) * n2 + (5 / 4) * n3) * (phi - phi0);
   const auto Mb = (3 * n + 3 * n * n + (21 / 8) * n3) * std::sin(phi - phi0) * std::cos(phi + phi0);
   const auto Mc = ((15 / 8) * n2 + (15 / 8) * n3) * std::sin(2 * (phi - phi0)) * std::cos(2 * (phi + phi0));
   const auto Md = (35 / 24) * n3 * std::sin(3 * (phi - phi0)) * std::cos(3 * (phi + phi0));
   const auto M = b * F0 * (Ma - Mb + Mc - Md);              // meridional arc

   const auto cos3phi = cosphi * cosphi * cosphi;
   const auto cos5phi = cos3phi * cosphi * cosphi;
   const auto tan2phi = std::tan(phi) * std::tan(phi);
   const auto tan4phi = tan2phi * tan2phi;

   const auto I = M + N0;
   const auto II = (nu / 2) * sinphi * cosphi;
   const auto III = (nu / 24) * sinphi * cos3phi * (5 - tan2phi + 9 * eta2);
   const auto IIIA = (nu / 720) * sinphi * cos5phi * (61 - 58 * tan2phi + tan4phi);
   const auto IV = nu * cosphi;
   const auto V = (nu / 6) * cos3phi * (nu / rho - tan2phi);
   const auto VI = (nu / 120) * cos5phi * (5 - 18 * tan2phi + tan4phi + 14 * eta2 - 58 * tan2phi * eta2);

   const auto DELTAlambda = lambda - lambda0;
   const auto DELTAlambda2 = DELTAlambda * DELTAlambda;
   const auto DELTAlambda3 = DELTAlambda2 * DELTAlambda;
   const auto DELTAlambda4 = DELTAlambda3 * DELTAlambda;
   const auto DELTAlambda5 = DELTAlambda4 * DELTAlambda;
   const auto DELTAlambda6 = DELTAlambda5 * DELTAlambda;

   auto N = I + II * DELTAlambda2 + III * DELTAlambda4 + IIIA * DELTAlambda6;
   auto E = E0 + IV * DELTAlambda + V * DELTAlambda3 + VI * DELTAlambda5;

   N = std::stod(toFixed(N, 3)); // round to mm precision
   E = std::stod(toFixed(E,3));

   try 
   {
      return OsGridRef(E, N); // note: gets truncated to SW corner of 1m grid square
   }
   catch (const std::exception& e) 
   {
      throw std::runtime_error(e.what());
   }
}

LatLonOsGridRef LatLonOsGridRef::convertDatum(const Datum& toDatum) const
{
   auto osgbED = LatLonEllipsoidalDatum::convertDatum(toDatum); // returns LatLonEllipsoidal_Datum
   return LatLonOsGridRef(osgbED.lat(), osgbED.lon(), osgbED.height(), osgbED.datum());
}
