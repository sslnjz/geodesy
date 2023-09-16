#include "utm_mgrs.h"

#include "mgrs.h"
#include "latlon_utm.h"

using namespace geodesy;

UtmMgrs::UtmMgrs(int zone, Hemisphere h, double easting, double northing, std::optional<Datum> datum)
   : Utm(zone, h, easting, northing, datum)
{
}

Mgrs UtmMgrs::toMgrs() const
{
   // MGRS zone is same as UTM zone
   const int zone = m_zone;

   // convert UTM to lat/long to get latitude to determine band
   const auto latlong = toLatLon();
   // grid zones are 8бу tall, 0буN is 10th band
   const auto band = latBands.at(static_cast<int>(std::floor(latlong.lat() / 8 + 10))); // latitude band

   // columns in zone 1 are A-H, zone 2 J-R, zone 3 S-Z, then repeating every 3rd zone
   const int col = static_cast<int>(std::floor(m_easting / 100e3));
   // (note -1 because eastings start at 166e3 due to 500km false origin)
   const auto e100k = e100kLetters[(zone - 1) % 3].at(col - 1);

   // rows in even zones are A-V, in odd zones are F-E
   const int row = static_cast<int>(std::fmod(std::floor(m_northing / 100e3), 20));
   const auto n100k = n100kLetters[(zone - 1) % 2].at(row);

   // truncate easting/northing to within 100km grid square
   auto easting = std::fmod(m_easting , 100e3);
   auto northing = std::fmod(m_northing, 100e3);

   // round to nm precision
   easting = std::stod(geodesy::toFixed(easting, 6));
   northing = std::stod(geodesy::toFixed(northing, 6));

   return Mgrs(zone, band, e100k, n100k, easting, northing);
}
