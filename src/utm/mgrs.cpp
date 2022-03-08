#include "mgrs.h"

#include "utm_mgrs.h"
#include "latlon_utm.h"
#include "strutil.h"

using namespace geodesy;

Mgrs::Mgrs(int zone, char band, char e100k, char n100k, double easting,
           double northing, Datum datum)
{
   if (!(1 <= zone && zone <= 60)) throw std::range_error("invalid MGRS zone");
   std::vector<std::string> errors; // check & report all other possible errors rather than reporting one-by-one

   if (latBands.find(band) == std::string::npos) 
      errors.emplace_back("invalid MGRS band");
   if (e100kLetters[(zone - 1) % 3].find(e100k) == std::string::npos) 
      errors.emplace_back("invalid MGRS 100km grid square column e100k for zone");
   if (n100kLetters[0].find(n100k) == std::string::npos) 
      errors.emplace_back("invalid MGRS 100km grid square row");
   if (std::isnan(easting)) 
      errors.emplace_back("invalid MGRS easting");
   if (std::isnan(northing)) 
      errors.emplace_back("invalid MGRS northing");
   if (!datum || !datum.ellipsoid)
      errors.emplace_back("unrecognised datum");

   if (!errors.empty())
   {
      auto collect_errors = [&]
      {
         std::string e;
         for (const auto& item : errors) e += item;
         return e;
      };
      throw std::runtime_error(collect_errors());
   }

   m_zone = zone;
   m_band = band;
   m_e100k = e100k;
   m_n100k = n100k;
   m_easting = easting;
   m_northing = northing;
   m_datum = datum;
}

UtmMgrs Mgrs::toUtm()
{
   const Utm::Hemisphere hemisphere = m_band >= 'N' ? Utm::Hemisphere::N : Utm::Hemisphere::S;

   // get easting specified by e100k (note +1 because eastings start at 166e3 due to 500km false origin)
   const auto col = e100kLetters[(m_zone - 1) % 3].find(m_e100k) + 1;
   const auto e100kNum = static_cast<double>(col) * 100e3; // e100k in metres

   // get northing specified by n100k
   const auto row = n100kLetters[(m_zone - 1) % 2].find(m_n100k);
   const auto n100kNum = static_cast<double>(row) * 100e3; // n100k in metres

   // get latitude of (bottom of) band
   const double latBand = (latBands.find(m_band) - 10) * 8;

   // get northing of bottom of band, extended to include entirety of bottom-most 100km square
   const auto nBand = std::floor(LatLonUtm(latBand, 3.0).toUtm().northing() / 100e3) * 100e3;

   // 100km grid square row letters repeat every 2,000km north; add enough 2,000km blocks to
   // get into required band
   auto n2M = 0.0; // northing of 2,000km block
   while (n2M + n100kNum + m_northing < nBand) n2M += 2000e3;

   return UtmMgrs(m_zone, hemisphere, e100kNum + m_easting, n2M + n100kNum + m_northing, m_datum);
}

std::string Mgrs::toString(unsigned int digits)
{
   if ( digits > 10 || (digits % 2) + 1 == 2 )
      throw std::range_error("invalid precision");

   // truncate to required precision
   const auto eRounded = std::floor(m_easting/std::pow(10, 5-digits/2));
   const auto nRounded = std::floor(m_northing/std::pow(10, 5-digits/2));

   // ensure leading zeros
   const auto zPadded = strutil::padLeft(std::to_string(m_zone), 2, '0');
   const auto ePadded = strutil::padLeft(std::to_string(eRounded), digits/2, '0');
   const auto nPadded = strutil::padLeft(std::to_string(nRounded), digits/2, '0');

   return zPadded + m_band + " " + m_e100k + m_n100k + " " + ePadded + " " + nPadded;
}

Mgrs Mgrs::parse(const std::string &mgrsGridRef)
{
   if (mgrsGridRef.empty()) throw std::invalid_argument("invalid MGRS grid reference");

   std::string mgrs = "";
   // check for military-style grid reference with no separators
   if (!strutil::strip(mgrsGridRef).empty()) {
      try
      {
         double num = std::stod(mgrsGridRef.substr(0, 2));
      }
      catch (const std::invalid_argument& e)
      {
         throw std::invalid_argument("invalid MGRS grid reference");
      }

      auto en = strutil::strip(mgrsGridRef).substr(5); // get easting/northing following zone/band/100ksq
      en = en.substr(0, en.length()/2)+' '+en.substr(-en.length()/2); // separate easting/northing
      mgrs = mgrsGridRef.substr(0, 3)+' '+mgrsGridRef.substr(3, 5)+' '+en; // insert spaces
   }

   // match separate elements (separated by whitespace)
   const auto ref = strutil::split(mgrs, ' ');
   if (ref.empty() || ref.size() != 4)
      throw std::invalid_argument("invalid MGRS grid reference");

   // split gzd into zone/band
   const auto gzd = ref[0];
   const auto zone = gzd.substr(0, 2);
   const auto band = gzd.substr(2, 3);

   // split 100km letter-pair into e/n
   const auto en100k = ref[1];
   const auto e100k = en100k.substr(0, 1);
   const auto n100k = en100k.substr(1, 2);

   auto e = ref[2], n = ref[3];

   // standardise to 10-digit refs - ie metres) (but only if < 10-digit refs, to allow decimals)
   e = e.length()>=5 ?  e : (e+"00000").substr(0, 5);
   n = n.length()>=5 ?  n : (n+"00000").substr(0, 5);

   return Mgrs(std::stod(zone), std::stod(band), std::stod(e100k), std::stod(n100k), std::stod(e), std::stod(n));
}
