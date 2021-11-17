#include "dms.h"

#include <iomanip>
#include <sstream>
#include <cmath>

#include "strutil.h"
#include "vector3d.hpp"

using namespace geodesy;

wchar_t Dms::_separator = L'\u202f';

wchar_t Dms::get_separator()
{
   return _separator;
}

void Dms::set_separator(wchar_t sep)
{
   _separator = sep;
}

double Dms::parse(const std::wstring& dms)
{
   // check for signed decimal degrees without NSEW, if so return it directly
   try
   {
      std::string::size_type sz;     // alias of size_t
      const double e = std::stod(dms, &sz);
      if (!isnan(e) && sz == dms.length()) return e;
   }
   catch (...)
   {
      return NAN;
   }
   
   // strip off any sign or compass dir'n & split out separate d/m/s
   const std::wstring trim_dms = strutil::strip(dms);
   std::vector<std::wstring> dms_parts = strutil::split_regex(
      std::regex_replace(trim_dms,std::wregex(L"(^-)|([NSEW]$)", std::regex_constants::icase), L""), 
      L"/[^0-9.,]+/");

   if (dms_parts[dms_parts.size() - 1].empty()) 
   {
      dms_parts.erase(dms_parts.begin() + (dms_parts.size() - 1));
   }

   if (dms_parts.empty())
       return NAN;

   std::vector<double> dms_parts_int{};
   for (auto elem: dms_parts)
   {
      try
      {
         const double d = std::stod(dms);
         dms_parts_int.emplace_back(d);
      }
      catch (...)
      {
         return NAN;
      }
   }

   // and convert to decimal degrees...
   double deg = NAN;
   switch (dms_parts.size())
   {
   case 3:
      deg = dms_parts_int[0]/1.000 + dms_parts_int[1]/60.000 + dms_parts_int[2]/3600.000;
      break;
   case 2:
      deg = dms_parts_int[0] / 1.000 + dms_parts_int[1] / 60.000;
      break;
   case 1:
      deg = dms_parts_int[0] / 1.000;
      break;
   default:
      return NAN;
   }

   if (std::regex_match(trim_dms.c_str(), std::wregex(L"(^-)|([NSEW]$"))) 
   {
      deg = -deg; // take '-', west and south as -ve
   }

   return deg;
}

std::wstring Dms::toDms(double deg, eFormat format)
{
   // give up here if we can't make a number from deg
   if (isnan(deg)) 
      return L"";

   // (unsigned result ready for appending compass dir'n)
   deg = std::abs(deg);  

   double dm = 0.0;
   double dd = 0.0;
   std::wstring dms;
   std::wstringstream stream;

   switch (format)
   {
   case D:
      {
         // left-pad with leading zeros (note may include decimals)
         stream << std::setprecision(4) << std::setfill(L'0') << deg << L"°";// round/right-pad degrees
         dms = stream.str();
      }
      break;
   case DM:
      {
         dd = std::floor(deg); // get component deg
         dm = std::fmod((deg * 60), 60); // get component min & round/right-pad

         if (std::fabs(60 - dm) < std::numeric_limits<double>::epsilon()) // check for rounding up
         {
            dm = 0.0;
            ++dd;
         }

         stream << std::setprecision(3) << std::setfill(L'0') << dd << L"°" << _separator
                << std::setprecision(2) << std::setfill(L'0') << dm << L"′";

         dms = stream.str();
      }
      break;
   case DMS:
      {
         dd = std::floor(deg); // get component deg
         dm = std::fmod((deg * 3600) / 60, 60); // get component min
         double ds = std::fmod(deg * 3600, 60); // get component sec & round/right-pad

         if(std::fabs(60 - ds) < std::numeric_limits<double>::epsilon()) // check for rounding up
         {
            ds = 0;
            ++dm;
         }

         // check for rounding up
         if (std::fabs(60 - dm) < std::numeric_limits<double>::epsilon()) // check for rounding up
         {
            dm = 0.0;
            ++dd;
         }

         stream << std::setprecision(3) << std::setfill(L'0') << dd << L"°" << _separator
                << std::setprecision(2) << std::setfill(L'0') << dm << L"′" << _separator
                << std::setprecision(2) << std::setfill(L'0') << ds << L"″";

         dms = stream.str();
      }
      break;
   case N:
      {
         stream << deg;
         dms = stream.str();
      }
      break;
   }
   return dms;
}

std::wstring Dms::toLatitude(double deg, eFormat format)
{
   const std::wstring lat = toDms(wrap90(deg), format);
   return lat.empty() ? L"-" : lat.substr(1) + _separator + (deg < 0 ? L"S" : L"N"); // knock off initial '0' for lat!
}

std::wstring Dms::toLongitude(double deg, eFormat format)
{
   const std::wstring lon = toDms(wrap180(deg), format);
   return lon.empty() ? L"-" : lon + _separator + (deg < 0 ? L"W" : L"E");
}

std::wstring Dms::toBearing(double deg, eFormat format)
{
   const std::wstring brng = toDms(wrap360(deg), format);
   // just in case rounding took us up to 360∼!
   return brng.empty() ? L"-" : std::regex_replace(brng, std::wregex(L"360"), L"0");
}

std::string Dms::compassPoint(double bearing, int precision)
{
   if (precision < 1 || precision > 3)
      throw std::range_error("invalid precision" + std::to_string(precision));

   // note precision could be extended to 4 for quarter-winds (eg NbNW), but I think they are little used
   bearing = wrap360(bearing); // normalize to range 0..360∼
   const std::string cardinals[] = {
         "N", "NNE", "NE", "ENE",
         "E", "ESE", "SE", "SSE",
         "S", "SSW", "SW", "WSW",
         "W", "WNW", "NW", "NNW" };

   // no of compass points at req＊d precision (1=>4, 2=>8, 3=>16)
   const int n = 4 * static_cast<int>(std::exp2(precision - 1)); 
   return cardinals[static_cast<int>(std::round(bearing * n / 360)) % n * 16 / n];
}

double Dms::wrap360(double degrees)
{
   if (0 <= degrees && degrees < 360)
      return degrees; // avoid rounding due to arithmetic ops if within range
    // bearing wrapping requires a sawtooth wave function with a vertical offset equal to the
    // amplitude and a corresponding phase shift; this changes the general sawtooth wave function from
    //     f(x) = (2ax/p - p/2) % p - a
    // to
    //     f(x) = (2ax/p) % p
    // where a = amplitude, p = period, % = modulo; however, c++ fmod is a remainder operator
    // not a modulo operator - for modulo, replace 'x%n' with '((x%n)+n)%n'
    const double x = degrees, a = 180, p = 360;
    return  std::fmod(std::fmod(2 * a* x / p, p) + p, p);
}

double Dms::wrap180(double degrees)
{
   if (-180 < degrees && degrees <= 180)
      return degrees; // avoid rounding due to arithmetic ops if within range
    // longitude wrapping requires a sawtooth wave function; a general sawtooth wave is
    //     f(x) = (2ax/p - p/2) % p - a
    // where a = amplitude, p = period, % = modulo; however, c++ fmod is a remainder operator
    // not a modulo operator - for modulo, replace 'x%n' with '((x%n)+n)%n'
    const double x = degrees, a = 180, p = 360;
    return std::fmod((std::fmod((2*a*x/p - p/2), p)+p) ,p) - a;
}

double Dms::wrap90(double degrees)
{
   if (-90 <= degrees && degrees <= 90)
      return degrees; // avoid rounding due to arithmetic ops if within range
    // latitude wrapping requires a triangle wave function; a general triangle wave is
    //     f(x) = 4a/p ⋅ | (x-p/4)%p - p/2 | - a
    // where a = amplitude, p = period, % = modulo; however, c++ fmod is a remainder operator
    // not a modulo operator - for modulo, replace 'x%n' with '((x%n)+n)%n'
    const double x = degrees, a = 90, p = 360;
    return 4*a/p * std::abs(std::fmod((std::fmod((x-p/4), p)+p) ,p) - p/2) - a;
}