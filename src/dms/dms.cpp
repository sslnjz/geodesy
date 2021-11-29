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
#include "dms.h"

#include <iomanip>
#include <sstream>
#include <cmath>
#include <locale>

#include "strutil.h"
#include "vector3d.h"

using namespace geodesy;

std::string& Dms::_separator = *new std::string("\u202f");

std::string Dms::get_separator()
{
   return _separator;
}

void Dms::set_separator(const std::string& sep)
{
   _separator = sep;
}

double Dms::parse(const std::string& dms)
{
   // check for signed decimal degrees without NSEW, if so return it directly
   if (std::regex_search(strutil::strip(dms), std::regex("^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$")))
   {
      return std::stod(dms);
   }

   // strip off any sign or compass dir'n & split out separate d/m/s
   const std::string trim_dms = strutil::strip(dms);
   std::vector<std::string> dms_parts = strutil::split_regex(
      std::regex_replace(trim_dms,std::regex("(^-)|([NSEW]$)",
                             std::regex_constants::icase), ""),
      "[^0-9.,]+");

   if (dms_parts[dms_parts.size() - 1].empty()) 
   {
      dms_parts.erase(dms_parts.begin() + (dms_parts.size() - 1));
   }

   if (dms_parts.empty())
       return NAN;

   std::vector<double> dms_parts_d{};
   for (auto elem: dms_parts)
   {
      try
      {
         const double d = std::stod(elem);
         dms_parts_d.emplace_back(d);
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
      deg = dms_parts_d[0]/1.000 + dms_parts_d[1]/60.000 + dms_parts_d[2]/3600.000;
      break;
   case 2:
      deg = dms_parts_d[0] / 1.000 + dms_parts_d[1] / 60.000;
      break;
   case 1:
      deg = dms_parts_d[0] / 1.000;
      break;
   default:
      return NAN;
   }

   if (std::regex_search(dms, std::regex("(^-)|([SW]$)")))
   {
      deg = -deg; // take '-', west and south as -ve
   }

   return deg;
}

std::string Dms::toDms(double deg, eFormat format, std::optional<int> dp)
{
   // give up here if we can't make a number from deg
   if (isnan(deg)) 
      return "";

   // default values
   if (dp == std::nullopt) 
   {
      switch (format)
      {
      case D: 
         dp = 4;
         break;
      case DM: 
         dp = 2;
         break;
      case DMS: 
         dp = 0;
         break;
      case N:
          dp = 4;
          break;
      }
   }

   // (unsigned result ready for appending compass dir'n)
   deg = std::fabs(deg);  

   std::string dms;
   switch (format)
   {
   case D:
      {
         std::string d;
         std::stringstream dss;
         // left-pad with leading zeros (note may include decimals)
         dss << std::fixed << std::setprecision(*dp) << std::setfill('0') << deg;// round/right-pad degrees
         d = dss.str();
         if (deg < 100) d = "0" + d;                    // left-pad with leading zeros (note may include decimals)
         if (deg < 10) d = "0" + d;
         dms = d + "°";
      }
      break;
   case DM:
      {
         std::string  m;
         std::stringstream dss, mss;

         int dd = static_cast<int>(std::floor(deg)); // get component deg
         double dm = std::fmod((deg * 60), 60); // get component min & round/right-pad

         if (std::fabs(60.0 - dm) < (*dp + 1) * 0.1) // check for rounding up
         {
            dm = 0.0;
            ++dd;
         }
         
         dss << std::setfill('0') << std::setw(3) << dd;
         mss << std::fixed << std::setprecision(*dp) << dm;
         m = mss.str();

         if (dm < 10)
            m = "0" + m;
         dms = dss.str() + "°" + _separator + m + "′";
      }
      break;
   case DMS:
      {
         std::string d, m, s;
         std::stringstream dss, mss, sss;

         int dd = static_cast<int>(std::floor(deg)); // get component deg
         int dm = static_cast<int>(std::fmod((deg * 3600) / 60, 60)); // get component min
         double ds = std::fmod(deg * 3600, 60); // get component sec & round/right-pad

         if(std::fabs(60.0 - ds) <= (*dp + 1) * 0.1 ) // check for rounding up
         {
            ds = 0;
            ++dm;
         }

         // check for rounding up
         if (dm == 60) // check for rounding up
         {
            dm = 0.0;
            ++dd;
         }

         dss << std::setfill('0') << std::setw(3) << dd;
         mss << std::setfill('0') << std::setw(2) << dm;
         sss << std::fixed << std::setprecision(*dp) << std::setfill('0') << ds;
         s = sss.str();

         if (ds < 10) s = "0" + s;

         dms = dss.str() + "°" + mss.str() + "′" + s + "″";
      }
      break;
   case N:
      {
         std::stringstream dmss;
         dmss << std::fixed << std::setprecision(*dp) << std::setfill('0') << deg;
         dms = dmss.str();
      }
      break;
   }
   return dms;
}

std::string Dms::toLat(double deg, eFormat format, std::optional<int> dp)
{
   const std::string lat = toDms(wrap90(deg), format, dp);
   return lat.empty() ? "-" : lat.substr(1) + _separator + (deg < 0 ? "S" : "N"); // knock off initial '0' for lat!
}

std::string Dms::toLon(double deg, eFormat format, std::optional<int> dp)
{
   const std::string lon = toDms(wrap180(deg), format, dp);
   return lon.empty() ? "-" : lon + _separator + (deg < 0 ? "W" : "E");
}

std::string Dms::toBearing(double deg, eFormat format, std::optional<int> dp)
{
   const std::string brng = toDms(wrap360(deg), format, dp);
   // just in case rounding took us up to 360∼!
   return brng.empty() ? "-" : std::regex_replace(brng, std::regex("360"), "0");
}

std::string Dms::fromLocale(const std::string& str)
{
   std::locale::global(std::locale("en_US.UTF-8"));
   std::stringstream ss;
   ss << str;
   return ss.str();
}

std::string Dms::toLocale(const std::string& str)
{
   std::locale::global(std::locale(""));
   std::stringstream ss;
   ss << str;
   return ss.str();
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