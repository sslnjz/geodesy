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
#ifndef DMS_H
#define DMS_H

#include <string>
#include <regex>
#include <optional>
#include <cmath>

#include "strutil.h"

namespace geodesy
{
   // @note: Unicode Degree = U+00B0. Prime = U+2032, Double prime = U+2033
   class Dms
   {
   public:
      Dms() = default;
      virtual ~Dms() = default;

      enum eFormat
      {
         N   = 0b1000,      // signed numeric
         D   = 0b0100,      // deg
         DM  = 0b0010,      // deg+min
         DMS = 0b0000       // deg+min+sec
      };

      /**
       * Separator character to be used to separate degrees, minutes, seconds, and cardinal directions.
       *
       * Default separator is U+202F ‘narrow no-break space’.
       *
       * To change this (e.g. to empty string or full space), set Dms.separator prior to invoking
       * formatting.
       *
       * @example
       *   import LatLon, { Dms } from '/js/geodesy/latlon-spherical.js';
       *   const p = new LatLon(51.2, 0.33).toString('dms');  // 51° 12′ 00″ N, 000° 19′ 48″ E
       *   Dms.separator = '';                                // no separator
       *   const pʹ = new LatLon(51.2, 0.33).toString('dms'); // 51°12′00″N, 000°19′48″E
       */
      static std::string separator();
      static void setSeparator(const std::string& sep);

      /**
       * Parses string representing degrees/minutes/seconds into numeric degrees.
       *
       * This is very flexible on formats, allowing signed decimal degrees, or deg-min-sec optionally
       * suffixed by compass direction (NSEW); a variety of separators are accepted. Examples -3.62,
       * '3 37 12W', '3°37′12″W'.
       *
       * Thousands/decimal separators must be comma/dot; use Dms.fromLocale to convert locale-specific
       * thousands/decimal separators.
       *
       * @param   {string|number} dms  - Degrees or deg/min/sec in variety of formats.
       * @returns {number}      - Degrees as decimal number.
       *
       * @example
       *   const double lat = Dms.parse('51° 28′ 40.37″ N');
       *   const double lon = Dms.parse('000° 00′ 05.29″ W');
       */

      template<class T>
      using ENABLE = std::enable_if_t<std::is_arithmetic_v<T> || std::is_convertible_v<T, std::string>>;
      template<class T, typename = ENABLE<T>>
      static double parse(const T& dms);

      /**
       * Converts decimal degrees to deg/min/sec format
       *  - degree, prime, double-prime symbols are added, but sign is discarded, though no compass
       *    direction is added.
       *  - degrees are zero-padded to 3 digits; for degrees latitude, use .slice(1) to remove leading
       *    zero.
       *
       * @private
       * @param   {number} deg - Degrees to be formatted as specified.
       * @param   {string} [format=dms] - Return value as 'd', 'dm', 'dms' for deg, deg+min, deg+min+sec.
       * @param   {number} [dp=4|2|0] - Number of decimal places to use – default 4 for d, 2 for dm, 0 for dms.
       * @returns {string} Degrees formatted as deg/min/secs according to specified format.
       */
      static std::string toDms(double deg, eFormat format = D, std::optional<int> dp = std::nullopt);
      
      /**
       * Converts numeric degrees to deg/min/sec latitude (2-digit degrees, suffixed with N/S).
       *
       * @param   {number} deg - Degrees to be formatted as specified.
       * @param   {string} [format=d] - Return value as 'd', 'dm', 'dms' for deg, deg+min, deg+min+sec.
       * @returns {string} Degrees formatted as deg/min/secs according to specified format.
       *
       * @example
       *   const std::string lat = Dms::toLat(-3.62, 'dms'); // 3°37′12″S
       */
      static std::string toLat(double deg, eFormat format = D, std::optional<int> dp = std::nullopt);

      /**
       * Convert numeric degrees to deg/min/sec longitude (3-digit degrees, suffixed with E/W).
       *
       * @param   {number} deg - Degrees to be formatted as specified.
       * @param   {string} [format=d] - Return value as 'd', 'dm', 'dms' for deg, deg+min, deg+min+sec.
       * @returns {string} Degrees formatted as deg/min/secs according to specified format.
       *
       * @example
       *   const std::string lon = Dms::toLon(-3.62, 'dms'); // 3°37′12″W
       */
      static std::string toLon(double deg, eFormat format = D, std::optional<int> dp = std::nullopt);

      /**
       * Converts numeric degrees to deg/min/sec as a bearing (0°..360°).
       *
       * @param   {number} deg - Degrees to be formatted as specified.
       * @param   {string} [format=d] - Return value as 'd', 'dm', 'dms' for deg, deg+min, deg+min+sec.
       * @returns {string} Degrees formatted as deg/min/secs according to specified format.
       *
       * @example
       *   const std::string lon = Dms::toBearing(-3.62, 'dms'); // 356°22′48″
       */
      [[nodiscard]] static std::string toBearing(double deg, eFormat format = D, std::optional<int> dp = std::nullopt);

      /**
       * Converts DMS string from locale thousands/decimal separators to JavaScript comma/dot separators
       * for subsequent parsing.
       *
       * Both thousands and decimal separators must be followed by a numeric character, to facilitate
       * parsing of single lat/long string (in which whitespace must be left after the comma separator).
       *
       * @param   {string} str - Degrees/minutes/seconds formatted with locale separators.
       * @returns {string} Degrees/minutes/seconds formatted with standard Javascript separators.
       *
       * @example
       *   const lat = Dms.fromLocale('51°28′40,12″N');                          // '51°28′40.12″N' in France
       *   const p = new LatLon(Dms.fromLocale('51°28′40,37″N, 000°00′05,29″W'); // '51.4779°N, 000.0015°W' in France
       */
      static std::string fromLocale(const std::string& str);

      /**
       * Converts DMS string from JavaScript comma/dot thousands/decimal separators to locale separators.
       *
       * Can also be used to format standard numbers such as distances.
       *
       * @param   {string} str - Degrees/minutes/seconds formatted with standard Javascript separators.
       * @returns {string} Degrees/minutes/seconds formatted with locale separators.
       *
       * @example
       *   const Dms.toLocale('123,456.789');                   // '123.456,789' in France
       *   const Dms.toLocale('51°28′40.12″N, 000°00′05.31″W'); // '51°28′40,12″N, 000°00′05,31″W' in France
       */
      static std::string toLocale(const std::string& str);

      /**
       * Returns compass point (to given precision) for supplied bearing.
       *
       * @param   {number} bearing - Bearing in degrees from north.
       * @param   {number} [precision=3] - Precision (1:cardinal / 2:intercardinal / 3:secondary-intercardinal).
       * @returns {string} Compass point for supplied bearing.
       *
       * @example
       *   const std::string point = Dms::compassPoint(24);    // point = 'NNE'
       *   const std::string point = Dms::compassPoint(24, 1); // point = 'N'
       */
      static std::string compassPoint(double bearing, int precision = 3);

      /**
       * Constrain degrees to range 0..360 (e.g. for bearings); -1 => 359, 361 => 1.
       *
       * @private
       * @param {number} degrees
       * @returns degrees within range 0..360.
       */
      static double wrap360(double degrees);

      /**
       * Constrain degrees to range -180..+180 (e.g. for longitude); -181 => 179, 181 => -179.
       *
       * @private
       * @param {number} degrees
       * @returns degrees within range -180..+180.
       */
      static double wrap180(double degrees);

      /**
       * Constrain degrees to range -90..+90 (e.g. for latitude); -91 => -89, 91 => 89.
       *
       * @private
       * @param {number} degrees
       * @returns degrees within range -90..+90.
       */
      static double wrap90(double degrees);

   private:
       /* Degree-minutes-seconds (& cardinal directions) separator character */
       static std::string& _separator; // U+202F = 'narrow no-break space'
   };

   template <class T, typename>
   double Dms::parse(const T& dms)
   {
      if constexpr (std::is_convertible_v<T, double>)
      {
         return dms;
      }

      // check for signed decimal degrees without NSEW, if so return it directly
      if (std::regex_search(strutil::strip(dms), std::regex("^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$")))
      {
         return std::stod(dms);
      }

      // strip off any sign or compass dir'n & split out separate d/m/s
      const std::string trim_dms = strutil::strip(dms);
      std::vector<std::string> dms_parts = strutil::split_regex(
         std::regex_replace(trim_dms, std::regex("(^-)|([NSEW]$)", std::regex_constants::icase), ""),
         "[^0-9.,]+");

      if (dms_parts[dms_parts.size() - 1].empty())
      {
         dms_parts.erase(dms_parts.begin() + (dms_parts.size() - 1));
      }

      if (dms_parts.empty())
         return NAN;

      std::vector<double> dms_parts_d{};
      for (auto elem : dms_parts)
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
         deg = dms_parts_d[0] / 1.000 + dms_parts_d[1] / 60.000 + dms_parts_d[2] / 3600.000;
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
}


#endif // DMS_H
