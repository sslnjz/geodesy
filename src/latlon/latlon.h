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
#ifndef LATLON_H
#define LATLON_H

#include <string>

#include "dms.h"

namespace geodesy
{
   /**
    * Base LatLon class
    */
   class LatLon
   {
   public:
      LatLon() noexcept;
      LatLon(double lat, double lon);
      LatLon(const std::string& lat, const std::string& lon);
      virtual ~LatLon();

      LatLon(const LatLon& rhs) noexcept;
      LatLon(LatLon&& rhs) noexcept;
      LatLon& operator=(const LatLon& rhs) noexcept;
      LatLon& operator=(LatLon&& rhs) noexcept;

      [[nodiscard]] double lat() const;
      [[nodiscard]] double latitude() const;

      template<class T>
      using ENABLE = std::enable_if_t<std::is_arithmetic_v<T> || std::is_convertible_v<T, std::string>>;

      template<class T, typename = ENABLE<T>>
      void setLat(const T& lat);
      template<class T, typename = ENABLE<T>>
      void setLatitude(const T& lat);

      /**
       * Parses a latitude/longitude point from a variety of formats.
       *
       * Latitude & longitude (in degrees) can be supplied as two separate parameters, as a single
       * comma-separated lat/lon string, or as a single object with { lat, lon } or GeoJSON properties.
       *
       * The latitude/longitude values may be numeric or strings; they may be signed decimal or
       * deg-min-sec (hexagesimal) suffixed by compass direction (NSEW); a variety of separators are
       * accepted. Examples -3.62, '3 37 12W', '3°37′12″W'.
       *
       * Thousands/decimal separators must be comma/dot; use Dms.fromLocale to convert locale-specific
       * thousands/decimal separators.
       *
       * @param   {number|string|Object} lat|latlon - Latitude (in degrees) or comma-separated lat/lon or lat/lon object.
       * @param   {number|string}        [lon]      - Longitude (in degrees).
       * @returns {LatLon} Latitude/longitude point.
       * @throws  {TypeError} Invalid point.
       *
       * @example
       *   const p1 = LatLon.parse(52.205, 0.119);                                    // numeric pair (≡ new LatLon)
       *   const p2 = LatLon.parse('52.205', '0.119');                                // numeric string pair (≡ new LatLon)
       *   const p3 = LatLon.parse('52.205, 0.119');                                  // single string numerics
       *   const p4 = LatLon.parse('52°12′18.0″N', '000°07′08.4″E');                  // DMS pair
       *   const p5 = LatLon.parse('52°12′18.0″N, 000°07′08.4″E');                    // single string DMS
       *   const p6 = LatLon.parse({ lat: 52.205, lon: 0.119 });                      // { lat, lon } object numeric
       *   const p7 = LatLon.parse({ lat: '52°12′18.0″N', lng: '000°07′08.4″E' });    // { lat, lng } object DMS
       *   const p8 = LatLon.parse({ type: 'Point', coordinates: [ 0.119, 52.205] }); // GeoJSON
       */
      static LatLon parse(double lat, double lon);
      static LatLon parse(const std::string& dms);
      static LatLon parse(const std::string& lat, const std::string& lon);

      /**
       * Longitude in degrees east from international reference meridian (including aliases lon, lng,
       * longitude): can be set as numeric or hexadecimal (deg-min-sec); returned as numeric.
       */
      [[nodiscard]] double lon() const;
      [[nodiscard]] double lng() const;
      [[nodiscard]] double longitude() const;

      template<class T, typename = ENABLE<T>>
      void setLon(const T& lon);
      template<class T, typename = ENABLE<T>>
      void setLng(const T& lon);
      template<class T, typename = ENABLE<T>>
      void setLongitude(const T& lon);

      /**
       * Returns a string representation of ‘this’ point, formatted as degrees, degrees+minutes, or
       * degrees+minutes+seconds.
       *
       * @param   {enum} [eFormat=D] - Format point as 'D', 'DM', 'DMS', or 'N' for signed numeric.
       * @returns {string} Comma-separated formatted latitude/longitude.
       *
       * @example
       *   const greenwich = new LatLon(51.47788, -0.00147);
       *   const d = greenwich.toString();                        // 51.4779°N, 000.0015°W
       *   const dms = greenwich.toString('dms', 2);              // 51°28′40.37″N, 000°00′05.29″W
       *   const [lat, lon] = greenwich.toString('n').split(','); // 51.4779, -0.0015
       */
      [[nodiscard]] std::string toString(Dms::eFormat e = Dms::D, std::optional<int> dp = std::nullopt) const;

      /**
       * Converts ‘this’ point to a GeoJSON object string.
       *
       * @returns {string} this point as a GeoJSON ‘Point’ string.
       *    { type: "Point", coordinates: [ lon, lat ] }
       */
      [[nodiscard]] std::string toGeoJSON() const;

      /**
       * Checks if another point is equal to ‘this’ point.
       *
       * @param   {LatLon} p1 - Point to be compared against p2.
       * @param   {LatLon} p2 - Point to be compared against p1.
       * @returns {bool}   True if points have identical latitude and longitude values.
       *
       */
      friend bool operator==(const LatLon& p1, const LatLon& p2);

      /**
       * Checks if another point is NOT equal to ‘this’ point.
       *
       * @param   {LatLon} p1 - Point to be compared against p2.
       * @param   {LatLon} p2 - Point to be compared against p1.
       * @returns {bool}   True if points have different latitude or longitude values.
       *
       */
      friend bool operator!=(const LatLon& p1, const LatLon& p2);

   protected:
      double m_lat; // Latitude in degrees north from equator
      double m_lon; // Longitude in degrees east from international reference meridian
   };

   template <class T, typename>
   void LatLon::setLat(const T& lat)
   {
      if constexpr (std::is_arithmetic_v<T>)
      {
         m_lat = Dms::wrap90(lat);
      }
      else
      {
         m_lat = Dms::wrap90(Dms::parse(lat));
      }
      
      if (std::isnan(m_lat))
         throw std::invalid_argument("invalid lat");
   }

   template <class T, typename>
   void LatLon::setLatitude(const T& lat)
   {
      setLat(lat);
   }

   template <class T, typename>
   void LatLon::setLon(const T& lon)
   {
      if constexpr (std::is_arithmetic_v<T>)
      {
         m_lon = Dms::wrap180(lon);
      }
      else
      {
         m_lon = Dms::wrap180(Dms::parse(lon));
      }
      
      if (std::isnan(m_lon))
         throw std::invalid_argument("invalid lon");
   }

   template <class T, typename>
   void LatLon::setLng(const T& lon)
   {
      setLon(lon);
   }

   template <class T, typename>
   void LatLon::setLongitude(const T& lon)
   {
      setLon(lon);
   }
}

#endif //LATLON_H
