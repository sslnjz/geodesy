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
#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <limits>
#include <cmath>
// π
constexpr auto pi = (3.141592653589793116);
constexpr auto epsilon = std::numeric_limits<double>::epsilon();

namespace geodesy
{

   /**
    * Conversion degrees to Radians.
    *
    * @param {number} degrees
    * @returns Radians.
    */
   [[maybe_unused]] static double toRadians(double degrees)
   {
       return degrees * pi / 180.000;
   }

   /**
    * Conversion radians to degrees.
    *
    * @param {number} radians
    * @returns degrees.
    */
   [[maybe_unused]] static double toDegrees(double radians)
   {
       return radians * 180 / pi;
   }

   /**
    * @brief Demonstrates the use of machine epsilon to compare floating-point values for equality
    * @tparam T float type
    * @param x number to compare
    * @param y number to be compare
    * @param ulp units in the last place
    * @return true if two float number almost equal with ulp
   */
   template<class T>
   typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
      almost_equal(T x, T y, int ulp)
   {
      // the machine epsilon has to be scaled to the magnitude of the values used
      // and multiplied by the desired precision in ULPs (units in the last place)
      return std::fabs(x - y) <= std::numeric_limits<T>::epsilon() * std::fabs(x + y) * ulp
         // unless the result is subnormal
         || std::fabs(x - y) < std::numeric_limits<T>::min();
   }

   template <class T>
   inline int sign(const T& z)
   {
      return (z == 0) ? 0 : std::signbit(z) ? -1 : 1;
   }

   /**
    * Conversion factor metres to kilometres.
    * Conversion factors; 1000 * LatLon.metresToKm gives 1.
    */
   static inline double getMetresToKm()
   {
      return 1.000 / 1000.000;
   }

   static inline double getMetresToMiles()
   {
      return 1.000 / 1609.344;
   }

   static inline double getMetresToNauticalMiles()
   {
      return 1.000 / 1852.000;
   }

}

#endif // ALGORITHM_H
