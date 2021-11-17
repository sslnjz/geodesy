
#ifndef ALGORITHM_H
#define ALGORITHM_H

// π
constexpr auto π = (3.141592653589793116);

namespace geodesy
{

   /**
    * Conversion degrees to Radians.
    *
    * @param {number} degrees
    * @returns Radians.
    */
   static double toRadians(double degrees)
   {
       return degrees * π / 180.000;
   }

   /**
    * Conversion radians to degrees.
    *
    * @param {number} radians
    * @returns degrees.
    */
   static double toDegrees(double radians)
   {
       return radians * 180 / π;
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

}

#endif // ALGORITHM_H
