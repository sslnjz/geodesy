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
#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <cmath>
#include <iomanip>
#include <limits>
#include <string>
#include <sstream>
#include <optional>

#include "algorithm.h"

namespace geodesy
{
   /**
    * Functions for manipulating generic 3-d vectors.
    *
    * Functions return vectors as return results, so that operations can be chained.
    *
    * @example
    *   const v = v1.cross(v2).dot(v3) // ≡ v1xv2⋅v3
    */
   class vector3d
   {
   public:

      /**
       * Creates a 3-d vector.
       *
       * @param {number} x - X component of vector.
       * @param {number} y - Y component of vector.
       * @param {number} z - Z component of vector.
       *
       * @example
       *   v = new vector3d(0.267, 0.535, 0.802);
       */
      constexpr vector3d() noexcept;
      vector3d(double x, double y, double z);

      virtual ~vector3d() = default;

      [[nodiscard]] constexpr inline double x() const noexcept;
      [[nodiscard]] constexpr inline double y() const noexcept;
      [[nodiscard]] constexpr inline double z() const noexcept;
      [[nodiscard]] inline double& rx() noexcept;
      [[nodiscard]] inline double& ry() noexcept;
      [[nodiscard]] inline double& rz() noexcept;


      /**
       * Length (magnitude or norm) of ‘this’ vector.
       *
       * @returns {number} Magnitude of this vector.
       */
      [[nodiscard]] inline double length() const;

      /**
       * Adds supplied vector to ‘this’ vector.
       *
       * @param   {Vector3d} v - Vector to be added to this vector.
       * @returns {Vector3d} Vector representing sum of this and v.
       */
      [[nodiscard]] vector3d plus(const vector3d& v) const
      {
         return vector3d(vx + v.x(), vy + v.y(), vz + v.z());
      }

      /**
       * Subtracts supplied vector from ‘this’ vector.
       *
       * @param   {Vector3d} v - Vector to be subtracted from this vector.
       * @returns {Vector3d} Vector representing difference between this and v.
       */
      [[nodiscard]] vector3d minus(const vector3d& v) const
      {
         return vector3d(vx - v.x(), vy - v.y(), vz - v.z());
      }

      /**
       * Multiplies ‘this’ vector by a scalar value.
       *
       * @param   {number}   x - Factor to multiply this vector by.
       * @returns {Vector3d} Vector scaled by x.
       */
      [[nodiscard]] vector3d times(double x) const
      {
         if (std::isnan(x)) throw std::invalid_argument("invalid scalar value：");
         return vector3d(vx * x, vy * x, vz * x);
      }

      /**
       * Divides ‘this’ vector by a scalar value.
       *
       * @param   {number}   x - Factor to divide this vector by.
       * @returns {Vector3d} Vector divided by x.
       */
      [[nodiscard]] vector3d dividedBy(double x) const
      {
         if (std::isnan(x) || essentiallyEqual(x, 0.0))
            throw std::invalid_argument("invalid scalar value");

         return vector3d(vx / x, vy / x, vz / x);
      }

      /**
       * Multiplies ‘this’ vector by the supplied vector using dot (scalar) product.
       *
       * @param   {Vector3d} v - Vector to be dotted with this vector.
       * @returns {number}   Dot product of ‘this’ and v.
       */
      [[nodiscard]] double dot(const vector3d& v) const
      {
         return vx * v.vx + vy * v.vy + vz * v.vz;
      }

      /**
       * Negates a vector to point in the opposite direction.
       *
       * @returns {Vector3d} Negated vector.
       */
      [[nodiscard]] vector3d negate() const
      {
         return vector3d(-vx, -vy, -vz);
      }


      /**
       * Multiplies ‘this’ vector by the supplied vector using cross (vector) product.
       *
       * @param   {Vector3d} v - Vector to be crossed with this vector.
       * @returns {Vector3d} Cross product of ‘this’ and v.
       */
      [[nodiscard]] vector3d cross(const vector3d& v) const
      {
         const double x = vy * v.vz - vz * v.vy;
         const double y = vz * v.vx - vx * v.vz;
         const double z = vx * v.vy - vy * v.vx;
         return vector3d(x, y, z);
      }

      /**
       * Normalizes a vector to its unit vector
       * – if the vector is already unit or is zero magnitude, this is a no-op.
       *
       * @returns {Vector3d} Normalised version of this vector.
       */
      [[nodiscard]] vector3d unit() const
      {
         const double norm = length();
         if (std::fabs(norm - 1) < std::numeric_limits<double>::epsilon()) 
         {
            return *this;
         }

         if (std::fabs(norm - 0) < std::numeric_limits<double>::epsilon()) 
         {
            return *this;
         }

         return vector3d(vx / norm, vy / norm, vz / norm);
      }

      /**
       * Calculates the angle between ‘this’ vector and supplied vector atan2(|p₁×p₂|, p₁·p₂) (or if
       * (extra-planar) ‘n’ supplied then atan2(n·p₁×p₂, p₁·p₂).
       *
       * @param   {Vector3d} v - Vector whose angle is to be determined from ‘this’ vector.
       * @param   {Vector3d} [n] - Plane normal: if supplied, angle is signed +ve if this->v is
       *                     clockwise looking along n, -ve in opposite direction.
       * @returns {number}   Angle (in radians) between this vector and supplied vector (in range 0..π
       *                     if n not supplied, range -π..+π if n supplied).
       */
      double angleTo(const vector3d& v, const std::optional<vector3d>& n = std::nullopt) const
      {
         // q.v. stackoverflow.com/questions/14066933#answer-16544330, but n·p₁×p₂ is numerically
         // ill-conditioned, so just calculate sign to apply to |p₁×p₂|
         // if n·p₁×p₂ is -ve, negate |p₁×p₂|
         const int sign = n == std::nullopt || cross(v).dot(*n) >= 0 ? 1 : -1;
         const double sintheta = cross(v).length() * sign;
         const double costheta = dot(v);
         return std::atan2(sintheta, costheta);
      }

      /**
       * Rotates ‘this’ point around an axis by a specified angle.
       *
       * @param   {Vector3d} axis - The axis being rotated around.
       * @param   {number}   angle - The angle of rotation (in degrees).
       * @returns {Vector3d} The rotated point.
       */
      vector3d rotateAround(const vector3d& axis, double angle) const
      {
         const double theta = toRadians(angle);
         // en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
         // en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
         const vector3d p = unit();
         const vector3d a = axis.unit();
         const double s = std::sin(theta);
         const double c = std::cos(theta);
         const double t = 1 - c;
         const double x = a.vx, y = a.vy, z = a.vz;
         const double r[3][3] = 
         { // rotation matrix for rotation about supplied axis
            { t * x * x + c, t * x * y - s * z, t * x * z + s * y } ,
            { t * x * y + s * z, t * y * y + c, t * y * z - s * x },
            { t * x * z - s * y, t * y * z + s * x, t * z * z + c },
         };
         // multiply r × p
         const double rp[3] = {
            r[0][0] * p.vx + r[0][1] * p.vy + r[0][2] * p.vz,
            r[1][0] * p.vx + r[1][1] * p.vy + r[1][2] * p.vz,
            r[2][0] * p.vx + r[2][1] * p.vy + r[2][2] * p.vz,
         };
         return vector3d(rp[0], rp[1], rp[2]);
         // https://en.wikipedia.org/wiki/Rodrigues'_rotation_formula
      }

      /**
       * String representation of vector.
       *
       * @param   {number} [dp=3] - Number of decimal places to be used.
       * @returns {string} Vector represented as [x,y,z].
       */
      [[nodiscard]] virtual std::string toString(int dp = 3) const
      {
         std::stringstream ss;
         ss << std::fixed << std::setprecision(dp) << "[" << vx << ", " << vy << ", " << vz << "]";
         return ss.str();
      }

      /**
       * Adds supplied vector to ‘this’ vector, same as add(v).
       *
       * @param   {Vector3d} v - Vector to be added to this vector.
       * @returns {Vector3d} Vector representing sum of this and v.
       */
      constexpr inline vector3d& operator+=(const vector3d& v)
      {
         vx += v.vx;
         vy += v.vy;
         vz += v.vz;
         return *this;
      }

      /**
       * Subtracts supplied vector from ‘this’ vector, same as minus(v).
       *
       * @param   {Vector3d} v - Vector to be subtracted from this vector.
       * @returns {Vector3d} Vector representing difference between this and v.
       */
      constexpr inline vector3d& operator-=(const vector3d& v)
      {
         vx -= v.vx;
         vy -= v.vy;
         vz -= v.vz;
         return *this;
      }

      /**
       * Multiplies ‘this’ vector by a scalar value, same as time(x).
       *
       * @param   {number}   x - Factor to multiply this vector by.
       * @returns {Vector3d} Vector scaled by x.
       */
      constexpr inline vector3d& operator*=(double x)
      {
         vx *= x;
         vy *= x;
         vz *= x;
         return *this;
      }

      /**
       * Multiplies ‘this’ vector by a scalar value, same as dividedBy(x).
       *
       * @param   {number}   x - Factor to multiply this vector by.
       * @returns {Vector3d} Vector scaled by x.
       */
      constexpr inline vector3d& operator/=(double x)
      {
         vx /= x;
         vy /= x;
         vz /= x;
         return *this;
      }

      constexpr inline bool operator==(const vector3d& v) const
      {
         if (this != &v)
         {
            return essentiallyEqual(vx, v.vx) &&
               essentiallyEqual(vy, v.vy) &&
               essentiallyEqual(vz, v.vz);
         }
         return true;
      }


      friend inline vector3d operator+(const vector3d& v1, const vector3d& v2)
      { return vector3d(v1.vx + v2.vx, v1.vy + v2.vy, v1.vz + v2.vz);}
      friend inline vector3d operator-(const vector3d& v1, const vector3d& v2)
      { return vector3d(v1.vx - v2.vx, v1.vy - v2.vy, v1.vz - v2.vz);}
      friend inline vector3d operator*(const vector3d& v, double x)
      { return vector3d(v.vx * x, v.vy * x, v.vz * x);}
      friend inline vector3d operator*(double x, const vector3d& v)
      { return vector3d(v.vx * x, v.vy * x, v.vz * x);}
      friend vector3d operator/(const vector3d& v, double x)
      { return vector3d(v.vx / x, v.vy / x, v.vz / x);}
      friend inline vector3d operator+(const vector3d& v)
      { return v;}
      friend inline vector3d operator-(const vector3d& v)
      { return vector3d(v.x(), v.y(), v.z());}

   private:
      double vx;
      double vy;
      double vz;
   };


   constexpr vector3d::vector3d() noexcept : vx(0), vy(0), vz(0) { }
   vector3d::vector3d(double x, double y, double z) : vx(x), vy(y), vz(z) 
   {
        if(std::isnan(x) || std::isnan(y) || std::isnan(z)){
            std::stringstream ss;
            ss << "invalid vector[" << x << ", " << y << ", " << z << "]";
            throw std::runtime_error(ss.str());
        }
   }

   constexpr inline double vector3d::x() const noexcept
   {
      return vx;
   }

   constexpr inline double vector3d::y() const noexcept
   {
      return vy;
   }

   constexpr inline double vector3d::z() const noexcept
   {
      return vz;
   }

   inline double& vector3d::rx() noexcept
   {
      return vx;
   }

   inline double& vector3d::ry() noexcept
   {
      return vy;
   }

   inline double& vector3d::rz() noexcept
   {
      return vz;
   }

   inline double vector3d::length() const
   {
      return std::sqrt(vx * vx + vy * vy + vz * vz);
   }
}

#endif // VECTOR3D_H
