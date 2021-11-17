
#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <cmath>
#include <iomanip>
#include <limits>
#include <string>
#include <sstream>

#include "algorithm.hpp"

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
      constexpr vector3d(double x, double y, double z) noexcept;

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
       * Multiplies ‘this’ vector by the supplied vector using dot (scalar) product.
       *
       * @param   {Vector3d} v - Vector to be dotted with this vector.
       * @returns {number}   Dot product of ‘this’ and v.
       */
      [[nodiscard]] double dot(const vector3d& v) const
      {
         return xv * v.xv + yv * v.yv + zv * v.zv;
      }

      /**
       * Negates a vector to point in the opposite direction.
       *
       * @returns {Vector3d} Negated vector.
       */
      [[nodiscard]] vector3d negate() const
      {
         return vector3d(-xv, -yv, -zv);
      }


      /**
       * Multiplies ‘this’ vector by the supplied vector using cross (vector) product.
       *
       * @param   {Vector3d} v - Vector to be crossed with this vector.
       * @returns {Vector3d} Cross product of ‘this’ and v.
       */
      [[nodiscard]] vector3d cross(const vector3d& v) const
      {
         const double x = yv * v.zv - zv * v.yv;
         const double y = zv * v.xv - xv * v.zv;
         const double z = xv * v.yv - yv * v.xv;
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

         return vector3d(xv / norm, yv / norm, zv / norm);
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
      double angleTo(const vector3d& v, vector3d* n = nullptr) const
      {
         // q.v. stackoverflow.com/questions/14066933#answer-16544330, but n·p₁×p₂ is numerically
         // ill-conditioned, so just calculate sign to apply to |p₁×p₂|
         // if n·p₁×p₂ is -ve, negate |p₁×p₂|
         const int sign = n == nullptr || cross(v).dot(*n) >= 0 ? 1 : -1;
         const double sinθ = cross(v).length() * sign;
         const double cosθ = dot(v);
         return std::atan2(sinθ, cosθ);
      }

      /**
       * Rotates ‘this’ point around an axis by a specified angle.
       *
       * @param   {Vector3d} axis - The axis being rotated around.
       * @param   {number}   angle - The angle of rotation (in degrees).
       * @returns {Vector3d} The rotated point.
       */
      vector3d rotateAround(const vector3d& axis, double angle)
      {
         const double θ = toRadians(angle);
         // en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
         // en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
         const vector3d p = unit();
         const vector3d a = axis.unit();
         const double s = std::sin(θ);
         const double c = std::cos(θ);
         const double t = 1 - c;
         const double x = a.xv, y = a.yv, z = a.zv;
         const double r[3][3] = 
         { // rotation matrix for rotation about supplied axis
            { t * x * x + c, t * x * y - s * z, t * x * z + s * y } ,
            { t * x * y + s * z, t * y * y + c, t * y * z - s * x },
            { t * x * z - s * y, t * y * z + s * x, t * z * z + c },
         };
         // multiply r × p
         const double rp[3] = {
            r[0][0] * p.xv + r[0][1] * p.yv + r[0][2] * p.zv,
            r[1][0] * p.xv + r[1][1] * p.yv + r[1][2] * p.zv,
            r[2][0] * p.xv + r[2][1] * p.yv + r[2][2] * p.zv,
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
      [[nodiscard]] std::string toString(int dp = 3) const
      {
         std::stringstream ss;
         ss << "[";
         ss << std::setprecision(dp) << xv << ", ";
         ss << std::setprecision(dp) << yv << ", ";
         ss << std::setprecision(dp) << zv;
         ss << "]";
         return ss.str();
      }

      /**
       * Adds supplied vector to ‘this’ vector.
       *
       * @param   {Vector3d} v - Vector to be added to this vector.
       * @returns {Vector3d} Vector representing sum of this and v.
       */
      friend constexpr inline vector3d operator+(const vector3d& v1, const vector3d& v2)
      {
         return vector3d(v1.xv + v2.xv, v1.yv + v2.yv, v1.zv + v2.zv);
      }

      /**
       * Subtracts supplied vector from ‘this’ vector.
       *
       * @param   {Vector3d} v - Vector to be subtracted from this vector.
       * @returns {Vector3d} Vector representing difference between this and v.
       */
      friend constexpr inline vector3d operator-(const vector3d& v1, const vector3d& v2)
      {
         return vector3d(v1.xv - v2.xv, v1.yv - v2.yv, v1.zv - v2.zv);
      }

      /**
       * Multiplies ‘this’ vector by a scalar value.
       *
       * @param   {number}   x - Factor to multiply this vector by.
       * @returns {Vector3d} Vector scaled by x.
       */
      friend constexpr vector3d operator*(const vector3d& v1, const vector3d& v2)
      {
         return vector3d(v1.xv * v2.xv, v1.yv * v2.xv, v1.zv * v2.xv);
      }

      /**
       * Divides ‘this’ vector by a scalar value.
       *
       * @param   {number}   x - Factor to divide this vector by.
       * @returns {Vector3d} Vector divided by x.
       */
      friend constexpr vector3d operator/(const vector3d& v1, const vector3d& v2)
      {
         return vector3d(v1.xv / v2.xv, v1.yv / v2.xv, v1.zv / v2.zv);
      }

   private:
      double xv;
      double yv;
      double zv;
   };


   constexpr vector3d::vector3d() noexcept : xv(0), yv(0), zv(0) { }
   constexpr vector3d::vector3d(double x, double y, double z) noexcept : xv(x), yv(y), zv(z) { }

   constexpr inline double vector3d::x() const noexcept
   {
      return xv;
   }

   constexpr inline double vector3d::y() const noexcept
   {
      return yv;
   }

   constexpr inline double vector3d::z() const noexcept
   {
      return zv;
   }

   inline double& vector3d::rx() noexcept
   {
      return xv;
   }

   inline double& vector3d::ry() noexcept
   {
      return yv;
   }

   inline double& vector3d::rz() noexcept
   {
      return zv;
   }

   inline double vector3d::length() const
   {
      return std::sqrt(xv * xv + yv * yv + zv * zv);
   }

}

#endif // VECTOR3D_H
