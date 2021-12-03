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

#include <gtest/gtest.h>

#include "geodesy/vector3d.h"

TEST(vector3d_unittest, Constructor)
{
   const geodesy::vector3d v3d(0.267, 0.535, 0.802);
   EXPECT_DOUBLE_EQ(v3d.x(), 0.267);
   EXPECT_DOUBLE_EQ(v3d.y(), 0.535);
   EXPECT_DOUBLE_EQ(v3d.z(), 0.802);
}

TEST(vector3d_unittest, methods)
{
   const auto v123 = geodesy::vector3d(1, 2, 3);
   const auto v321 = geodesy::vector3d(3, 2, 1);
   EXPECT_EQ(v123.plus(v321), geodesy::vector3d(4, 4, 4));
   EXPECT_EQ(v123.minus(v321), geodesy::vector3d(-2, 0, 2));
   EXPECT_EQ(v123.times(2), geodesy::vector3d(2, 4, 6));
   EXPECT_EQ(v123.dividedBy(2), geodesy::vector3d(0.5, 1, 1.5));
   EXPECT_EQ(v123.dot(v321), 10);
   EXPECT_EQ(v123.cross(v321), geodesy::vector3d(-4, 8, -4));
   EXPECT_EQ(v123.negate(), geodesy::vector3d(-1, -2, -3));
   EXPECT_EQ(v123.length(), 3.7416573867739413);
   EXPECT_EQ(v123.unit().toString(), "[0.267,0.535,0.802]");
   EXPECT_EQ(geodesy::toFixed(geodesy::toDegrees(v123.angleTo(v321)), 3), "44.415");
   EXPECT_EQ(geodesy::toFixed(geodesy::toDegrees(v123.angleTo(v321, v123.cross(v321))), 3), "44.415");
   EXPECT_EQ(geodesy::toFixed(geodesy::toDegrees(v123.angleTo(v321, v321.cross(v123))), 3), "-44.415");
   EXPECT_EQ(geodesy::toFixed(geodesy::toDegrees(v123.angleTo(v321, v123)), 3), "44.415");
   EXPECT_EQ(v123.rotateAround(geodesy::vector3d(0, 0, 1), 90).toString(), "[-0.535,0.267,0.802]");
   EXPECT_EQ(v123.toString(), "[1.000,2.000,3.000]");
   EXPECT_EQ(v123.toString(6), "[1.000000,2.000000,3.000000]");
}

TEST(vector3d_unittest, operators)
{
   auto v123 = geodesy::vector3d(1, 2, 3);
   const auto v321 = geodesy::vector3d(3, 2, 1);
   EXPECT_EQ((v123 += v321), geodesy::vector3d(4, 4, 4));
   v123 = geodesy::vector3d(1, 2, 3);
   EXPECT_EQ(v123 -= v321, geodesy::vector3d(-2, 0, 2));
   v123 = geodesy::vector3d(1, 2, 3);
   EXPECT_EQ(v123 *= 2,    geodesy::vector3d(2, 4, 6));
   v123 = geodesy::vector3d(1, 2, 3);
   EXPECT_EQ(v123 /= 2,    geodesy::vector3d(0.5, 1, 1.5));
}