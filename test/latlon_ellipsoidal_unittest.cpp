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

#include "geodesy/latlon_ellipsoidal.h"
#include "geodesy/algorithm.h"

using LatLonE = geodesy::LatLonEllipsoidal;
class latlon_ellipsoidal_unittest : public testing::Test
{
protected:
   void SetUp() override
   {
      geodesy::Dms::setSeparator("");

   }

   void TearDown() override
   {

   }

   const double R = 6371e3;
};

TEST_F(latlon_ellipsoidal_unittest, constructor)
{
   EXPECT_EQ(LatLonE(51.47788, -0.00147, 17).toString(geodesy::Dms::D, 4, 2), "51.4779°N, 000.0015°W +17.00m");
}

TEST_F(latlon_ellipsoidal_unittest, parse)
{
   EXPECT_EQ(LatLonE::parse(51.47788, -0.00147).toString(), "51.4779°N, 000.0015°W");
   EXPECT_EQ(LatLonE::parse("51°28′40″N, 000°00′05″W", 17).toString(), "51.4778°N, 000.0014°W");
   EXPECT_EQ(LatLonE::parse(51.47788, -0.00147).toString(), "51.4779°N, 000.0015°W");
   EXPECT_EQ(LatLonE::parse(51.47788, -0.00147, 99).toString(geodesy::Dms::D, 4, 0), "51.4779°N, 000.0015°W +99m");
   EXPECT_EQ(LatLonE::parse("51.47788", "-0.00147", 0).toString(), "51.4779°N, 000.0015°W");
   EXPECT_EQ(LatLonE::parse("51.47788", "-0.00147", 99).toString(geodesy::Dms::D, 4, 0), "51.4779°N, 000.0015°W +99m");
   EXPECT_EQ(LatLonE::parse("51°28.67′N, 000°00.09′E").toString(geodesy::Dms::DM), "51°28.67′N, 000°00.09′E");
   EXPECT_EQ(LatLonE::parse("51°28.67′N, 000°00.09′E", 99).toString(geodesy::Dms::DM, 2, 0), "51°28.67′N, 000°00.09′E +99m");
   EXPECT_EQ(LatLonE::parse("51°28′40″N, 000°00′05″E").toString(geodesy::Dms::DMS), "51°28′40″N, 000°00′05″E");
   EXPECT_EQ(LatLonE::parse("51°28′40″N, 000°00′05″E", 99).toString(geodesy::Dms::DMS, 0, 0), "51°28′40″N, 000°00′05″E +99m");
   EXPECT_EQ(LatLonE::parse("51.47788, -0.00147").toString(), "51.4779°N, 000.0015°W");
   EXPECT_EQ(LatLonE::parse("51.47788, -0.00147", 99).toString(geodesy::Dms::D, 4, 0), "51.4779°N, 000.0015°W +99m");
   EXPECT_EQ(LatLonE::parse("51.47788, -0.00147", 99).toString(geodesy::Dms::D, 4, 0),"51.4779°N, 000.0015°W +99m");
}

TEST_F(latlon_ellipsoidal_unittest, toString)
{
   EXPECT_EQ(LatLonE(1, -2).toString(), "01.0000°N, 002.0000°W");
   EXPECT_EQ(LatLonE(1, -2).toString(geodesy::Dms::D), "01.0000°N, 002.0000°W");
   EXPECT_EQ(LatLonE(1, -2).toString(geodesy::Dms::DM), "01°00.00′N, 002°00.00′W");
   EXPECT_EQ(LatLonE(1, -2).toString(geodesy::Dms::DMS), "01°00′00″N, 002°00′00″W");
   EXPECT_EQ(LatLonE(1, -2).toString(geodesy::Dms::D, 6), "01.000000°N, 002.000000°W");
   EXPECT_EQ(LatLonE(1, -2).toString(geodesy::Dms::DM, 4), "01°00.0000′N, 002°00.0000′W");
   EXPECT_EQ(LatLonE(1, -2).toString(geodesy::Dms::DMS, 2), "01°00′00.00″N, 002°00′00.00″W");
   EXPECT_EQ(LatLonE(1, -2).toString(geodesy::Dms::D, 6, 2), "01.000000°N, 002.000000°W +0.00m");
   EXPECT_EQ(LatLonE(1, -2).toString(geodesy::Dms::DM, 4, 2), "01°00.0000′N, 002°00.0000′W +0.00m");
   EXPECT_EQ(LatLonE(1, -2).toString(geodesy::Dms::DMS, 2, 2), "01°00′00.00″N, 002°00′00.00″W +0.00m");
   EXPECT_EQ(LatLonE(1, -2, 99).toString(geodesy::Dms::D, 6, 2), "01.000000°N, 002.000000°W +99.00m");
   EXPECT_EQ(LatLonE(1, -2, 99).toString(geodesy::Dms::DM, 4, 2), "01°00.0000′N, 002°00.0000′W +99.00m");
   EXPECT_EQ(LatLonE(1, -2, 99).toString(geodesy::Dms::DMS, 2, 2), "01°00′00.00″N, 002°00′00.00″W +99.00m");
   EXPECT_EQ(LatLonE(1, -2, -99).toString(geodesy::Dms::D, 6, 2), "01.000000°N, 002.000000°W -99.00m");
   EXPECT_EQ(LatLonE(1, -2, -99).toString(geodesy::Dms::DM, 4, 2), "01°00.0000′N, 002°00.0000′W -99.00m");
   EXPECT_EQ(LatLonE(1, -2, -99).toString(geodesy::Dms::DMS, 2, 2), "01°00′00.00″N, 002°00′00.00″W -99.00m");
   EXPECT_EQ(LatLonE(1, -2).toString(geodesy::Dms::N), "1.0000, -2.0000");
   EXPECT_EQ(LatLonE(1, -2).toString(geodesy::Dms::N, 6), "1.000000, -2.000000");
   EXPECT_EQ(LatLonE(1, -2).toString(geodesy::Dms::N, 6, 0), "1.000000, -2.000000 +0m");
   EXPECT_EQ(LatLonE(1, -2).toString(geodesy::Dms::N, 6, 2), "1.000000, -2.000000 +0.00m");
   EXPECT_EQ(LatLonE(1, -2, 99).toString(geodesy::Dms::N, 6, 2), "1.000000, -2.000000 +99.00m");
   EXPECT_EQ(LatLonE(1, -2, -99).toString(geodesy::Dms::N, 6, 2), "1.000000, -2.000000 -99.00m");
}

TEST_F(latlon_ellipsoidal_unittest, getters_setters)
{
   auto p = LatLonE(51.47788, -0.00147, 99);
   EXPECT_DOUBLE_EQ(p.lat(), (51.47788));
   EXPECT_DOUBLE_EQ(p.latitude(), (51.47788));
   EXPECT_DOUBLE_EQ(p.lon(), (-0.00147));
   EXPECT_DOUBLE_EQ(p.lng(), (-0.00147));
   EXPECT_DOUBLE_EQ(p.longitude(), (-0.00147));
   EXPECT_DOUBLE_EQ(p.height(), (99));
   p.setLat(48.8584);
   EXPECT_DOUBLE_EQ(p.lat(), 48.8584);
   p.setLatitude(48.8584);
   EXPECT_DOUBLE_EQ(p.latitude(), (48.8584));
   p.setLon(2.2945);
   EXPECT_DOUBLE_EQ(p.lon(), 2.2945);
   p.setLng(2.2945);
   EXPECT_DOUBLE_EQ(p.lng(), (2.2945));
   p.setLongitude( 2.2945);
   EXPECT_DOUBLE_EQ(p.longitude(), (2.2945));
   p.setHeight(9);
   EXPECT_DOUBLE_EQ(p.height(), (9));
}