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

#include "geodesy/cartesian.h"
#include "geodesy/latlon_ellipsoidal.h"
#include "geodesy/algorithm.h"

#include <limits>
#include <stdexcept>

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

namespace
{
constexpr double metreTolerance = 1e-3;
constexpr double angularTolerance = 1e-9;
}

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

TEST_F(latlon_ellipsoidal_unittest, wgs84_ellipsoid_constants)
{
   const auto wgs84 = LatLonE::ellipsoids().WGS84;

   EXPECT_DOUBLE_EQ(wgs84.a, 6378137.0);
   EXPECT_NEAR(wgs84.b, 6356752.314245, 1e-6);
   EXPECT_DOUBLE_EQ(wgs84.f, 1.0 / 298.257223563);
}

TEST_F(latlon_ellipsoidal_unittest, converts_equator_origin_to_cartesian)
{
   const auto cartesian = LatLonE(0.0, 0.0).toCartesian();

   EXPECT_NEAR(cartesian.x(), 6378137.0, metreTolerance);
   EXPECT_NEAR(cartesian.y(), 0.0, metreTolerance);
   EXPECT_NEAR(cartesian.z(), 0.0, metreTolerance);
}

TEST_F(latlon_ellipsoidal_unittest, converts_height_above_equator_to_cartesian)
{
   const auto cartesian = LatLonE(0.0, 0.0, 1000.0).toCartesian();

   EXPECT_NEAR(cartesian.x(), 6379137.0, metreTolerance);
   EXPECT_NEAR(cartesian.y(), 0.0, metreTolerance);
   EXPECT_NEAR(cartesian.z(), 0.0, metreTolerance);
}

TEST_F(latlon_ellipsoidal_unittest, converts_greenwich_reference_point_to_cartesian)
{
   const auto cartesian = LatLonE(51.47788, -0.00147, 17.0).toCartesian();

   EXPECT_NEAR(cartesian.x(), 3980584.835, metreTolerance);
   EXPECT_NEAR(cartesian.y(), -102.127, metreTolerance);
   EXPECT_NEAR(cartesian.z(), 4966843.366, metreTolerance);
}

TEST_F(latlon_ellipsoidal_unittest, converts_cartesian_equator_to_latlon)
{
   const auto point = geodesy::Cartesian(6378137.0, 0.0, 0.0).toLatLon();

   EXPECT_NEAR(point.lat(), 0.0, angularTolerance);
   EXPECT_NEAR(point.lon(), 0.0, angularTolerance);
   EXPECT_NEAR(point.height(), 0.0, metreTolerance);
}

TEST_F(latlon_ellipsoidal_unittest, converts_poles_between_latlon_and_cartesian)
{
   const auto northCartesian = LatLonE(90.0, 0.0).toCartesian();
   EXPECT_NEAR(northCartesian.x(), 0.0, metreTolerance);
   EXPECT_NEAR(northCartesian.y(), 0.0, metreTolerance);
   EXPECT_NEAR(northCartesian.z(), 6356752.314245, metreTolerance);

   const auto northPole = geodesy::Cartesian(0.0, 0.0, 6356752.314245).toLatLon();
   EXPECT_NEAR(northPole.lat(), 90.0, angularTolerance);
   EXPECT_NEAR(northPole.lon(), 0.0, angularTolerance);
   EXPECT_NEAR(northPole.height(), 0.0, metreTolerance);

   const auto southPole = geodesy::Cartesian(0.0, 0.0, -6356752.314245).toLatLon();
   EXPECT_NEAR(southPole.lat(), -90.0, angularTolerance);
   EXPECT_NEAR(southPole.lon(), 0.0, angularTolerance);
   EXPECT_NEAR(southPole.height(), 0.0, metreTolerance);
}

TEST_F(latlon_ellipsoidal_unittest, round_trips_latlon_height_through_cartesian)
{
   const auto source = LatLonE(50.7978, 4.3592, 123.45);
   const auto point = source.toCartesian().toLatLon();

   EXPECT_NEAR(point.lat(), source.lat(), angularTolerance);
   EXPECT_NEAR(point.lon(), source.lon(), angularTolerance);
   EXPECT_NEAR(point.height(), source.height(), metreTolerance);
}

TEST_F(latlon_ellipsoidal_unittest, rejects_non_finite_cartesian_coordinates)
{
   const double infinity = std::numeric_limits<double>::infinity();

   EXPECT_THROW(static_cast<void>(geodesy::Cartesian(infinity, 0.0, 0.0).toLatLon()), std::invalid_argument);
   EXPECT_THROW(static_cast<void>(geodesy::Cartesian(0.0, infinity, 0.0).toLatLon()), std::invalid_argument);
   EXPECT_THROW(static_cast<void>(geodesy::Cartesian(0.0, 0.0, infinity).toLatLon()), std::invalid_argument);
}

TEST_F(latlon_ellipsoidal_unittest, rejects_invalid_cartesian_origin_and_ellipsoid)
{
   EXPECT_THROW(static_cast<void>(geodesy::Cartesian(0.0, 0.0, 0.0).toLatLon()), std::domain_error);

   const geodesy::Ellipsoid invalidEllipsoid { 6378137.0, 6356752.314245, 1.0 };
   EXPECT_THROW(static_cast<void>(geodesy::Cartesian(6378137.0, 0.0, 0.0).toLatLon(invalidEllipsoid)),
                std::invalid_argument);
}

TEST_F(latlon_ellipsoidal_unittest, rejects_non_finite_geodetic_height)
{
   const double infinity = std::numeric_limits<double>::infinity();

   EXPECT_THROW(static_cast<void>(LatLonE(0.0, 0.0, infinity).toCartesian()), std::invalid_argument);
}

TEST_F(latlon_ellipsoidal_unittest, round_trips_non_wgs84_ellipsoid_through_cartesian)
{
   const auto datum = LatLonE::datums().OSGB36;
   const auto source = LatLonE(52.65798, 1.71605, 24.7, datum);
   const auto point = source.toCartesian().toLatLon(datum.ellipsoid);

   EXPECT_NEAR(point.lat(), source.lat(), angularTolerance);
   EXPECT_NEAR(point.lon(), source.lon(), angularTolerance);
   EXPECT_NEAR(point.height(), source.height(), metreTolerance);
}
