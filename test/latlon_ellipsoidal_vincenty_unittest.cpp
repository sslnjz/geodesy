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

#include "geodesy/latlon_ellipsoidal_vincenty.h"
#include "geodesy/latlon_ellipsoidal_datum.h"
#include "geodesy/dms.h"

#include <cmath>
#include <limits>
#include <stdexcept>

using Vincenty = geodesy::LatLonEllipsoidalVincenty;

class latlon_ellipsoidal_vincenty_unittest : public testing::Test
{
protected:
   void SetUp() override
   {
      geodesy::Dms::setSeparator("");
   }
   void TearDown() override
   {

   }

public:
   const double m_circEquatorial = 40075016.686; // eslint-disable-line no-unused-vars
   const double m_circMeridional = 40007862.918;
};

TEST_F(latlon_ellipsoidal_vincenty_unittest, example)
{
   EXPECT_EQ("", geodesy::Dms::separator());
   EXPECT_NEAR(Vincenty(50.06632, -5.71475).distanceTo({58.64402, -3.07009}), 969954.166, 0.001);
   EXPECT_NEAR(Vincenty(50.06632, -5.71475).initialBearingTo({58.64402, -3.07009}), 9.1419, 0.0001);
   EXPECT_NEAR(Vincenty(50.06632, -5.71475).finalBearingTo({58.64402, -3.07009}), 11.2972, 0.0001);
   EXPECT_NEAR(Vincenty(-37.95103, 144.42487).finalBearingOn(54972.271, 306.86816), 307.1736, 0.0001);

   EXPECT_EQ(Vincenty(-37.95103, 144.42487).destinationPoint(54972.271, 306.86816).toString(), "37.6528°S, 143.9265°E");
   EXPECT_EQ(Vincenty(50.06632, -5.71475).intermediatePointTo({58.64402, -3.07009}, 0.5).toString(), "54.3639°N, 004.5304°W");
}

TEST_F(latlon_ellipsoidal_vincenty_unittest, Uk)
{
   const Vincenty le(50.06632, -5.71475);
   const Vincenty jog(58.64402, -3.07009);
   const auto dist = 969954.166, brngInit = 9.1418775, brngFinal = 11.2972204;
   EXPECT_NEAR(le.distanceTo(jog), dist, 0.001);
   EXPECT_NEAR(le.initialBearingTo(jog), brngInit, 0.0000001);
   EXPECT_NEAR(le.finalBearingTo(jog), brngFinal, 0.0000001);
   EXPECT_EQ(le.destinationPoint(dist, brngInit).toString(geodesy::Dms::D), jog.toString(geodesy::Dms::D));
   EXPECT_NEAR(le.finalBearingOn(dist, brngInit), brngFinal, 0.0000001);
   EXPECT_EQ(le.intermediatePointTo(jog, 0), le);
   EXPECT_TRUE(le.intermediatePointTo(jog, 1).equals(jog));
}

TEST_F(latlon_ellipsoidal_vincenty_unittest, inverse_returns_nan_bearings_for_coincident_points)
{
   const Vincenty point(50.06632, -5.71475);

   EXPECT_NEAR(point.distanceTo(point), 0.0, 1e-9);
   EXPECT_TRUE(std::isnan(point.initialBearingTo(point)));
   EXPECT_TRUE(std::isnan(point.finalBearingTo(point)));
}

TEST_F(latlon_ellipsoidal_vincenty_unittest, inverse_handles_exact_antipodal_meridional_path)
{
   const Vincenty start(0.0, 0.0);
   const Vincenty antipode(0.0, 180.0);

   EXPECT_NEAR(start.distanceTo(antipode), 20003931.459, 0.001);
   EXPECT_NEAR(start.initialBearingTo(antipode), 0.0, 1e-12);
   EXPECT_NEAR(start.finalBearingTo(antipode), 180.0, 1e-12);
}

TEST_F(latlon_ellipsoidal_vincenty_unittest, inverse_reports_nan_for_near_antipodal_non_convergence)
{
   const Vincenty start(0.0, 0.0);
   const Vincenty nearAntipode(0.5, 179.7);

   EXPECT_TRUE(std::isnan(start.distanceTo(nearAntipode)));
   EXPECT_TRUE(std::isnan(start.initialBearingTo(nearAntipode)));
   EXPECT_TRUE(std::isnan(start.finalBearingTo(nearAntipode)));
}

TEST_F(latlon_ellipsoidal_vincenty_unittest, direct_returns_same_point_and_nan_bearing_for_zero_distance)
{
   const Vincenty start(12.345678, 98.765432);

   const auto destination = start.destinationPoint(0.0, 45.0);

   EXPECT_TRUE(destination.equals(start));
   EXPECT_TRUE(std::isnan(start.finalBearingOn(0.0, 45.0)));
}

TEST_F(latlon_ellipsoidal_vincenty_unittest, direct_rejects_non_finite_distance_and_bearing)
{
   const Vincenty start(0.0, 0.0);
   const double infinity = std::numeric_limits<double>::infinity();

   EXPECT_THROW(start.destinationPoint(std::numeric_limits<double>::quiet_NaN(), 90.0), std::invalid_argument);
   EXPECT_THROW(start.destinationPoint(infinity, 90.0), std::invalid_argument);
   EXPECT_THROW(start.destinationPoint(1000.0, std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
   EXPECT_THROW(start.destinationPoint(1000.0, infinity), std::invalid_argument);
}

TEST_F(latlon_ellipsoidal_vincenty_unittest, direct_rejects_starting_point_above_ellipsoid)
{
   const Vincenty elevated(51.47788, -0.00147, 10.0);

   EXPECT_THROW(elevated.destinationPoint(1000.0, 90.0), std::invalid_argument);
   EXPECT_THROW(elevated.finalBearingOn(1000.0, 90.0), std::invalid_argument);
}

TEST_F(latlon_ellipsoidal_vincenty_unittest, direct_normalizes_destination_across_antimeridian)
{
   const Vincenty start(0.0, 179.9);

   const auto destination = start.destinationPoint(50000.0, 90.0);

   EXPECT_NEAR(destination.lat(), 0.0, 1e-12);
   EXPECT_NEAR(destination.lon(), -179.6508423579403, 1e-12);
   EXPECT_NEAR(start.finalBearingOn(50000.0, 90.0), 90.0, 1e-12);
}

TEST_F(latlon_ellipsoidal_vincenty_unittest, direct_handles_near_antipodal_distance)
{
   const Vincenty start(0.0, 0.0);

   const auto destination = start.destinationPoint(19900000.0, 45.0);

   EXPECT_NEAR(destination.lat(), 0.5572358170221987, 1e-12);
   EXPECT_NEAR(destination.lon(), 179.01990444725766, 1e-12);
   EXPECT_NEAR(start.finalBearingOn(19900000.0, 45.0), 134.99730824394538, 1e-12);
}
