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

#include <sstream>
#include <cmath>

#include "geodesy/latlon_spherical.h"
#include "geodesy/algorithm.h"

class latlon_spherical_unittest : public testing::Test
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

TEST_F(latlon_spherical_unittest, constructor)
{
   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical(52.205, 0.119).distanceTo(geodesy::LatLonSpherical(48.857, 2.351))), "404279");
   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical(52.205, 0.119).distanceTo(geodesy::LatLonSpherical(48.857, 2.351), 3959), 1), "251.2");
   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical(52.205, 0.119).initialBearingTo(geodesy::LatLonSpherical(48.857, 2.351)), 1), "156.2");
   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical(52.205, 0.119).finalBearingTo(geodesy::LatLonSpherical(48.857, 2.351)), 1), "157.9");
   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical(53.2611, -0.7972).crossTrackDistanceTo(geodesy::LatLonSpherical(53.3206, -1.7297), geodesy::LatLonSpherical(53.1887, 0.1334)), 1), "-307.5");
   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical(53.2611, -0.7972).alongTrackDistanceTo(geodesy::LatLonSpherical(53.3206, -1.7297), geodesy::LatLonSpherical(53.1887, 0.1334))), "62331");
   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical(51.127, 1.338).rhumbDistanceTo(geodesy::LatLonSpherical(50.964, 1.853))), "40308");
   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical(51.127, 1.338).rhumbBearingTo(geodesy::LatLonSpherical(50.964, 1.853)), 1), "116.7");

   EXPECT_EQ(geodesy::LatLonSpherical(51.127, 1.338).rhumbDestinationPoint(40300, 116.7).toString(), "50.9642°N, 001.8530°E");
   EXPECT_EQ(geodesy::LatLonSpherical(52.205, 0.119).midpointTo(geodesy::LatLonSpherical(48.857, 2.351)).toString(), "50.5363°N, 001.2746°E");
   EXPECT_EQ(geodesy::LatLonSpherical(52.205, 0.119).intermediatePointTo(geodesy::LatLonSpherical(48.857, 2.351), 0.25).toString(), "51.3721°N, 000.7073°E");
   EXPECT_EQ(geodesy::LatLonSpherical(51.47788, -0.00147).destinationPoint(7794, 300.7).toString(), "51.5136°N, 000.0983°W");
   EXPECT_EQ(geodesy::LatLonSpherical::intersection(geodesy::LatLonSpherical(51.8853, 0.2545), 108.547, geodesy::LatLonSpherical(49.0034, 2.5735), 32.435).toString(), "50.9078°N, 004.5084°E");
   EXPECT_EQ(geodesy::LatLonSpherical(51.127, 1.338).rhumbMidpointTo(geodesy::LatLonSpherical(50.964, 1.853)).toString(), "51.0455°N, 001.5957°E");

   std::vector vec = { geodesy::LatLonSpherical(0, 0), geodesy::LatLonSpherical(1, 0), geodesy::LatLonSpherical(0, 1) };
   EXPECT_EQ(geodesy::toExponential(geodesy::LatLonSpherical::areaOf(vec), 2), "6.18e+09");
}

TEST_F(latlon_spherical_unittest, alternate_point_formats)
{
   const auto cambg = geodesy::LatLonSpherical(52.205, 0.119);
   EXPECT_EQ(geodesy::toExponential(cambg.distanceTo(geodesy::LatLonSpherical(48.857, 2.351)), 3), "4.043e+05");
   EXPECT_EQ(geodesy::toExponential(cambg.distanceTo(geodesy::LatLonSpherical("48.857", "2.351")), 3), "4.043e+05");
   EXPECT_EQ(geodesy::toExponential(cambg.distanceTo(geodesy::LatLonSpherical("48°51′25.2″N", "002°21′03.6″E")), 3), "4.043e+05");
}

TEST_F(latlon_spherical_unittest, dist_brng_dest)
{
   const auto cambg = geodesy::LatLonSpherical(52.205, 0.119), paris = geodesy::LatLonSpherical(48.857, 2.351);
   EXPECT_EQ(geodesy::toExponential(cambg.distanceTo(paris), 3), "4.043e+05");
   EXPECT_EQ(geodesy::toFixed(cambg.distanceTo(paris, 3959), 1), "251.2");
   EXPECT_EQ(geodesy::toFixed(cambg.initialBearingTo(paris), 1), "156.2");
   EXPECT_EQ(geodesy::toFixed(cambg.finalBearingTo(paris), 1), "157.9");
   EXPECT_TRUE(std::isnan(cambg.initialBearingTo(cambg)));
   EXPECT_TRUE(std::isnan(cambg.finalBearingTo(cambg)));
   EXPECT_EQ(geodesy::toFixed(paris.initialBearingTo(cambg), 1), "337.9");
   EXPECT_EQ(cambg.midpointTo(paris).toString(), "50.5363°N, 001.2746°E");
   EXPECT_EQ(cambg.intermediatePointTo(paris, 0.25).toString(), "51.3721°N, 000.7073°E");
   EXPECT_EQ(cambg.intermediatePointTo(cambg, 0.25).toString(), "52.2050°N, 000.1190°E");
   const auto greenwich = geodesy::LatLonSpherical(51.47788, -0.00147);
   const double dist = 7794, brng = 300.7;
   EXPECT_EQ(greenwich.destinationPoint(dist, brng).toString(), "51.5136°N, 000.0983°W");
   EXPECT_EQ(greenwich.destinationPoint(dist, brng, 6371e3).toString(), "51.5136°N, 000.0983°W");
}

TEST_F(latlon_spherical_unittest, dist_brng_dest_fails)
{
   const auto cambg = geodesy::LatLonSpherical(52.205, 0.119), paris = geodesy::LatLonSpherical(48.857, 2.351);
   EXPECT_THROW(cambg.distanceTo({ "xxx", "xxx" }), std::invalid_argument);
};

TEST_F(latlon_spherical_unittest, intersection)
{
   const double N = 0, E = 90, S = 180, W = 270;
   EXPECT_EQ(geodesy::LatLonSpherical::intersection(geodesy::LatLonSpherical(0, 1), N, geodesy::LatLonSpherical(1, 0), E).toString(), "00.9998°N, 001.0000°E");
   EXPECT_EQ(geodesy::LatLonSpherical::intersection(geodesy::LatLonSpherical(1, 0), E, geodesy::LatLonSpherical(0, 1), N).toString(), "00.9998°N, 001.0000°E");
   EXPECT_THROW(geodesy::LatLonSpherical::intersection(geodesy::LatLonSpherical(2, 1), N, geodesy::LatLonSpherical(1, 0), E), std::runtime_error);
   EXPECT_THROW(geodesy::LatLonSpherical::intersection(geodesy::LatLonSpherical(0, 1), N, geodesy::LatLonSpherical(1, 0), W), std::runtime_error);
   EXPECT_THROW(geodesy::LatLonSpherical::intersection(geodesy::LatLonSpherical(1, 0), W, geodesy::LatLonSpherical(0, 1), N), std::runtime_error);
   EXPECT_THROW(geodesy::LatLonSpherical::intersection(geodesy::LatLonSpherical(0, 1), S, geodesy::LatLonSpherical(1, 0), E), std::runtime_error);
   EXPECT_THROW(geodesy::LatLonSpherical::intersection(geodesy::LatLonSpherical(1, 0), E, geodesy::LatLonSpherical(0, 1), S), std::runtime_error);
   EXPECT_EQ(geodesy::LatLonSpherical::intersection(geodesy::LatLonSpherical(0, 1), S, geodesy::LatLonSpherical(1, 0), W).toString(), "00.9998°S, 179.0000°W");
   EXPECT_EQ(geodesy::LatLonSpherical::intersection(geodesy::LatLonSpherical(1, 0), W, geodesy::LatLonSpherical(0, 1), S).toString(), "00.9998°S, 179.0000°W");

   EXPECT_THROW(geodesy::LatLonSpherical::intersection(geodesy::LatLonSpherical(0, 1), N, geodesy::LatLonSpherical(1, 90), E), std::runtime_error);
   EXPECT_EQ(geodesy::LatLonSpherical::intersection(geodesy::LatLonSpherical(0, 1), N, geodesy::LatLonSpherical(1, 92), E).toString(), "00.0175°N, 179.0000°W");

   EXPECT_EQ(geodesy::LatLonSpherical::intersection(geodesy::LatLonSpherical(1, 1), N, geodesy::LatLonSpherical(1, 1), E).toString(), "01.0000°N, 001.0000°E");

   const auto stn = geodesy::LatLonSpherical(51.8853, 0.2545), cdg = geodesy::LatLonSpherical(49.0034, 2.5735);
   EXPECT_EQ(geodesy::LatLonSpherical::intersection(stn, 108.547, cdg, 32.435).toString(), "50.9078°N, 004.5084°E");

   EXPECT_EQ(geodesy::LatLonSpherical::intersection(geodesy::LatLonSpherical(51, 0), 120, geodesy::LatLonSpherical(50, 0), 60).toString(), "50.4921°N, 001.3612°E");

   //TODO: Windows OK, MAC failed, needs check
   EXPECT_EQ(geodesy::LatLonSpherical::intersection(geodesy::LatLonSpherical(-77.6966041375563, 18.28125), 179.99999999999994, geodesy::LatLonSpherical(89, 180), 180).toString(), "90.0000°S, 163.9902°W");
}

TEST_F(latlon_spherical_unittest, cross_track_along_track)
{
   EXPECT_EQ(geodesy::toExponential(geodesy::LatLonSpherical(10, 1).crossTrackDistanceTo(geodesy::LatLonSpherical(0, 0), geodesy::LatLonSpherical(0, 2)), 3), "-1.112e+06");

   const auto bradwell = geodesy::LatLonSpherical(53.3206, -1.7297);
   const auto dunham   = geodesy::LatLonSpherical(53.2611, -0.7972);
   const auto partney  = geodesy::LatLonSpherical(53.1887, 0.1334);
   EXPECT_EQ(geodesy::toPrecision(dunham.crossTrackDistanceTo(bradwell, partney), 4), "-307.5");
   EXPECT_EQ(geodesy::toPrecision(dunham.alongTrackDistanceTo(bradwell, partney), 4), "6.233e+04");

   EXPECT_EQ(geodesy::toPrecision(geodesy::LatLonSpherical(1, 1).crossTrackDistanceTo(geodesy::LatLonSpherical(0, 0), geodesy::LatLonSpherical(0, 2)), 4), "-1.112e+05");
   EXPECT_EQ(geodesy::toPrecision(geodesy::LatLonSpherical(-1, 1).crossTrackDistanceTo(geodesy::LatLonSpherical(0, 0), geodesy::LatLonSpherical(0, 2)), 4), "1.112e+05");
   EXPECT_EQ(geodesy::toPrecision(geodesy::LatLonSpherical(-1, -1).crossTrackDistanceTo(geodesy::LatLonSpherical(0, 0), geodesy::LatLonSpherical(0, 2)), 4), "1.112e+05");
   EXPECT_EQ(geodesy::toPrecision(geodesy::LatLonSpherical(1, -1).crossTrackDistanceTo(geodesy::LatLonSpherical(0, 0), geodesy::LatLonSpherical(0, 2)), 4), "-1.112e+05"); // eslint-disable-line space-in-parens

   EXPECT_EQ(geodesy::toPrecision(geodesy::LatLonSpherical(1, 1).alongTrackDistanceTo(geodesy::LatLonSpherical(0, 0), geodesy::LatLonSpherical(0, 2)), 4), "1.112e+05"); // eslint-disable-line space-in-parens
   EXPECT_EQ(geodesy::toPrecision(geodesy::LatLonSpherical(-1, 1).alongTrackDistanceTo(geodesy::LatLonSpherical(0, 0), geodesy::LatLonSpherical(0, 2)), 4), "1.112e+05");
   EXPECT_EQ(geodesy::toPrecision(geodesy::LatLonSpherical(-1, -1).alongTrackDistanceTo(geodesy::LatLonSpherical(0, 0), geodesy::LatLonSpherical(0, 2)), 4), "-1.112e+05");
   EXPECT_EQ(geodesy::toPrecision(geodesy::LatLonSpherical(1, -1).alongTrackDistanceTo(geodesy::LatLonSpherical(0, 0), geodesy::LatLonSpherical(0, 2)), 4), "-1.112e+05"); // eslint-disable-line space-in-parens

   EXPECT_EQ(geodesy::LatLonSpherical(10, 0).crossTrackDistanceTo(geodesy::LatLonSpherical(10, 0), geodesy::LatLonSpherical(0, 2)), 0);
   EXPECT_EQ(geodesy::LatLonSpherical(10, 0).alongTrackDistanceTo(geodesy::LatLonSpherical(10, 0), geodesy::LatLonSpherical(0, 2)), 0);
}

TEST_F(latlon_spherical_unittest, misc)
{
   EXPECT_EQ(geodesy::LatLonSpherical(0, 0).maxLatitude(0), 90);
   EXPECT_EQ(geodesy::LatLonSpherical(0, 0).maxLatitude(1), 89);
   EXPECT_EQ(geodesy::LatLonSpherical(0, 0).maxLatitude(90), 0);
   
   const auto& [lon1, lon2] = geodesy::LatLonSpherical::crossingParallels(geodesy::LatLonSpherical(0, 0), geodesy::LatLonSpherical(60, 30), 30);
   EXPECT_EQ(geodesy::LatLonSpherical(30, lon1).toString(geodesy::Dms::DMS), "30°00′00″N, 009°35′39″E");
   EXPECT_EQ(geodesy::LatLonSpherical(30, lon2).toString(geodesy::Dms::DMS), "30°00′00″N, 170°24′21″E");

   const std::pair<double, double> empty = {};
   EXPECT_EQ(geodesy::LatLonSpherical::crossingParallels(geodesy::LatLonSpherical(0, 0), geodesy::LatLonSpherical(10, 60), 60), empty);
   EXPECT_EQ(geodesy::LatLonSpherical::crossingParallels(geodesy::LatLonSpherical(0, 0), geodesy::LatLonSpherical(0, 0), 0), empty);
}

TEST_F(latlon_spherical_unittest, area_polygon_based)
{
   const auto epsilon = std::numeric_limits<double>::epsilon();
   const auto R = 6371e3;

   std::vector polyTriangle = { geodesy::LatLonSpherical(1, 1), geodesy::LatLonSpherical(2, 1), geodesy::LatLonSpherical(1, 2) };
   std::vector polySquareCw = { geodesy::LatLonSpherical(1, 1), geodesy::LatLonSpherical(2, 1), geodesy::LatLonSpherical(2, 2), geodesy::LatLonSpherical(1, 2) };
   std::vector polySquareCcw = { geodesy::LatLonSpherical(1, 1), geodesy::LatLonSpherical(1, 2), geodesy::LatLonSpherical(2, 2), geodesy::LatLonSpherical(2, 1) };
   std::vector polyOctant = { geodesy::LatLonSpherical(0, epsilon), geodesy::LatLonSpherical(90, 0), geodesy::LatLonSpherical(0, 90 - epsilon) };
   std::vector polyOctantS = { geodesy::LatLonSpherical(-epsilon, epsilon), geodesy::LatLonSpherical(90, 0), geodesy::LatLonSpherical(-epsilon, 90 - epsilon) };
   std::vector polyQuadrant = { geodesy::LatLonSpherical(epsilon, epsilon), geodesy::LatLonSpherical(90, epsilon), geodesy::LatLonSpherical(epsilon, 180 - epsilon), geodesy::LatLonSpherical(epsilon, 90) };
   std::vector polyHemiE = { geodesy::LatLonSpherical(epsilon, epsilon), geodesy::LatLonSpherical(90 - epsilon, 0), geodesy::LatLonSpherical(90 - epsilon, 180), geodesy::LatLonSpherical(epsilon, 180), geodesy::LatLonSpherical(-epsilon, 180), geodesy::LatLonSpherical(-90 + epsilon, 180), geodesy::LatLonSpherical(-90 + epsilon, 0), geodesy::LatLonSpherical(-epsilon, epsilon) };
   std::vector polyPole = { geodesy::LatLonSpherical(89, 0), geodesy::LatLonSpherical(89, 120), geodesy::LatLonSpherical(89, -120) };
   std::vector polyConcave = { geodesy::LatLonSpherical(1, 1), geodesy::LatLonSpherical(5, 1), geodesy::LatLonSpherical(5, 3), geodesy::LatLonSpherical(1, 3), geodesy::LatLonSpherical(3, 2) };

   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical::areaOf(polyTriangle)), "6181527888");
   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical::areaOf(polyTriangle, 6371e3)), "6181527888");
   polyTriangle.emplace_back(polyTriangle[0]);
   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical::areaOf(polyTriangle)), "6181527888");
   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical::areaOf(polySquareCw)), "12360230987");
   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical::areaOf(polySquareCcw)), "12360230987");
   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical::areaOf(polyOctant), 1), geodesy::toFixed((pi * R * R / 2), 1));
   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical::areaOf(polyOctantS), 1), geodesy::toFixed((pi * R * R / 2), 1));
   EXPECT_EQ(geodesy::LatLonSpherical::areaOf(polyQuadrant), pi * R * R);
   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical::areaOf(polyHemiE), 1), geodesy::toFixed((2 * pi * R * R), 1));
   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical::areaOf(polyPole)), "16063139192");
   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical::areaOf(polyConcave)), "74042699236");
}

TEST_F(latlon_spherical_unittest, Ed_Williams)
{ // www.edwilliams.org/avform.htm
   const auto lax = geodesy::LatLonSpherical(geodesy::Dms::parse("33° 57′N"), geodesy::Dms::parse("118° 24′W"));
   const auto jfk = geodesy::LatLonSpherical(geodesy::Dms::parse("40° 38′N"), geodesy::Dms::parse("073° 47′W"));
   const auto r = 180 * 60 / pi; // earth radius in nautical miles
   EXPECT_EQ(geodesy::toPrecision(lax.distanceTo(jfk, r), 4), "2144");
   EXPECT_EQ(geodesy::toPrecision(lax.initialBearingTo(jfk), 2), "66");

   std::string dm = lax.intermediatePointTo(jfk, 100.0/2144.0).toString(geodesy::Dms::DM, 0);
   EXPECT_EQ(lax.intermediatePointTo(jfk, 100.0/2144.0).toString(geodesy::Dms::DM, 0), "34°37′N, 116°33′W");

   const auto d = geodesy::LatLonSpherical(geodesy::Dms::parse("34:30N"), geodesy::Dms::parse("116:30W"));
   EXPECT_EQ(geodesy::toPrecision(d.crossTrackDistanceTo(lax, jfk, r), 5), "7.4523");
   EXPECT_EQ(geodesy::toPrecision(d.alongTrackDistanceTo(lax, jfk, r), 5), "99.588");
   EXPECT_EQ(lax.intermediatePointTo(jfk, 0.4).toString(geodesy::Dms::DM, 3), "38°40.167′N, 101°37.570′W");
   const auto reo = geodesy::LatLonSpherical(geodesy::Dms::parse("42.600N"), geodesy::Dms::parse("117.866W"));
   const auto bke = geodesy::LatLonSpherical(geodesy::Dms::parse("44.840N"), geodesy::Dms::parse("117.806W"));
   EXPECT_EQ(geodesy::LatLonSpherical::intersection(reo, 51, bke, 137).toString(geodesy::Dms::D, 3), "43.572°N, 116.189°W");
}

TEST_F(latlon_spherical_unittest, rhumb_lines)
{
   const auto dov = geodesy::LatLonSpherical(51.127, 1.338), cal = geodesy::LatLonSpherical(50.964, 1.853);
   EXPECT_EQ(geodesy::toPrecision(dov.rhumbDistanceTo(cal), 4), "4.031e+04");
   EXPECT_EQ(geodesy::toPrecision(dov.rhumbDistanceTo(cal, 6371e3), 4), "4.031e+04");
   EXPECT_EQ(geodesy::toFixed(geodesy::LatLonSpherical(1, -179).rhumbDistanceTo(geodesy::LatLonSpherical(1, 179)), 6), 
      geodesy::toFixed(geodesy::LatLonSpherical(1, 1).rhumbDistanceTo(geodesy::LatLonSpherical(1, -1)), 6));
   EXPECT_EQ(geodesy::toFixed(dov.rhumbBearingTo(cal), 1), "116.7");
   EXPECT_EQ(geodesy::LatLonSpherical(1, -179).rhumbBearingTo(geodesy::LatLonSpherical(1, 179)), 270);
   EXPECT_EQ(geodesy::LatLonSpherical(1, 179).rhumbBearingTo(geodesy::LatLonSpherical(1, -179)), 90);
   EXPECT_TRUE(std::isnan(dov.rhumbBearingTo(dov)));
   EXPECT_EQ(dov.rhumbDestinationPoint(40310, 116.7).toString(), "50.9641°N, 001.8531°E");
   EXPECT_EQ(dov.rhumbDestinationPoint(40310, 116.7, 6371e3).toString(), "50.9641°N, 001.8531°E");
   EXPECT_EQ(geodesy::LatLonSpherical(1, 1).rhumbDestinationPoint(111178, 90).toString(), "01.0000°N, 002.0000°E");
   EXPECT_EQ(geodesy::LatLonSpherical(1, 179).rhumbDestinationPoint(222356, 90).toString(), "01.0000°N, 179.0000°W");
   EXPECT_EQ(geodesy::LatLonSpherical(1, -179).rhumbDestinationPoint(222356, 270).toString(), "01.0000°N, 179.0000°E");
   EXPECT_EQ(dov.rhumbMidpointTo(cal).toString(), "51.0455°N, 001.5957°E");
   EXPECT_EQ(geodesy::LatLonSpherical(1, -179).rhumbMidpointTo(geodesy::LatLonSpherical(1, 178)).toString(), "01.0000°N, 179.5000°E");
}