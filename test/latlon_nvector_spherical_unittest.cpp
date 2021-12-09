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

#include "geodesy/latlon_nvector_spherical.h"
#include "geodesy/nvector_spherical.h"
#include "geodesy/algorithm.h"
#include "geodesy/strutil.h"

namespace 
{
   using LatLonN = geodesy::LatLonNvectorSpherical;
   using NS = geodesy::NvectorSpherical;
   class latlon_nvector_spherical_unittest : public testing::Test
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


   TEST_F(latlon_nvector_spherical_unittest, examples)
   {
      EXPECT_EQ(LatLonN(52.205, 0.119).toString(), "52.2050°N, 000.1190°E");
      EXPECT_EQ(LatLonN(45, 45).toNvector().toString(), "[0.500, 0.500, 0.707]");
      EXPECT_EQ(LatLonN(53.3206, -1.7297).greatCircle(96.0).toString(), "[-0.794, 0.129, 0.594]");
      EXPECT_EQ(geodesy::toFixed(LatLonN(52.205, 0.119).distanceTo(LatLonN(48.857, 2.351))), "404279");
      EXPECT_EQ(geodesy::toFixed(LatLonN(52.205, 0.119).distanceTo(LatLonN(48.857, 2.351), 3959), 1), "251.2");
      EXPECT_EQ(geodesy::toFixed(LatLonN(52.205, 0.119).initialBearingTo(LatLonN(48.857, 2.351)), 1), "156.2");
      EXPECT_EQ(geodesy::toFixed(LatLonN(52.205, 0.119).finalBearingTo(LatLonN(48.857, 2.351)), 1), "157.9");
      EXPECT_EQ(LatLonN(52.205, 0.119).midpointTo(LatLonN(48.857, 2.351)).toString(), "50.5363°N, 001.2746°E");
      EXPECT_EQ(LatLonN(52.205, 0.119).intermediatePointTo(LatLonN(48.857, 2.351), 0.25).toString(), "51.3721°N, 000.7073°E");
      EXPECT_EQ(LatLonN(52.205, 0.119).intermediatePointOnChordTo(LatLonN(48.857, 2.351), 0.25).toString(), "51.3723°N, 000.7072°E");
      EXPECT_EQ(LatLonN(51.47788, -0.00147).destinationPoint(7794, 300.7).toString(), "51.5136°N, 000.0983°W");
      EXPECT_EQ(LatLonN::intersection(LatLonN(51.8853, 0.2545), 108.547, LatLonN(49.0034, 2.5735), 32.435).toString(), "50.9078°N, 004.5084°E");
      EXPECT_EQ(geodesy::toFixed(LatLonN(53.2611, -0.7972).crossTrackDistanceTo(LatLonN(53.3206, -1.7297), LatLonN(53.1887, 0.1334)), 1), "-307.5");
      EXPECT_EQ(geodesy::toFixed(LatLonN(53.2611, -0.7972).alongTrackDistanceTo(LatLonN(53.3206, -1.7297), LatLonN(53.1887, 0.1334))), "62331");
      EXPECT_EQ(LatLonN(51.0, 1.9).nearestPointOnSegment(LatLonN(51.0, 1.0), LatLonN(51.0, 2.0)).toString(), "51.0004°N, 001.9000°E");
      EXPECT_EQ(LatLonN(51.0, 2.1).nearestPointOnSegment(LatLonN(51.0, 1.0), LatLonN(51.0, 2.0)).toString(), "51.0000°N, 002.0000°E");
      EXPECT_EQ(LatLonN(10, -140).nearestPointOnSegment(LatLonN(0, 20), LatLonN(0, 40)).toString(), "00.0000°N, 020.0000°E");
      EXPECT_TRUE(LatLonN(52, 1).isWithinExtent(LatLonN(51, 1), LatLonN(52, 2)));
      EXPECT_FALSE(LatLonN(51, 0).isWithinExtent(LatLonN(51, 1), LatLonN(52, 2)));
      EXPECT_FALSE(LatLonN(51, 0).isWithinExtent(LatLonN(51, 1), LatLonN(51, 1)));
      EXPECT_EQ(LatLonN::triangulate(LatLonN(50.7175, 1.65139), 333.3508, LatLonN(50.9250, 1.7094), 310.1414).toString(), "51.1297°N, 001.3214°E");
      EXPECT_EQ(LatLonN::trilaterate(LatLonN(0, 0), 157e3, LatLonN(0, 1), 111e3, LatLonN(1, 0), 111e3).toString(), "00.9985°N, 000.9986°E");
      std::vector bounds({ LatLonN(45, 1), LatLonN(45, 2), LatLonN(46, 2), LatLonN(46, 1) });
      EXPECT_TRUE(LatLonN(45.1, 1.1).isEnclosedBy(bounds));
      std::vector A({ LatLonN(0, 0), LatLonN(1, 0), LatLonN(0, 1) });
      EXPECT_EQ(geodesy::toExponential(LatLonN::areaOf(A), 2), "6.18e+09");
      std::vector B({ LatLonN(0, 0), LatLonN(0, 1), LatLonN(1, 0) });
      EXPECT_EQ(geodesy::toExponential(LatLonN::areaOf(B), 2), "6.18e+09");
      EXPECT_EQ(LatLonN::meanOf({ LatLonN(1, 1), LatLonN(4, 2), LatLonN(1, 3) }).toString(), "02.0001°N, 002.0000°E");
      EXPECT_TRUE(LatLonN(52.205, 0.119).equals(LatLonN(52.205, 0.119)));
      const auto greenwich = LatLonN(51.47788, -0.00147);
      EXPECT_EQ(greenwich.toString(), "51.4779°N, 000.0015°W");
      EXPECT_EQ(greenwich.toString(geodesy::Dms::DMS, 2), "51°28′40.37″N, 000°00′05.29″W");
      EXPECT_EQ(greenwich.toString(geodesy::Dms::N), "51.4779, -0.0015");
   }


   TEST_F(latlon_nvector_spherical_unittest, examples_NvectorSpherical)
   {
      EXPECT_EQ(NS(0.5000, 0.5000, 0.7071).toString(), "[0.500, 0.500, 0.707]");
      EXPECT_EQ(NS(0.5000, 0.5000, 0.7071).toLatLon().toString(geodesy::Dms::D, 1), "45.0°N, 045.0°E");
   }

   TEST_F(latlon_nvector_spherical_unittest, constructor_with_strings)
   {
      EXPECT_EQ(geodesy::toFixed(LatLonN("52.205", "0.119").distanceTo(LatLonN("48.857", "2.351"))), "404279");
   }

   TEST_F(latlon_nvector_spherical_unittest, constructor_fail)
   {
      EXPECT_THROW(LatLonN("x", "x"), std::invalid_argument);
   }

   TEST_F(latlon_nvector_spherical_unittest, getters_setters)
   {
        auto camb = LatLonN(0, 0);
        camb.setLat("52°12′18″N");
        camb.setLon("000°07′08″E");
        EXPECT_NEAR(camb.lat(), 52.205, 0.001);
        EXPECT_NEAR(camb.latitude(), 52.205, 0.001);
        EXPECT_NEAR(camb.lon(), 0.119, 0.001);
        EXPECT_NEAR(camb.lng(), 0.119, 0.001);
        EXPECT_NEAR(camb.longitude(), 0.119, 0.001);
        camb.setLat(52.205);
        camb.setLon(0.119);
        EXPECT_DOUBLE_EQ(camb.lat(), 52.205);
        EXPECT_DOUBLE_EQ(camb.latitude(), 52.205);
        EXPECT_DOUBLE_EQ(camb.lon(), 0.119);
        EXPECT_DOUBLE_EQ(camb.lng(), 0.119);
        EXPECT_DOUBLE_EQ(camb.longitude(), 0.119);
    }

   
   TEST_F(latlon_nvector_spherical_unittest, setters_fail)
   {
       auto camb = LatLonN(0, 0);
       EXPECT_THROW(camb.setLat("xxx"), std::invalid_argument);
       EXPECT_THROW(camb.setLatitude("xxx"), std::invalid_argument);
       EXPECT_THROW(camb.setLon("xxx"), std::invalid_argument);
       EXPECT_THROW(camb.setLng("xxx"), std::invalid_argument);
       EXPECT_THROW(camb.setLongitude("xxx"), std::invalid_argument);
   }

   TEST_F(latlon_nvector_spherical_unittest, toString)
   {
        auto btTower = LatLonN(51.521470, -0.138833);
        EXPECT_EQ(btTower.toString(), "51.5215°N, 000.1388°W");
        EXPECT_EQ(btTower.toString(geodesy::Dms::D, 6), "51.521470°N, 000.138833°W");
        EXPECT_EQ(btTower.toString(geodesy::Dms::DM, 4), "51°31.2882′N, 000°08.3300′W");
        EXPECT_EQ(btTower.toString(geodesy::Dms::DMS, 2), "51°31′17.29″N, 000°08′19.80″W");
        EXPECT_EQ(btTower.toString(geodesy::Dms::N), "51.5215, -0.1388");
        EXPECT_EQ(btTower.toString(geodesy::Dms::N, 6), "51.521470, -0.138833");
   }

   TEST_F(latlon_nvector_spherical_unittest, great_circle)
   {
       EXPECT_EQ(LatLonN(53.3206, -1.7297).greatCircle(96.0).toString(), "[-0.794, 0.129, 0.594]");
       EXPECT_EQ(LatLonN(53.3206, -1.7297).toNvector().greatCircle(96.0).toString(), "[-0.794, 0.129, 0.594]");
   }

   TEST_F(latlon_nvector_spherical_unittest, dist_brng_dest)
   {
      const auto cambg = LatLonN(52.205, 0.119), paris = LatLonN(48.857, 2.351);
      EXPECT_EQ(geodesy::toPrecision(cambg.distanceTo(paris), 4), "4.043e+05");
      EXPECT_EQ(geodesy::toPrecision(cambg.distanceTo(paris, 3959), 4), "251.2");
      EXPECT_EQ(geodesy::toFixed(cambg.initialBearingTo(paris), 1), "156.2");
      EXPECT_EQ(geodesy::toFixed(cambg.finalBearingTo(paris), 1), "157.9");
      EXPECT_TRUE(std::isnan(cambg.initialBearingTo(cambg)));
      EXPECT_TRUE(std::isnan(cambg.finalBearingTo(cambg)));
      EXPECT_EQ(geodesy::toFixed(paris.initialBearingTo(cambg), 1), "337.9");
      EXPECT_EQ(cambg.midpointTo(paris).toString(), "50.5363°N, 001.2746°E");
      EXPECT_EQ(cambg.intermediatePointTo(paris, 0.25).toString(), "51.3721°N, 000.7073°E");
      EXPECT_EQ(cambg.intermediatePointTo(cambg, 0.25).toString(), "52.2050°N, 000.1190°E");
      EXPECT_EQ(LatLonN(90, 0).intermediatePointTo(LatLonN(0, 90), 0.75).toString(), "22.5000°N, 090.0000°E");
      EXPECT_EQ(LatLonN(90, 0).intermediatePointOnChordTo(LatLonN(0, 90), 0.75).toString(), "18.4349°N, 090.0000°E");
      const auto greenwich = LatLonN(51.47788, -0.00147);
      double dist = 7794, brng = 300.7;
      EXPECT_EQ(greenwich.destinationPoint(dist, brng).toString(), "51.5136°N, 000.0983°W");
      EXPECT_EQ(greenwich.destinationPoint(dist, brng, 6371e3).toString(), "51.5136°N, 000.0983°W");
   }

   TEST_F(latlon_nvector_spherical_unittest, intersection)
   {
      const double N = 0, E = 90, S = 180, W = 270;
      EXPECT_EQ(LatLonN::intersection(LatLonN(0, 1), N, LatLonN(1, 0), E).toString(), "00.9998°N, 001.0000°E");
      EXPECT_EQ(LatLonN::intersection(LatLonN(1, 0), E, LatLonN(0, 1), N).toString(), "00.9998°N, 001.0000°E");
      EXPECT_EQ(LatLonN::intersection(LatLonN(2, 1), N, LatLonN(1, 0), E).toString(), "00.9998°S, 179.0000°W");
      EXPECT_EQ(LatLonN::intersection(LatLonN(0, 1), N, LatLonN(1, 0), W).toString(), "00.9998°S, 179.0000°W");
      EXPECT_EQ(LatLonN::intersection(LatLonN(1, 0), W, LatLonN(0, 1), N).toString(), "00.9998°S, 179.0000°W");
      EXPECT_EQ(LatLonN::intersection(LatLonN(0, 1), S, LatLonN(1, 0), E).toString(), "00.9998°S, 179.0000°W");
      EXPECT_EQ(LatLonN::intersection(LatLonN(1, 0), E, LatLonN(0, 1), S).toString(), "00.9998°S, 179.0000°W");
      EXPECT_EQ(LatLonN::intersection(LatLonN(0, 1), S, LatLonN(1, 0), W).toString(), "00.9998°S, 179.0000°W");
      EXPECT_EQ(LatLonN::intersection(LatLonN(1, 0), W, LatLonN(0, 1), S).toString(), "00.9998°S, 179.0000°W");
      EXPECT_EQ(LatLonN::intersection(LatLonN(0, 1), N, LatLonN(1, 90), E).toString(), "00.0175°S, 179.0000°W");
      EXPECT_EQ(LatLonN::intersection(LatLonN(0, 1), N, LatLonN(1, 92), E).toString(), "00.0175°N, 179.0000°W");
      EXPECT_EQ(LatLonN::intersection(LatLonN(1, 0), LatLonN(1, 3), LatLonN(2, 2), S).toString(), "01.0003°N, 002.0000°E");
      EXPECT_EQ(LatLonN::intersection(LatLonN(2, 2), S, LatLonN(1, 0), LatLonN(1, 3)).toString(), "01.0003°N, 002.0000°E");
      EXPECT_EQ(LatLonN::intersection(LatLonN(1, 0), LatLonN(1, 3), LatLonN(2, 2), N).toString(), "01.0003°S, 178.0000°W");
      EXPECT_EQ(LatLonN::intersection(LatLonN(2, 2), N, LatLonN(1, 0), LatLonN(1, 3)).toString(), "01.0003°S, 178.0000°W");
      EXPECT_EQ(LatLonN::intersection(LatLonN(1, 1), LatLonN(2, 2), LatLonN(1, 4), LatLonN(2, 3)).toString(), "02.4994°N, 002.5000°E");
      EXPECT_EQ(LatLonN::intersection(LatLonN(1, 1), N, LatLonN(1, 1), E).toString(), "01.0000°N, 001.0000°E");

      const auto stn = LatLonN(51.8853, 0.2545), cdg = LatLonN(49.0034, 2.5735);
      EXPECT_EQ(LatLonN::intersection(stn, 108.547, cdg, 32.435).toString(), "50.9078°N, 004.5084°E");
      EXPECT_EQ(LatLonN::intersection(LatLonN(51, 0), 120, LatLonN(50, 0), 60).toString(), "50.4921°N, 001.3612°E");
   }

   TEST_F(latlon_nvector_spherical_unittest, crosstrack_alongtrack)
   {
      EXPECT_EQ(geodesy::toPrecision(LatLonN(10, 1).crossTrackDistanceTo(LatLonN(0, 0), LatLonN(0, 2)), 4), "-1.112e+06");
      EXPECT_EQ(geodesy::toPrecision(LatLonN(10, 0).crossTrackDistanceTo(LatLonN(0, 0), 90), 4), "-1.112e+06");
      EXPECT_EQ(geodesy::toPrecision(LatLonN(10, 0).crossTrackDistanceTo(LatLonN(0, 0), 270), 4), "1.112e+06");

      const auto bradwell = LatLonN(53.3206, -1.7297), dunham = LatLonN(53.2611, -0.7972), partney = LatLonN(53.1887,  0.1334);
      EXPECT_EQ(geodesy::toPrecision(dunham.crossTrackDistanceTo(bradwell, partney), 4), "-307.5");
      EXPECT_EQ(geodesy::toPrecision(dunham.alongTrackDistanceTo(bradwell, partney), 4), "6.233e+04");
      EXPECT_EQ(geodesy::toPrecision(dunham.alongTrackDistanceTo(bradwell, 96.0), 4), "6.233e+04");

      EXPECT_EQ(geodesy::toPrecision(LatLonN(1, 1).crossTrackDistanceTo(LatLonN(0, 0), LatLonN(0, 2)), 4), "-1.112e+05");
      EXPECT_EQ(geodesy::toPrecision(LatLonN(-1,  1).crossTrackDistanceTo(LatLonN(0, 0), LatLonN(0, 2)), 4), "1.112e+05");
      EXPECT_EQ(geodesy::toPrecision(LatLonN(-1, -1).crossTrackDistanceTo(LatLonN(0, 0), LatLonN(0, 2)), 4), "1.112e+05");
      EXPECT_EQ(geodesy::toPrecision(LatLonN( 1, -1).crossTrackDistanceTo(LatLonN(0, 0), LatLonN(0, 2)), 4), "-1.112e+05"); // eslint-disable-line space-in-parens
      EXPECT_EQ(geodesy::toPrecision(LatLonN( 1,  1).alongTrackDistanceTo(LatLonN(0, 0), LatLonN(0, 2)), 4), "1.112e+05"); // eslint-disable-line space-in-parens
      EXPECT_EQ(geodesy::toPrecision(LatLonN(-1,  1).alongTrackDistanceTo(LatLonN(0, 0), LatLonN(0, 2)), 4), "1.112e+05");
      EXPECT_EQ(geodesy::toPrecision(LatLonN(-1, -1).alongTrackDistanceTo(LatLonN(0, 0), LatLonN(0, 2)), 4), "-1.112e+05");
      EXPECT_EQ(geodesy::toPrecision(LatLonN( 1, -1).alongTrackDistanceTo(LatLonN(0, 0), LatLonN(0, 2)), 4), "-1.112e+05"); // eslint-disable-line space-in-parens

      EXPECT_EQ(geodesy::toPrecision(LatLonN(1, 0).crossTrackDistanceTo(LatLonN(0, 0), 90), 4), "-1.112e+05");
      EXPECT_EQ(geodesy::toPrecision(LatLonN(1, 0).crossTrackDistanceTo(LatLonN(0, 0), 270), 4), "1.112e+05");

      EXPECT_EQ(LatLonN(10, 0).crossTrackDistanceTo(LatLonN(10, 0), 0), 0);
      EXPECT_EQ(LatLonN(10, 0).crossTrackDistanceTo(LatLonN(10, 0), 0), 0);
   }

   TEST_F(latlon_nvector_spherical_unittest, trilaterate)
   { // http://gis.stackexchange.com/a/415/41129

      const double d1 = 265.710701754, d2 = 234.592423446, d3 = 54.8954278262;
      const auto p1 = LatLonN(37.418436, -121.963477);
      const auto p2 = LatLonN(37.417243, -121.961889);
      const auto p3 = LatLonN(37.418692, -121.960194);
      EXPECT_EQ(LatLonN::trilaterate(p1, d1, p2, d2, p3, d3).toString(geodesy::Dms::D, 6), "37.419078°N, 121.960579°W");
      EXPECT_EQ(LatLonN::trilaterate(p1, d1, p1, d2, p1, d3), geodesy::LatLonNvectorSpherical({}));
   }

   TEST_F(latlon_nvector_spherical_unittest, area_enclosed_polygon_based) 
   {
      std::vector polyTriangle = { LatLonN(1, 1), LatLonN(2, 1), LatLonN(1, 2) };
      std::vector polySquareCw = { LatLonN(1, 1), LatLonN(2, 1), LatLonN(2, 2), LatLonN(1, 2) };
      std::vector polySquareCcw = { LatLonN(1, 1), LatLonN(1, 2), LatLonN(2, 2), LatLonN(2, 1) };
      std::vector polyOctant = { LatLonN(0, epsilon), LatLonN(90, 0), LatLonN(0, 90-epsilon) };
      std::vector polyOctantS = { LatLonN(-epsilon, epsilon), LatLonN(90, 0), LatLonN(-epsilon, 90-epsilon) };
      // const polyQuadrant = { LatLonN(epsilon, epsilon), LatLonN(90, epsilon), LatLonN(epsilon, 180-epsilon), LatLonN(epsilon, 90) };
      std::vector polyHemiE = { LatLonN(epsilon, epsilon), LatLonN(90-epsilon, 0), LatLonN(90-epsilon, 180), LatLonN(epsilon, 180), LatLonN(-epsilon, 180), LatLonN(-90+epsilon, 180), LatLonN(-90+epsilon, 0), LatLonN(-epsilon, epsilon) };
      std::vector polyGc = { LatLonN(10, 0), LatLonN(10, 90), LatLonN(0, 45) };
      std::vector polyPole = { LatLonN(89, 0), LatLonN(89, 120), LatLonN(89, -120) };
      std::vector polyPoleEdge = { LatLonN(85, 90), LatLonN(85, 0), LatLonN(85, -90) };
      std::vector polyConcave = { LatLonN(1, 1), LatLonN(5, 1), LatLonN(5, 3), LatLonN(1, 3), LatLonN(3, 2) };

      EXPECT_EQ(geodesy::toFixed(LatLonN::areaOf(polyTriangle), 0), "6181527888");
      EXPECT_EQ(geodesy::toFixed(LatLonN::areaOf(polyTriangle, 6371e3), 0), "6181527888");
      polyTriangle.emplace_back(polyTriangle[0]);
      EXPECT_EQ(geodesy::toFixed(LatLonN::areaOf(polyTriangle), 0), "6181527888");
      EXPECT_EQ(geodesy::toFixed(LatLonN::areaOf(polySquareCw), 0), "12360230987");
      EXPECT_EQ(geodesy::toFixed(LatLonN::areaOf(polySquareCcw), 0), "12360230987");
      EXPECT_EQ(LatLonN::areaOf(polyOctant), pi*R*R/2);
      EXPECT_EQ(LatLonN::areaOf(polyOctantS), pi*R*R/2);
      // TODO: fails: test('quadrant area', () => LatLon.areaOf(polyQuadrant).should.equal(π*R*R));
      EXPECT_EQ(geodesy::toFixed(LatLonN::areaOf(polyHemiE), 1), geodesy::toFixed(2*pi*R*R, 1));
      EXPECT_EQ(geodesy::toFixed(LatLonN::areaOf(polyPole)), "16063139192");
      EXPECT_EQ(geodesy::toFixed(LatLonN::areaOf(polyConcave)), "74042699236");

      EXPECT_TRUE(LatLonN(45, 1).isEnclosedBy(polyHemiE));
      // TODO: fails: test('hemisphere enclosed n', () => LatLonN(45, -1).isEnclosedBy(polyHemiE).should.be.false);
      EXPECT_TRUE(LatLonN(14, 45).isEnclosedBy(polyGc));
      EXPECT_FALSE(LatLonN(15, 45).isEnclosedBy(polyGc));
      EXPECT_TRUE(LatLonN(90, 0).isEnclosedBy(polyPole));
      EXPECT_TRUE(LatLonN(90, 0).isEnclosedBy(polyPoleEdge));
      EXPECT_TRUE(LatLonN(4, 2).isEnclosedBy(polyConcave));
      EXPECT_FALSE(LatLonN(2, 2).isEnclosedBy(polyConcave));
   }

   TEST_F(latlon_nvector_spherical_unittest, Ed_Williams)
   { // www.edwilliams.org/avform.htm
      const auto lax = LatLonN(geodesy::Dms::parse("33° 57′N"), geodesy::Dms::parse("118° 24′W"));
      const auto jfk = LatLonN(geodesy::Dms::parse("40° 38′N"), geodesy::Dms::parse("073° 47′W"));
      const auto r = 180*60/pi; // earth radius in nautical miles
      EXPECT_EQ(geodesy::toPrecision(lax.distanceTo(jfk, r), 4), "2144");
      EXPECT_EQ(geodesy::toPrecision(lax.initialBearingTo(jfk), 2), "66");
      EXPECT_EQ(lax.intermediatePointTo(jfk, 100.0/2144.0).toString(geodesy::Dms::DM, 0), "34°37′N, 116°33′W");
      const auto d = LatLonN(geodesy::Dms::parse("34:30N"), geodesy::Dms::parse("116:30W"));
      EXPECT_EQ(geodesy::toPrecision(d.crossTrackDistanceTo(lax, jfk, r), 5), "7.4523");
      EXPECT_EQ(geodesy::toPrecision(d.alongTrackDistanceTo(lax, jfk, r), 5), "99.588");
      EXPECT_EQ(lax.intermediatePointTo(jfk, 0.4).toString(geodesy::Dms::DM, 3), "38°40.167′N, 101°37.570′W");
      const auto reo = LatLonN(geodesy::Dms::parse("42.600N"), geodesy::Dms::parse("117.866W"));
      const auto bke = LatLonN(geodesy::Dms::parse("44.840N"), geodesy::Dms::parse("117.806W"));
      EXPECT_EQ(LatLonN::intersection(reo, 51, bke, 137).toString(geodesy::Dms::D, 3), "43.572°N, 116.189°W");
   }

   TEST_F(latlon_nvector_spherical_unittest, misc)
   {
      EXPECT_TRUE(LatLonN(52.205, 0.119).equals(LatLonN(52.205, 0.119)));
      EXPECT_FALSE(LatLonN(52.206, 0.119).equals(LatLonN(52.205, 0.119)));
      EXPECT_EQ(LatLonN(52.205, 0.119).toGeoJSON(), "{ type: \"Point\", coordinates: [0.119, 52.205] }");
   }

   TEST_F(latlon_nvector_spherical_unittest, example_5_surface_distance)
   {
      const auto a = LatLonN(88, 0);
      const auto b = LatLonN(89, -170);
      const auto dist = a.distanceTo(b); // 332.5 km
      EXPECT_EQ(geodesy::toFixed(dist), "332456");
   }

   TEST_F(latlon_nvector_spherical_unittest,example_6_Interpolated_position)
   {
      const auto a = LatLonN(89, 0);
      const auto b = LatLonN(89, 180);
      const auto i = a.intermediatePointOnChordTo(b, (double)(16-10)/(20-10));
      EXPECT_EQ(i.toString(geodesy::Dms::D, 7), "89.7999805°N, 180.0000000°E");
   }

   TEST_F(latlon_nvector_spherical_unittest, example_7_Mean_position)
   {
      const std::vector points = {
         LatLonN(90, 0),
         LatLonN(60, 10),
         LatLonN(50, -20),
      };
      const auto mean = LatLonN::meanOf(points); // 67.2362°N, 006.9175°W
      EXPECT_EQ(mean.toString(), "67.2362°N, 006.9175°W");
   }

   TEST_F(latlon_nvector_spherical_unittest, example_8_A_and_azimuth_distance_to_B)
   {
      const auto a = LatLonN(80, -90);
      const auto b = a.destinationPoint(1000, 200); // 79.9915°N, 090.0177°W
      EXPECT_EQ(b.toString(), "79.9915°N, 090.0177°W");
   }

   TEST_F(latlon_nvector_spherical_unittest, example_9_Intersection_of_two_paths)
   {
      const auto a1 = LatLonN(10, 20);
      const auto a2 = LatLonN(30, 40);
      const auto b1 = LatLonN(50, 60);
      const auto b2 = LatLonN(70, 80);
      const auto c = LatLonN::intersection(a1, a2, b1, b2); // 40.3186°N, 055.9019°E
      EXPECT_EQ(c.toString(), "40.3186°N, 055.9019°E");
   }

   TEST_F(latlon_nvector_spherical_unittest, example_10_Cross_track_distance)
   {
      const auto a1 = LatLonN(0, 0);
      const auto a2 = LatLonN(10, 0);
      const auto b = LatLonN(1, 0.1);
      const auto c = b.crossTrackDistanceTo(a1, a2); // 11.12 km
      EXPECT_EQ(geodesy::toFixed(c), "11118");
   }
}

