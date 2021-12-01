// /*****************************************************************************
// *                                                                            *
// *  COPYRIGHT 2021 ACCEL FLIGHT SIMULATION. ALL RIGHTS RESERVED               *
// *                                                                            *
// *  This file is part of ACCEL Platform Software                              *
// *                                                                            *
// *  All information and content in this document is the Property              *
// *  to ACCEL (Tianjin) Flight Simulation Co., Ltd. and shall not              *
// *  be disclosed, disseminated, copied, or used except for purposes           *
// *  expressly authorized in writing by ACCEL. Any unauthorized use            *
// *  of the content will be considered an infringement of                      *
// *  ACCEL's intellectual property rights.                                     *
// *                                                                            *
// *  @file     latlon_spherical_unittest.cpp                                   *
// *  @brief                                                                    *
// *                                                                            *
// *                                                                            *
// *  @author   Song, Binbin                                                    *
// *  @email    binbin.song@accelflightsimulation.com                           *
// *  @version  0.1.0(version)                                                  *
// *  @date     2021/11/30                                                      *
// *                                                                            *
// *----------------------------------------------------------------------------*
// *  Remark         : Description                                              *
// *----------------------------------------------------------------------------*
// *  Change History :                                                          *
// *  <Date>     | <Version> | <Author>       | <Description>                   *
// *----------------------------------------------------------------------------*
// *  2021/11/30  | 0.1.0   | Song, Binbin   | Create file                      *
// *----------------------------------------------------------------------------*
// *                                                                            *
// ******************************************************************************/
#include <gtest/gtest.h>

#include <sstream>
#include <cmath>

#include "geodesy/latlon_spherical.h"

std::string toFixed(double value, int p = 0)
{
   std::stringstream ss;
   ss << std::fixed << std::setprecision(p) << value;
   return ss.str();
}

std::string toPrecision(double value, int p = 0)
{
   std::stringstream ss;
   ss << std::setprecision(p) << value;
   return ss.str();
}


std::string toExponential(double value, int p = 0)
{
   std::stringstream ss;
   ss << std::setprecision(p) << std::scientific  << value;
   return ss.str();
}


TEST(latlon_spherical_unittest, constructor)
{
   EXPECT_EQ(toFixed(geodesy::LatLonSpherical(52.205, 0.119).distanceTo(geodesy::LatLonSpherical(48.857, 2.351))), "404279");
   EXPECT_EQ(toFixed(geodesy::LatLonSpherical(52.205, 0.119).distanceTo(geodesy::LatLonSpherical(48.857, 2.351), 3959), 1), "251.2");
   EXPECT_EQ(toFixed(geodesy::LatLonSpherical(52.205, 0.119).initialBearingTo(geodesy::LatLonSpherical(48.857, 2.351)), 1), "156.2");
   EXPECT_EQ(toFixed(geodesy::LatLonSpherical(52.205, 0.119).finalBearingTo(geodesy::LatLonSpherical(48.857, 2.351)), 1), "157.9");
   EXPECT_EQ(toFixed(geodesy::LatLonSpherical(53.2611, -0.7972).crossTrackDistanceTo(geodesy::LatLonSpherical(53.3206, -1.7297), geodesy::LatLonSpherical(53.1887, 0.1334)), 1), "-307.5");
   EXPECT_EQ(toFixed(geodesy::LatLonSpherical(53.2611, -0.7972).alongTrackDistanceTo(geodesy::LatLonSpherical(53.3206, -1.7297), geodesy::LatLonSpherical(53.1887, 0.1334))), "62331");
   EXPECT_EQ(toFixed(geodesy::LatLonSpherical(51.127, 1.338).rhumbDistanceTo(geodesy::LatLonSpherical(50.964, 1.853))), "40308");
   EXPECT_EQ(toFixed(geodesy::LatLonSpherical(51.127, 1.338).rhumbBearingTo(geodesy::LatLonSpherical(50.964, 1.853)), 1), "116.7");

   EXPECT_EQ(geodesy::LatLonSpherical(51.127, 1.338).rhumbDestinationPoint(40300, 116.7).toString(), "50.9642°N, 001.8530°E");
   EXPECT_EQ(geodesy::LatLonSpherical(52.205, 0.119).midpointTo(geodesy::LatLonSpherical(48.857, 2.351)).toString(), "50.5363°N, 001.2746°E");
   EXPECT_EQ(geodesy::LatLonSpherical(52.205, 0.119).intermediatePointTo(geodesy::LatLonSpherical(48.857, 2.351), 0.25).toString(), "51.3721°N, 000.7073°E");
   EXPECT_EQ(geodesy::LatLonSpherical(51.47788, -0.00147).destinationPoint(7794, 300.7).toString(), "51.5136°N, 000.0983°W");
   EXPECT_EQ(geodesy::LatLonSpherical::intersection(geodesy::LatLonSpherical(51.8853, 0.2545), 108.547, geodesy::LatLonSpherical(49.0034, 2.5735), 32.435).toString(), "50.9078°N, 004.5084°E");
   EXPECT_EQ(geodesy::LatLonSpherical(51.127, 1.338).rhumbMidpointTo(geodesy::LatLonSpherical(50.964, 1.853)).toString(), "51.0455°N, 001.5957°E");

   std::vector vec = { geodesy::LatLonSpherical(0, 0), geodesy::LatLonSpherical(1, 0), geodesy::LatLonSpherical(0, 1) };
   EXPECT_EQ(toExponential(geodesy::LatLonSpherical::areaOf(vec), 2), "6.18e+09");
}

TEST(latlon_spherical_unittest, alternate_point_formats)
{
   const auto cambg = geodesy::LatLonSpherical(52.205, 0.119);
   EXPECT_EQ(toExponential(cambg.distanceTo(geodesy::LatLonSpherical(48.857, 2.351)), 3), "4.043e+05");
   EXPECT_EQ(toExponential(cambg.distanceTo(geodesy::LatLonSpherical("48.857", "2.351")), 3), "4.043e+05");
   EXPECT_EQ(toExponential(cambg.distanceTo(geodesy::LatLonSpherical("48°51′25.2″N", "002°21′03.6″E")), 3), "4.043e+05");
}

TEST(latlon_spherical_unittest, dist_brng_dest)
{
   const auto cambg = geodesy::LatLonSpherical(52.205, 0.119), paris = geodesy::LatLonSpherical(48.857, 2.351);
   EXPECT_EQ(toExponential(cambg.distanceTo(paris), 3), "4.043e+05");
   EXPECT_EQ(toFixed(cambg.distanceTo(paris, 3959), 1), "251.2");
   EXPECT_EQ(toFixed(cambg.initialBearingTo(paris), 1), "156.2");
   EXPECT_EQ(toFixed(cambg.finalBearingTo(paris), 1), "157.9");
   EXPECT_TRUE(std::isnan(cambg.initialBearingTo(cambg)));
   EXPECT_TRUE(std::isnan(cambg.finalBearingTo(cambg)));
   EXPECT_EQ(toFixed(paris.initialBearingTo(cambg), 1), "337.9");
   EXPECT_EQ(cambg.midpointTo(paris).toString(), "50.5363°N, 001.2746°E");
   EXPECT_EQ(cambg.intermediatePointTo(paris, 0.25).toString(), "51.3721°N, 000.7073°E");
   EXPECT_EQ(cambg.intermediatePointTo(cambg, 0.25).toString(), "52.2050°N, 000.1190°E");
   const auto greenwich = geodesy::LatLonSpherical(51.47788, -0.00147);
   const double dist = 7794, brng = 300.7;
   EXPECT_EQ(greenwich.destinationPoint(dist, brng).toString(), "51.5136°N, 000.0983°W");
   EXPECT_EQ(greenwich.destinationPoint(dist, brng, 6371e3).toString(), "51.5136°N, 000.0983°W");
}

TEST(latlon_spherical_unittest, dist_brng_dest_fails)
{
   const auto cambg = geodesy::LatLonSpherical(52.205, 0.119), paris = geodesy::LatLonSpherical(48.857, 2.351);
   EXPECT_THROW(cambg.distanceTo({ "xxx", "xxx" }), std::invalid_argument);
};

TEST(latlon_spherical_unittest, intersection)
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
   EXPECT_EQ(geodesy::LatLonSpherical::intersection(geodesy::LatLonSpherical(-77.6966041375563, 18.2812500000000), 179.99999999999995, geodesy::LatLonSpherical(89, 180), 180).toString(), "90.0000°S, 163.9902°W");
}

TEST(latlon_spherical_unittest, cross_track_along_track)
{
   EXPECT_EQ(toExponential(geodesy::LatLonSpherical(10, 1).crossTrackDistanceTo(geodesy::LatLonSpherical(0, 0), geodesy::LatLonSpherical(0, 2)), 3), "-1.112e+06");
//
//   const bradwell = new LatLon(53.3206, -1.7297), dunham = new LatLon(53.2611, -0.7972), partney = new LatLon(53.1887, 0.1334);
//   test('cross-track', () = > dunham.crossTrackDistanceTo(bradwell, partney).toPrecision(4).should.equal('-307.5'));
//   test('cross-track (fail)', () = > should.Throw(function() { new LatLon(10, 1).crossTrackDistanceTo(null, new LatLon(0, 2)); }, TypeError, 'invalid (null) point'));
//   test('cross-track (fail)', () = > should.Throw(function() { new LatLon(10, 1).crossTrackDistanceTo(new LatLon(0, 0), null); }, TypeError, 'invalid (null) point'));
//   test('along-track', () = > dunham.alongTrackDistanceTo(bradwell, partney).toPrecision(4).should.equal('6.233e+4'));
//
//   test('cross-track NE', () = > new LatLon(1, 1).crossTrackDistanceTo(new LatLon(0, 0), new LatLon(0, 2)).toPrecision(4).should.equal('-1.112e+5'));
//   test('cross-track SE', () = > new LatLon(-1, 1).crossTrackDistanceTo(new LatLon(0, 0), new LatLon(0, 2)).toPrecision(4).should.equal('1.112e+5'));
//   test('cross-track SW?', () = > new LatLon(-1, -1).crossTrackDistanceTo(new LatLon(0, 0), new LatLon(0, 2)).toPrecision(4).should.equal('1.112e+5'));
//   test('cross-track NW?', () = > new LatLon(1, -1).crossTrackDistanceTo(new LatLon(0, 0), new LatLon(0, 2)).toPrecision(4).should.equal('-1.112e+5')); // eslint-disable-line space-in-parens
//
//   test('along-track NE', () = > new LatLon(1, 1).alongTrackDistanceTo(new LatLon(0, 0), new LatLon(0, 2)).toPrecision(4).should.equal('1.112e+5')); // eslint-disable-line space-in-parens
//   test('along-track SE', () = > new LatLon(-1, 1).alongTrackDistanceTo(new LatLon(0, 0), new LatLon(0, 2)).toPrecision(4).should.equal('1.112e+5'));
//   test('along-track SW', () = > new LatLon(-1, -1).alongTrackDistanceTo(new LatLon(0, 0), new LatLon(0, 2)).toPrecision(4).should.equal('-1.112e+5'));
//   test('along-track NW', () = > new LatLon(1, -1).alongTrackDistanceTo(new LatLon(0, 0), new LatLon(0, 2)).toPrecision(4).should.equal('-1.112e+5')); // eslint-disable-line space-in-parens
//
//   test('cross-track coinc', () = > new LatLon(10, 0).crossTrackDistanceTo(new LatLon(10, 0), new LatLon(0, 2)).should.equal(0));
//   test('along-track coinc', () = > new LatLon(10, 0).alongTrackDistanceTo(new LatLon(10, 0), new LatLon(0, 2)).should.equal(0));
//   test('cross-track (fail)', () = > should.Throw(function() { new LatLon(0, 0).crossTrackDistanceTo(null, 0); }, TypeError, 'invalid (null) point'));
//   test('cross-track (fail)', () = > should.Throw(function() { new LatLon(0, 0).crossTrackDistanceTo(new LatLon(0, 0), 'x'); }, TypeError, 'invalid point ‘x’'));
//   test('along-track (fail)', () = > should.Throw(function() { new LatLon(0, 0).alongTrackDistanceTo(null, 0); }, TypeError, 'invalid (null) point'));
//   test('along-track (fail)', () = > should.Throw(function() { new LatLon(0, 0).alongTrackDistanceTo(new LatLon(0, 0), 'x'); }, TypeError, 'invalid point ‘x’'));
}
//
//describe('misc', function() {
//   test('Clairaut 0°', () = > new LatLon(0, 0).maxLatitude(0).should.equal(90));
//   test('Clairaut 1°', () = > new LatLon(0, 0).maxLatitude(1).should.equal(89));
//   test('Clairaut 90°', () = > new LatLon(0, 0).maxLatitude(90).should.equal(0));
//
//   const parallels = LatLon.crossingParallels(new LatLon(0, 0), new LatLon(60, 30), 30);
//   test('parallels 1', () = > new LatLon(30, parallels.lon1).toString('dms').should.equal('30°00′00″N, 009°35′39″E'));
//   test('parallels 2', () = > new LatLon(30, parallels.lon2).toString('dms').should.equal('30°00′00″N, 170°24′21″E'));
//   test('parallels -', () = > should.equal(LatLon.crossingParallels(new LatLon(0, 0), new LatLon(10, 60), 60), null));
//   test('parallels coinc', () = > should.equal(LatLon.crossingParallels(new LatLon(0, 0), new LatLon(0, 0), 0), null));
//});
//
//describe('area (polygon-based)', function() {
//   const polyTriangle = [new LatLon(1, 1), new LatLon(2, 1), new LatLon(1, 2)];
//   const polySquareCw = [new LatLon(1, 1), new LatLon(2, 1), new LatLon(2, 2), new LatLon(1, 2)];
//   const polySquareCcw = [new LatLon(1, 1), new LatLon(1, 2), new LatLon(2, 2), new LatLon(2, 1)];
//   const polyOctant = [new LatLon(0, ε), new LatLon(90, 0), new LatLon(0, 90 - ε)];
//   const polyOctantS = [new LatLon(-ε, ε), new LatLon(90, 0), new LatLon(-ε, 90 - ε)];
//   const polyQuadrant = [new LatLon(ε, ε), new LatLon(90, ε), new LatLon(ε, 180 - ε), new LatLon(ε, 90)];
//   const polyHemiE = [new LatLon(ε, ε), new LatLon(90 - ε, 0), new LatLon(90 - ε, 180), new LatLon(ε, 180), new LatLon(-ε, 180), new LatLon(-90 + ε, 180), new LatLon(-90 + ε, 0), new LatLon(-ε, ε)];
//   const polyPole = [new LatLon(89, 0), new LatLon(89, 120), new LatLon(89, -120)];
//   const polyConcave = [new LatLon(1, 1), new LatLon(5, 1), new LatLon(5, 3), new LatLon(1, 3), new LatLon(3, 2)];
//
//   test('triangle area', () = > LatLon.areaOf(polyTriangle).toFixed(0).should.equal('6181527888'));
//   test('triangle area radius', () = > LatLon.areaOf(polyTriangle, 6371e3).toFixed(0).should.equal('6181527888'));
//   test('triangle area closed', () = > LatLon.areaOf(polyTriangle.concat(polyTriangle[0])).toFixed(0).should.equal('6181527888'));
//   test('square cw area', () = > LatLon.areaOf(polySquareCw).toFixed(0).should.equal('12360230987'));
//   test('square ccw area', () = > LatLon.areaOf(polySquareCcw).toFixed(0).should.equal('12360230987'));
//   test('octant area', () = > LatLon.areaOf(polyOctant).toFixed(1).should.equal((π * R * R / 2).toFixed(1)));
//   test('super-octant area', () = > LatLon.areaOf(polyOctantS).toFixed(1).should.equal((π * R * R / 2).toFixed(1)));
//   test('quadrant area', () = > LatLon.areaOf(polyQuadrant).should.equal(π * R * R));
//   test('hemisphere area', () = > LatLon.areaOf(polyHemiE).toFixed(1).should.equal((2 * π * R * R).toFixed(1)));
//   test('pole area', () = > LatLon.areaOf(polyPole).toFixed(0).should.equal('16063139192'));
//   test('concave area', () = > LatLon.areaOf(polyConcave).toFixed(0).should.equal('74042699236'));
//});
//
//describe('Ed Williams', function() { // www.edwilliams.org/avform.htm
//   const lax = new LatLon(Dms.parse('33° 57′N'), Dms.parse('118° 24′W'));
//   const jfk = new LatLon(Dms.parse('40° 38′N'), Dms.parse('073° 47′W'));
//   const r = 180 * 60 / π; // earth radius in nautical miles
//   test('distance nm', () = > lax.distanceTo(jfk, r).toPrecision(4).should.equal('2144'));
//   test('bearing', () = > lax.initialBearingTo(jfk).toPrecision(2).should.equal('66'));
//   test('intermediate', () = > lax.intermediatePointTo(jfk, 100 / 2144).toString('dm', 0).should.equal('34°37′N, 116°33′W'));
//   const d = new LatLon(Dms.parse('34:30N'), Dms.parse('116:30W'));
//   test('cross-track', () = > d.crossTrackDistanceTo(lax, jfk, r).toPrecision(5).should.equal('7.4523'));
//   test('along-track', () = > d.alongTrackDistanceTo(lax, jfk, r).toPrecision(5).should.equal('99.588'));
//   test('intermediate', () = > lax.intermediatePointTo(jfk, 0.4).toString('dm', 3).should.equal('38°40.167′N, 101°37.570′W'));
//   const reo = new LatLon(Dms.parse('42.600N'), Dms.parse('117.866W'));
//   const bke = new LatLon(Dms.parse('44.840N'), Dms.parse('117.806W'));
//   test('intersection', () = > LatLon.intersection(reo, 51, bke, 137).toString('d', 3).should.equal('43.572°N, 116.189°W'));
//});
//
//describe('rhumb lines', function() {
//   const dov = new LatLon(51.127, 1.338), cal = new LatLon(50.964, 1.853);
//   test('distance', () = > dov.rhumbDistanceTo(cal).toPrecision(4).should.equal('4.031e+4'));
//   test('distance r', () = > dov.rhumbDistanceTo(cal, 6371e3).toPrecision(4).should.equal('4.031e+4'));
//   test('distance dateline E-W', () = > new LatLon(1, -179).rhumbDistanceTo(new LatLon(1, 179)).toFixed(6).should.equal(new LatLon(1, 1).rhumbDistanceTo(new LatLon(1, -1)).toFixed(6)));
//   test('distance err', () = > should.Throw(function() { dov.rhumbDistanceTo(false); }, TypeError, 'invalid point ‘false’'));
//   test('bearing', () = > dov.rhumbBearingTo(cal).toFixed(1).should.equal('116.7'));
//   test('bearing dateline', () = > new LatLon(1, -179).rhumbBearingTo(new LatLon(1, 179)).should.equal(270));
//   test('bearing dateline', () = > new LatLon(1, 179).rhumbBearingTo(new LatLon(1, -179)).should.equal(90));
//   test('bearing coincident', () = > dov.rhumbBearingTo(dov).should.be.NaN);
//   test('bearing err', () = > should.Throw(function() { dov.rhumbBearingTo(false); }, TypeError, 'invalid point ‘false’'));
//   test('dest’n', () = > dov.rhumbDestinationPoint(40310, 116.7).toString().should.equal('50.9641°N, 001.8531°E'));
//   test('dest’n', () = > dov.rhumbDestinationPoint(40310, 116.7, 6371e3).toString().should.equal('50.9641°N, 001.8531°E'));
//   test('dest’n', () = > new LatLon(1, 1).rhumbDestinationPoint(111178, 90).toString().should.equal('01.0000°N, 002.0000°E'));
//   test('dest’n dateline', () = > new LatLon(1, 179).rhumbDestinationPoint(222356, 90).toString().should.equal('01.0000°N, 179.0000°W'));
//   test('dest’n dateline', () = > new LatLon(1, -179).rhumbDestinationPoint(222356, 270).toString().should.equal('01.0000°N, 179.0000°E'));
//   test('midpoint', () = > dov.rhumbMidpointTo(cal).toString().should.equal('51.0455°N, 001.5957°E'));
//   test('midpoint dateline', () = > new LatLon(1, -179).rhumbMidpointTo(new LatLon(1, 178)).toString().should.equal('01.0000°N, 179.5000°E'));
//   test('midpoint err', () = > should.Throw(function() { dov.rhumbMidpointTo(false); }, TypeError, 'invalid point ‘false’'));
//});