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

#include <cmath>

#include "geodesy/latlon.h"
#include "geodesy/dms.h"

TEST(latlon_unittest, examples)
{
   geodesy::Dms::set_separator(""); // tests are easier without any DMS separator

   EXPECT_EQ(geodesy::LatLon(52.205, 0.119).toString(),                          "52.2050°N, 000.1190°E");
   EXPECT_EQ(geodesy::LatLon::parse(52.205, 0.119).toString(),                   "52.2050°N, 000.1190°E");
   EXPECT_EQ(geodesy::LatLon::parse("52.205", "0.119").toString(),               "52.2050°N, 000.1190°E");
   EXPECT_EQ(geodesy::LatLon::parse("52.205, 0.119").toString(),                 "52.2050°N, 000.1190°E");
   EXPECT_EQ(geodesy::LatLon::parse("52°12′18.0″N", "000°07′08.4″E").toString(), "52.2050°N, 000.1190°E");
   EXPECT_EQ(geodesy::LatLon::parse("52°12′18.0″N, 000°07′08.4″E").toString(),   "52.2050°N, 000.1190°E");

   EXPECT_EQ(geodesy::LatLon(52.205, 0.119), geodesy::LatLon(52.205, 0.119));

   const auto greenwich = geodesy::LatLon(51.47788, -0.00147);
   EXPECT_EQ(greenwich.toString(),                       "51.4779°N, 000.0015°W");
   EXPECT_EQ(greenwich.toString(geodesy::Dms::DMS, 2),   "51°28′40.37″N, 000°00′05.29″W");
   EXPECT_EQ(greenwich.toString(geodesy::Dms::N),        "51.4779,-0.0015");
}

TEST(latlon_unittest, constructor_with_strings)
{
   EXPECT_EQ(geodesy::LatLon("52.205", "0.119"), geodesy::LatLon(52.205, 0.119));
}

TEST(latlon_unittest, constructor_fail)
{
   EXPECT_THROW(geodesy::LatLon("x", "x"), std::invalid_argument);
   EXPECT_THROW(geodesy::LatLon("x", "x"), std::invalid_argument);
};

TEST(latlon_unittest, parse_fail)
{
   EXPECT_THROW(geodesy::LatLon::parse("cam", "bridge"), std::invalid_argument);
   EXPECT_THROW(geodesy::LatLon::parse("cambridge"),     std::invalid_argument);
   EXPECT_THROW(geodesy::LatLon::parse(NAN, NAN), std::invalid_argument);
}

TEST(latlon_unittest, getters_setters)
{
   auto camb = geodesy::LatLon();
   //camb.lat = camb.latitude = "52° 12′ 18″ N";
   //camb.lon = camb.lng = camb.longitude = "000° 07′ 08″ E";
   camb.setLat("52° 12′ 18″ N");
   camb.setLon("000° 07′ 08″ E");
   EXPECT_NEAR(camb.lat(),       52.205, 0.001);
   EXPECT_NEAR(camb.latitude(),  52.205, 0.001);
   EXPECT_NEAR(camb.lon(),       0.119,  0.001);
   EXPECT_NEAR(camb.lng(),       0.119,  0.001);
   EXPECT_NEAR(camb.longitude(), 0.119,  0.001);
                                             
   camb.setLatitude(52.205);                 
   camb.setLongitude(0.119);                 
   EXPECT_NEAR(camb.lat(),       52.205, 0.001);
   EXPECT_NEAR(camb.latitude(),  52.205, 0.001);
   EXPECT_NEAR(camb.lon(),       0.119,  0.001);
   EXPECT_NEAR(camb.lng(),       0.119,  0.001);
   EXPECT_NEAR(camb.longitude(), 0.119,  0.001);
}

TEST(latlon_unittest, setters_fail)
{
   auto camb = geodesy::LatLon(0, 0);
   EXPECT_THROW(camb.setLat("xxx"),       std::invalid_argument);
   EXPECT_THROW(camb.setLatitude("xxx"),  std::invalid_argument);
   EXPECT_THROW(camb.setLon("xxx"),       std::invalid_argument);
   EXPECT_THROW(camb.setLng("xxx"),       std::invalid_argument);
   EXPECT_THROW(camb.setLongitude("xxx"), std::invalid_argument);
};

TEST(latlon_unittest, toString)
{
   const auto btTower = geodesy::LatLon(51.521470, -0.138833);
   EXPECT_EQ(btTower.toString(),                      "51.5215°N, 000.1388°W");
   EXPECT_EQ(btTower.toString(geodesy::Dms::D, 6),    "51.521470°N, 000.138833°W");
   EXPECT_EQ(btTower.toString(geodesy::Dms::DM, 4),   "51°31.2882′N, 000°08.3300′W");
   EXPECT_EQ(btTower.toString(geodesy::Dms::DMS, 2),  "51°31′17.29″N, 000°08′19.80″W");
   EXPECT_EQ(btTower.toString(geodesy::Dms::N),       "51.5215,-0.1388");
   EXPECT_EQ(btTower.toString(geodesy::Dms::N, 6),    "51.521470,-0.138833");
}

TEST(latlon_unittest, misc)
{
   EXPECT_TRUE(geodesy::LatLon(52.205, 0.119) == geodesy::LatLon(52.205, 0.119));
   EXPECT_FALSE(geodesy::LatLon(52.206, 0.119) == geodesy::LatLon(52.205, 0.119));
   EXPECT_EQ(geodesy::LatLon(52.205, 0.119).toGeoJSON(), "{ type: \"Point\", coordinates : [0.119, 52.205] }");
};