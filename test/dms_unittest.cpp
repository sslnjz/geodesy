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

#include "geodesy/dms.h"

TEST(dms_unittest, _0_)
{
    EXPECT_EQ(0, geodesy::Dms::parse("0.0°"));
    EXPECT_EQ(0, geodesy::Dms::parse("0°"));
    EXPECT_EQ(0, geodesy::Dms::parse("000 00 00 "));
    EXPECT_EQ(0, geodesy::Dms::parse("000°00′00″"));
    EXPECT_EQ(0, geodesy::Dms::parse("000°00′00.0″"));

    EXPECT_EQ(geodesy::Dms::toDms(0), "000.0000°");
    EXPECT_EQ(geodesy::Dms::toDms(0, geodesy::Dms::D, 0), "000°");
    EXPECT_EQ(geodesy::Dms::toDms(0, geodesy::Dms::DMS), "000°00′00″");
    EXPECT_EQ(geodesy::Dms::toDms(0, geodesy::Dms::DMS, 2), "000°00′00.00″");
}

TEST(dms_unittest, parse_variations)
{
   const std::string variations[] = {
         "45.76260",
         "45.76260 ",
         "45.76260°",
         "45°45.756′",
         "45° 45.756′",
         "45 45.756",
         "45°45′45.36″",
         "45º45'45.36\"",
         "45°45’45.36”",
         "45 45 45.36 ",
         "45° 45′ 45.36″",
         "45º 45' 45.36\"",
         "45° 45’ 45.36”",
   };

   for (const auto& v : variations) EXPECT_DOUBLE_EQ(45.76260,  geodesy::Dms::parse(v));
   for (const auto& v : variations) EXPECT_DOUBLE_EQ(-45.76260, geodesy::Dms::parse("-" + v));
   for (const auto& v : variations) EXPECT_DOUBLE_EQ(45.76260,  geodesy::Dms::parse(v + "N"));
   for (const auto& v : variations) EXPECT_DOUBLE_EQ(-45.76260, geodesy::Dms::parse(v + "S"));
   for (const auto& v : variations) EXPECT_DOUBLE_EQ(45.76260,  geodesy::Dms::parse(v + "E"));
   for (const auto& v : variations) EXPECT_DOUBLE_EQ(-45.76260, geodesy::Dms::parse(v + "W"));
   for (const auto& v : variations) EXPECT_DOUBLE_EQ(45.76260,  geodesy::Dms::parse(v + " N"));
   for (const auto& v : variations) EXPECT_DOUBLE_EQ(-45.76260, geodesy::Dms::parse(v + " S"));
   for (const auto& v : variations) EXPECT_DOUBLE_EQ(45.76260,  geodesy::Dms::parse(v + " E"));
   for (const auto& v : variations) EXPECT_DOUBLE_EQ(-45.76260, geodesy::Dms::parse(v + " W"));

   EXPECT_DOUBLE_EQ(45.76260, geodesy::Dms::parse(" 45°45′45.36″ "));
}

TEST(dms_unittest, parse_out_of_range)
{
   EXPECT_DOUBLE_EQ(geodesy::Dms::parse("185"), 185);
   EXPECT_DOUBLE_EQ(geodesy::Dms::parse("365"), 365);
   EXPECT_DOUBLE_EQ(geodesy::Dms::parse("-185"), -185);
   EXPECT_DOUBLE_EQ(geodesy::Dms::parse("-365"), -365);
}

TEST(dms_unittest, output_variations)
{
   geodesy::Dms::set_separator("");
   EXPECT_EQ(geodesy::Dms::toDms(9.1525),                         "009.1525°");
   EXPECT_EQ(geodesy::Dms::toDms(9.1525, geodesy::Dms::D),        "009.1525°");
   EXPECT_EQ(geodesy::Dms::toDms(9.1525, geodesy::Dms::DM),       "009°09.15′");
   EXPECT_EQ(geodesy::Dms::toDms(9.1525, geodesy::Dms::DMS),      "009°09′09″");
   EXPECT_EQ(geodesy::Dms::toDms(9.1525, geodesy::Dms::D, 6),     "009.152500°");
   EXPECT_EQ(geodesy::Dms::toDms(9.1525, geodesy::Dms::DM, 4),    "009°09.1500′");
   EXPECT_EQ(geodesy::Dms::toDms(9.1525, geodesy::Dms::DMS, 2),   "009°09′09.00″");
   EXPECT_EQ(geodesy::Dms::toDms(9.1525, geodesy::Dms::D),        "009.1525°");
   EXPECT_EQ(geodesy::Dms::toDms(9.1525, geodesy::Dms::D, 6),     "009.152500°");
}

TEST(dms_unittest, compass_points)
{
   EXPECT_EQ(geodesy::Dms::compassPoint(1),      "N");
   EXPECT_EQ(geodesy::Dms::compassPoint(0),      "N");
   EXPECT_EQ(geodesy::Dms::compassPoint(-1),     "N");
   EXPECT_EQ(geodesy::Dms::compassPoint(359),    "N");
   EXPECT_EQ(geodesy::Dms::compassPoint(24),     "NNE");
   EXPECT_EQ(geodesy::Dms::compassPoint(24, 1),  "N");
   EXPECT_EQ(geodesy::Dms::compassPoint(24, 2),  "NE");
   EXPECT_EQ(geodesy::Dms::compassPoint(24, 3),  "NNE");
   EXPECT_EQ(geodesy::Dms::compassPoint(226),    "SW");
   EXPECT_EQ(geodesy::Dms::compassPoint(226, 1), "W");
   EXPECT_EQ(geodesy::Dms::compassPoint(226, 2), "SW");
   EXPECT_EQ(geodesy::Dms::compassPoint(226, 3), "SW");
   EXPECT_EQ(geodesy::Dms::compassPoint(237),    "WSW");
   EXPECT_EQ(geodesy::Dms::compassPoint(237, 1), "W");
   EXPECT_EQ(geodesy::Dms::compassPoint(237, 2), "SW");
   EXPECT_EQ(geodesy::Dms::compassPoint(237, 3), "WSW");
   EXPECT_THROW(geodesy::Dms::compassPoint(0, 0), std::range_error);
}

TEST(dms_unittest, misc)
{
   //needs to align separator
   geodesy::Dms::set_separator("");
   EXPECT_EQ(geodesy::Dms::toLat(51.2, geodesy::Dms::DMS),              "51°12′00″N");
   EXPECT_EQ(geodesy::Dms::toLat(51.19999999999999, geodesy::Dms::DM),  "51°12.00′N");
   EXPECT_EQ(geodesy::Dms::toLat(51.19999999999999, geodesy::Dms::DMS), "51°12′00″N");
   EXPECT_EQ(geodesy::Dms::toLon(0.33, geodesy::Dms::DMS),              "000°19′48″E");
   EXPECT_EQ(geodesy::Dms::toDms(51.99999999999999, geodesy::Dms::D),   "052.0000°");
   EXPECT_EQ(geodesy::Dms::toDms(51.99999999999999, geodesy::Dms::DM),  "052°00.00′");
   EXPECT_EQ(geodesy::Dms::toDms(51.99999999999999, geodesy::Dms::DMS), "052°00′00″");
   EXPECT_EQ(geodesy::Dms::toDms(51.19999999999999, geodesy::Dms::D),   "051.2000°");
   EXPECT_EQ(geodesy::Dms::toDms(51.19999999999999, geodesy::Dms::DM),  "051°12.00′");
   EXPECT_EQ(geodesy::Dms::toDms(51.19999999999999, geodesy::Dms::DMS), "051°12′00″");
   EXPECT_EQ(geodesy::Dms::toBearing(1),                                "001.0000°");
   EXPECT_EQ(geodesy::Dms::toLocale("123,456.789"),                     "123,456.789");
   EXPECT_EQ(geodesy::Dms::fromLocale("51°28′40.12″N"),                 "51°28′40.12″N");
   EXPECT_EQ(geodesy::Dms::fromLocale("51°28′40.12″N, 000°00′05.31″W"), "51°28′40.12″N, 000°00′05.31″W");
}

TEST(dms_unittest, parse_failures)
{
   EXPECT_TRUE(std::isnan(geodesy::Dms::parse("0 0 0 0")));
   EXPECT_TRUE(std::isnan(geodesy::Dms::parse("xxx")));
   EXPECT_TRUE(std::isnan(geodesy::Dms::parse("")));
}


TEST(dms_unittest, wrap360)
{
   EXPECT_EQ(geodesy::Dms::wrap360(-450),270);
   EXPECT_EQ(geodesy::Dms::wrap360(-405),315);
   EXPECT_EQ(geodesy::Dms::wrap360(-360),0);
   EXPECT_EQ(geodesy::Dms::wrap360(-315),45);
   EXPECT_EQ(geodesy::Dms::wrap360(-270),90);
   EXPECT_EQ(geodesy::Dms::wrap360(-225),135);
   EXPECT_EQ(geodesy::Dms::wrap360(-180),180);
   EXPECT_EQ(geodesy::Dms::wrap360(-135),225);
   EXPECT_EQ(geodesy::Dms::wrap360(-90),270);
   EXPECT_EQ(geodesy::Dms::wrap360(-45),315);
   EXPECT_EQ(geodesy::Dms::wrap360(0),0);
   EXPECT_EQ(geodesy::Dms::wrap360(45),45);
   EXPECT_EQ(geodesy::Dms::wrap360(90),90);
   EXPECT_EQ(geodesy::Dms::wrap360(135),135);
   EXPECT_EQ(geodesy::Dms::wrap360(180),180);
   EXPECT_EQ(geodesy::Dms::wrap360(225),225);
   EXPECT_EQ(geodesy::Dms::wrap360(270),270);
   EXPECT_EQ(geodesy::Dms::wrap360(315),315);
   EXPECT_EQ(geodesy::Dms::wrap360(360),0);
   EXPECT_EQ(geodesy::Dms::wrap360(405),45);
   EXPECT_EQ(geodesy::Dms::wrap360(450),90);
}

TEST(dms_unittest, wrap180)
{
   EXPECT_EQ(geodesy::Dms::wrap180(-450),  -90);
   EXPECT_EQ(geodesy::Dms::wrap180(-405),  -45);
   EXPECT_EQ(geodesy::Dms::wrap180(-360),    0);
   EXPECT_EQ(geodesy::Dms::wrap180(-315),   45);
   EXPECT_EQ(geodesy::Dms::wrap180(-270),   90);
   EXPECT_EQ(geodesy::Dms::wrap180(-225),  135);
   EXPECT_EQ(geodesy::Dms::wrap180(-180), -180);
   EXPECT_EQ(geodesy::Dms::wrap180(-135), -135);
   EXPECT_EQ(geodesy::Dms::wrap180(-90),   -90);
   EXPECT_EQ(geodesy::Dms::wrap180(-45),   -45);
   EXPECT_EQ(geodesy::Dms::wrap180(0),       0);
   EXPECT_EQ(geodesy::Dms::wrap180(45),     45);
   EXPECT_EQ(geodesy::Dms::wrap180(90),     90);
   EXPECT_EQ(geodesy::Dms::wrap180(135),   135);
   EXPECT_EQ(geodesy::Dms::wrap180(180),   180);
   EXPECT_EQ(geodesy::Dms::wrap180(225),  -135);
   EXPECT_EQ(geodesy::Dms::wrap180(270),   -90);
   EXPECT_EQ(geodesy::Dms::wrap180(315),   -45);
   EXPECT_EQ(geodesy::Dms::wrap180(360),     0);
   EXPECT_EQ(geodesy::Dms::wrap180(405),    45);
   EXPECT_EQ(geodesy::Dms::wrap180(450),    90);
}


TEST(dms_unittest, wrap90)
{
   EXPECT_EQ(geodesy::Dms::wrap90(-450), -90);
   EXPECT_EQ(geodesy::Dms::wrap90(-405), -45);
   EXPECT_EQ(geodesy::Dms::wrap90(-360),   0);
   EXPECT_EQ(geodesy::Dms::wrap90(-315),  45);
   EXPECT_EQ(geodesy::Dms::wrap90(-270),  90);
   EXPECT_EQ(geodesy::Dms::wrap90(-225),  45);
   EXPECT_EQ(geodesy::Dms::wrap90(-180),   0);
   EXPECT_EQ(geodesy::Dms::wrap90(-135), -45);
   EXPECT_EQ(geodesy::Dms::wrap90(-90),  -90);
   EXPECT_EQ(geodesy::Dms::wrap90(-45),  -45);
   EXPECT_EQ(geodesy::Dms::wrap90(0),      0);
   EXPECT_EQ(geodesy::Dms::wrap90(45),    45);
   EXPECT_EQ(geodesy::Dms::wrap90(90),    90);
   EXPECT_EQ(geodesy::Dms::wrap90(135),   45);
   EXPECT_EQ(geodesy::Dms::wrap90(180),    0);
   EXPECT_EQ(geodesy::Dms::wrap90(225),  -45);
   EXPECT_EQ(geodesy::Dms::wrap90(270),  -90);
   EXPECT_EQ(geodesy::Dms::wrap90(315),  -45);
   EXPECT_EQ(geodesy::Dms::wrap90(360),    0);
   EXPECT_EQ(geodesy::Dms::wrap90(405),   45);
   EXPECT_EQ(geodesy::Dms::wrap90(450),   90);
}
