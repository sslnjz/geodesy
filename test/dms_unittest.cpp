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

#include "geodesy/dms.h"

TEST(dms_unittest, parse_0)
{
    EXPECT_EQ(0, geodesy::Dms::parse("0.0°"));
    EXPECT_EQ(0, geodesy::Dms::parse("0°"));
    EXPECT_EQ(0, geodesy::Dms::parse("000 00 00 "));
    EXPECT_EQ(0, geodesy::Dms::parse("000°00′00″"));
    EXPECT_EQ(0, geodesy::Dms::parse("000°00′00.0″"));
}

TEST(dms_unittest, toDms_0)
{
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
         "45º45\'45.36\"",
         "45°45’45.36”",
         "45 45 45.36 ",
         "45° 45′ 45.36″",
         "45º 45\' 45.36\"",
         "45° 45’ 45.36”",
   };

   for (auto var: variations) EXPECT_DOUBLE_EQ(45.76260, geodesy::Dms::parse(var));
   for (auto var: variations) EXPECT_DOUBLE_EQ(-45.76260, geodesy::Dms::parse("-" + var));
   for (auto var: variations) EXPECT_DOUBLE_EQ(45.76260, geodesy::Dms::parse(var + "N"));
   for (auto var: variations) EXPECT_DOUBLE_EQ(-45.76260, geodesy::Dms::parse(var + "S"));
   for (auto var: variations) EXPECT_DOUBLE_EQ(45.76260, geodesy::Dms::parse(var + "E"));
   for (auto var: variations) EXPECT_DOUBLE_EQ(-45.76260, geodesy::Dms::parse(var + "W"));
   for (auto var: variations) EXPECT_DOUBLE_EQ(45.76260, geodesy::Dms::parse(var + " N"));
   for (auto var: variations) EXPECT_DOUBLE_EQ(-45.76260, geodesy::Dms::parse(var + " S"));
   for (auto var: variations) EXPECT_DOUBLE_EQ(45.76260, geodesy::Dms::parse(var + " E"));
   for (auto var: variations) EXPECT_DOUBLE_EQ(-45.76260, geodesy::Dms::parse(var + " W"));
}
