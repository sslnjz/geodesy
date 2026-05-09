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

#include "geodesy/latlon_osgridref.h"
#include "geodesy/osgridref.h"

#include <limits>
#include <stdexcept>

#include <gtest/gtest.h>

namespace
{

using geodesy::LatLonEllipsoidalDatum;
using geodesy::LatLonOsGridRef;
using geodesy::OsGridRef;

TEST(osgridref_unittest, formats_reference_example_at_multiple_precisions)
{
    const OsGridRef grid(651409.0, 313177.0);

    EXPECT_EQ(grid.toString(), "TG 51409 13177");
    EXPECT_EQ(grid.toString(8), "TG 5140 1317");
    EXPECT_EQ(grid.toString(6), "TG 514 131");
    EXPECT_EQ(grid.toString(0), "651409,313177");
}

TEST(osgridref_unittest, parse_accepts_standard_compact_and_numeric_references)
{
    EXPECT_EQ(OsGridRef::parse("TG 51409 13177").toString(), "TG 51409 13177");
    EXPECT_EQ(OsGridRef::parse("tg5140913177").toString(), "TG 51409 13177");
    EXPECT_EQ(OsGridRef::parse("651409,313177").toString(), "TG 51409 13177");
    EXPECT_EQ(OsGridRef::parse("SU 387 148").toString(), "SU 38700 14800");
}

TEST(osgridref_unittest, os_grid_to_latlon_matches_reference_in_wgs84_and_osgb36)
{
    const OsGridRef grid(651409.903, 313177.270);

    const LatLonOsGridRef wgs84 = grid.toLatLon();
    EXPECT_NEAR(wgs84.lat(), 52.65797861, 2e-7);
    EXPECT_NEAR(wgs84.lon(), 1.71605194, 2e-7);
    EXPECT_EQ(wgs84.datum(), LatLonEllipsoidalDatum::datums().WGS84);

    const LatLonOsGridRef osgb36 = grid.toLatLon(LatLonEllipsoidalDatum::datums().OSGB36);
    EXPECT_NEAR(osgb36.lat(), 52.65757028, 2e-7);
    EXPECT_NEAR(osgb36.lon(), 1.71792167, 2e-7);
    EXPECT_EQ(osgb36.datum(), LatLonEllipsoidalDatum::datums().OSGB36);
}

TEST(osgridref_unittest, latlon_to_os_grid_matches_reference_example)
{
    const LatLonOsGridRef point(52.65798, 1.71605);

    EXPECT_EQ(point.toOsGrid().toString(), "TG 51409 13177");
}

TEST(osgridref_unittest, osgb36_latlon_to_os_grid_does_not_apply_second_datum_conversion)
{
    const LatLonOsGridRef point(52.65757028, 1.71792167, 0.0, LatLonEllipsoidalDatum::datums().OSGB36);

    EXPECT_EQ(point.toOsGrid().toString(), "TG 51409 13177");
}

TEST(osgridref_unittest, rejects_invalid_constructor_parse_and_format_inputs)
{
    EXPECT_THROW(OsGridRef(-1.0, 313177.0), std::invalid_argument);
    EXPECT_THROW(OsGridRef(651409.0, 1300000.001), std::invalid_argument);
    EXPECT_THROW(OsGridRef(std::numeric_limits<double>::infinity(), 313177.0), std::invalid_argument);

    EXPECT_THROW((void)OsGridRef::parse(""), std::invalid_argument);
    EXPECT_THROW((void)OsGridRef::parse("TG 51409 1317"), std::invalid_argument);
    EXPECT_THROW((void)OsGridRef::parse("TI 51409 13177"), std::invalid_argument);
    EXPECT_THROW((void)OsGridRef::parse("XX 51409 13177"), std::invalid_argument);
    EXPECT_THROW((void)OsGridRef::parse("700001,313177"), std::invalid_argument);

    EXPECT_THROW((void)OsGridRef(651409.0, 313177.0).toString(1), std::invalid_argument);
    EXPECT_THROW((void)OsGridRef(651409.0, 313177.0).toString(18), std::invalid_argument);
}

}
