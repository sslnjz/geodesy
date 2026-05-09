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

#include "geodesy/latlon_utm_mgrs.h"
#include "geodesy/mgrs.h"
#include "geodesy/utm_mgrs.h"

#include <cmath>
#include <limits>
#include <stdexcept>

#include <gtest/gtest.h>

namespace
{

using geodesy::LatLonUtmMgrs;
using geodesy::Mgrs;
using geodesy::Utm;
using geodesy::UtmMgrs;

TEST(mgrs_unittest, formats_reference_example_at_multiple_precisions)
{
    const Mgrs mgrs(31, 'U', 'D', 'Q', 48251.795, 11932.678);

    EXPECT_EQ(mgrs.toString(), "31U DQ 48251 11932");
    EXPECT_EQ(mgrs.toString(4), "31U DQ 48 11");
    EXPECT_EQ(mgrs.toString(0), "31U DQ");
}

TEST(mgrs_unittest, parse_accepts_spaced_and_military_style_references)
{
    EXPECT_EQ(Mgrs::parse("31U DQ 48251 11932").toString(), "31U DQ 48251 11932");
    EXPECT_EQ(Mgrs::parse("31UDQ4825111932").toString(), "31U DQ 48251 11932");
    EXPECT_EQ(Mgrs::parse("31U DQ 48 11").toString(), "31U DQ 48000 11000");
}

TEST(mgrs_unittest, mgrs_to_utm_matches_reference_example)
{
    const UtmMgrs utm = Mgrs::parse("31U DQ 48251 11932").toUtm();

    EXPECT_EQ(utm.zone(), 31);
    EXPECT_EQ(utm.hemisphere(), Utm::Hemisphere::N);
    EXPECT_NEAR(utm.easting(), 448251.0, 1e-9);
    EXPECT_NEAR(utm.northing(), 5411932.0, 1e-9);
}

TEST(mgrs_unittest, utm_to_mgrs_matches_reference_example)
{
    const Mgrs mgrs = UtmMgrs(31, Utm::Hemisphere::N, 448251.795, 5411932.678).toMgrs();

    EXPECT_EQ(mgrs.toString(), "31U DQ 48251 11932");
}

TEST(mgrs_unittest, latlon_to_mgrs_matches_greenwich_reference_example)
{
    const Mgrs mgrs = LatLonUtmMgrs(51.47788055555556, -0.0014694444444444444).toUtm().toMgrs();

    EXPECT_EQ(mgrs.toString(), "30U YC 08215 07233");
}

TEST(mgrs_unittest, round_trip_preserves_grid_square_southwest_corner)
{
    const Mgrs mgrs = LatLonUtmMgrs(-33.857, 151.215).toUtm().toMgrs();
    const UtmMgrs utm = Mgrs::parse(mgrs.toString()).toUtm();

    EXPECT_EQ(utm.zone(), 56);
    EXPECT_EQ(utm.hemisphere(), Utm::Hemisphere::S);
    EXPECT_NEAR(std::fmod(utm.easting(), 1.0), 0.0, 1e-9);
    EXPECT_NEAR(std::fmod(utm.northing(), 1.0), 0.0, 1e-9);
    EXPECT_EQ(utm.toMgrs().toString(), mgrs.toString());
}

TEST(mgrs_unittest, rejects_invalid_constructor_inputs)
{
    EXPECT_THROW(Mgrs(0, 'U', 'D', 'Q', 48251.0, 11932.0), std::invalid_argument);
    EXPECT_THROW(Mgrs(61, 'U', 'D', 'Q', 48251.0, 11932.0), std::invalid_argument);
    EXPECT_THROW(Mgrs(31, 'I', 'D', 'Q', 48251.0, 11932.0), std::invalid_argument);
    EXPECT_THROW(Mgrs(31, 'U', 'I', 'Q', 48251.0, 11932.0), std::invalid_argument);
    EXPECT_THROW(Mgrs(31, 'U', 'D', 'O', 48251.0, 11932.0), std::invalid_argument);
    EXPECT_THROW(Mgrs(32, 'X', 'P', 'F', 0.0, 0.0), std::invalid_argument);
    EXPECT_THROW(Mgrs(31, 'U', 'D', 'Q', -1.0, 11932.0), std::invalid_argument);
    EXPECT_THROW(Mgrs(31, 'U', 'D', 'Q', 100000.0, 11932.0), std::invalid_argument);
    EXPECT_THROW(Mgrs(31, 'U', 'D', 'Q', 48251.0, std::numeric_limits<double>::quiet_NaN()),
        std::invalid_argument);
}

TEST(mgrs_unittest, rejects_invalid_parse_inputs)
{
    EXPECT_THROW((void)Mgrs::parse(""), std::invalid_argument);
    EXPECT_THROW((void)Mgrs::parse("31U DQ 48251"), std::invalid_argument);
    EXPECT_THROW((void)Mgrs::parse("31U DI 48251 11932"), std::invalid_argument);
    EXPECT_THROW((void)Mgrs::parse("31U DQ 4825 11932"), std::invalid_argument);
    EXPECT_THROW((void)Mgrs::parse("31UDQ482511932"), std::invalid_argument);
    EXPECT_THROW((void)Mgrs::parse("31u DQ 48251 11932"), std::invalid_argument);
    EXPECT_THROW((void)Mgrs::parse("32X PF 00000 00000"), std::invalid_argument);
}

TEST(mgrs_unittest, rejects_invalid_format_precision)
{
    const Mgrs mgrs(31, 'U', 'D', 'Q', 48251.0, 11932.0);

    EXPECT_THROW((void)mgrs.toString(1), std::invalid_argument);
    EXPECT_THROW((void)mgrs.toString(11), std::invalid_argument);
}

}
