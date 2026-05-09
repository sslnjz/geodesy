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

#include "geodesy/utm.h"

#include <cmath>
#include <limits>
#include <stdexcept>

#include <gtest/gtest.h>

namespace
{

using geodesy::Utm;

TEST(utm_unittest, value_type_preserves_fields_and_optional_projection_metadata)
{
    const Utm utm(31, Utm::Hemisphere::N, 448251.795, 5411932.678,
        geodesy::LatLonEllipsoidal::datums().WGS84, -0.531672641, 0.999632879);

    EXPECT_EQ(utm.zone(), 31);
    EXPECT_EQ(utm.hemisphere(), Utm::Hemisphere::N);
    EXPECT_NEAR(utm.easting(), 448251.795, 1e-9);
    EXPECT_NEAR(utm.northing(), 5411932.678, 1e-9);
    ASSERT_TRUE(utm.convergence().has_value());
    ASSERT_TRUE(utm.scale().has_value());
    EXPECT_NEAR(*utm.convergence(), -0.531672641, 1e-12);
    EXPECT_NEAR(*utm.scale(), 0.999632879, 1e-12);
}

TEST(utm_unittest, parse_accepts_space_separated_coordinate)
{
    const Utm utm = Utm::parse("31 N 448251.795 5411932.678");

    EXPECT_EQ(utm.zone(), 31);
    EXPECT_EQ(utm.hemisphere(), Utm::Hemisphere::N);
    EXPECT_NEAR(utm.easting(), 448251.795, 1e-9);
    EXPECT_NEAR(utm.northing(), 5411932.678, 1e-9);
}

TEST(utm_unittest, parse_accepts_lowercase_hemisphere)
{
    const Utm utm = Utm::parse("48 s 377298.745 8516965.206");

    EXPECT_EQ(utm.zone(), 48);
    EXPECT_EQ(utm.hemisphere(), Utm::Hemisphere::S);
    EXPECT_NEAR(utm.easting(), 377298.745, 1e-9);
    EXPECT_NEAR(utm.northing(), 8516965.206, 1e-9);
}

TEST(utm_unittest, to_string_formats_zone_hemisphere_easting_and_northing)
{
    const Utm utm(31, Utm::Hemisphere::N, 448251.0, 5411932.0);

    EXPECT_EQ(utm.toString(), "31 N 448251 5411932");
    EXPECT_EQ(utm.toString(4), "31 N 448251.0000 5411932.0000");
}

TEST(utm_unittest, to_string_pads_single_digit_zone)
{
    const Utm utm(3, Utm::Hemisphere::S, 500000.0, 5000000.0);

    EXPECT_EQ(utm.toString(3), "03 S 500000.000 5000000.000");
}

TEST(utm_unittest, rejects_invalid_value_type_inputs)
{
    EXPECT_THROW(Utm(0, Utm::Hemisphere::N, 448251.0, 5411932.0), std::invalid_argument);
    EXPECT_THROW(Utm(61, Utm::Hemisphere::N, 448251.0, 5411932.0), std::invalid_argument);
    EXPECT_THROW(Utm(31, Utm::Hemisphere::N, -1.0, 5411932.0), std::invalid_argument);
    EXPECT_THROW(Utm(31, Utm::Hemisphere::N, 1000000.001, 5411932.0), std::invalid_argument);
    EXPECT_THROW(Utm(31, Utm::Hemisphere::N, 448251.0, -1.0), std::invalid_argument);
    EXPECT_THROW(Utm(31, Utm::Hemisphere::N, 448251.0, 9329006.0), std::invalid_argument);
    EXPECT_THROW(Utm(31, Utm::Hemisphere::S, 448251.0, 1116913.999), std::invalid_argument);
    EXPECT_THROW(Utm(31, Utm::Hemisphere::S, 448251.0, 10000000.0), std::invalid_argument);
    EXPECT_THROW(Utm(31, Utm::Hemisphere::N, std::numeric_limits<double>::quiet_NaN(), 5411932.0),
        std::invalid_argument);
}

TEST(utm_unittest, rejects_invalid_parse_inputs)
{
    EXPECT_THROW(Utm::parse("31 448251 5411932"), std::invalid_argument);
    EXPECT_THROW(Utm::parse("31 Q 448251 5411932"), std::invalid_argument);
    EXPECT_THROW(Utm::parse("31 N x 5411932"), std::invalid_argument);
    EXPECT_THROW(Utm::parse("31 N 448251 y"), std::invalid_argument);
}

TEST(utm_unittest, verify_flag_allows_extended_coordinate_values)
{
    const Utm utm(31, Utm::Hemisphere::N, -10.0, -20.0, geodesy::LatLonEllipsoidal::datums().WGS84,
        std::nullopt, std::nullopt, false);

    EXPECT_NEAR(utm.easting(), -10.0, 1e-12);
    EXPECT_NEAR(utm.northing(), -20.0, 1e-12);
}

}
