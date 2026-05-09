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

#include "geodesy/latlon_utm.h"
#include "geodesy/utm.h"

#include <cmath>
#include <limits>
#include <stdexcept>

#include <gtest/gtest.h>

namespace
{

using geodesy::Utm;
using geodesy::LatLonUtm;

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
    EXPECT_THROW(Utm(31, Utm::Hemisphere::N, 448251.0, 5411932.0,
        geodesy::LatLonEllipsoidal::datums().WGS84, std::numeric_limits<double>::quiet_NaN(), 0.9996),
        std::invalid_argument);
    EXPECT_THROW(Utm(31, Utm::Hemisphere::N, 448251.0, 5411932.0,
        geodesy::LatLonEllipsoidal::datums().WGS84, 0.0, 0.0), std::invalid_argument);
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

TEST(utm_unittest, latlon_to_utm_matches_reference_example)
{
    const Utm utm = LatLonUtm(48.8582, 2.2945).toUtm();

    EXPECT_EQ(utm.zone(), 31);
    EXPECT_EQ(utm.hemisphere(), Utm::Hemisphere::N);
    EXPECT_NEAR(utm.easting(), 448251.795, 0.001);
    EXPECT_NEAR(utm.northing(), 5411932.678, 0.001);
    ASSERT_TRUE(utm.convergence().has_value());
    ASSERT_TRUE(utm.scale().has_value());
    EXPECT_NEAR(*utm.convergence(), -0.531312209, 1e-9);
    EXPECT_NEAR(*utm.scale(), 0.99963289743, 1e-12);
}

TEST(utm_unittest, utm_to_latlon_matches_reference_example)
{
    const LatLonUtm point = Utm(31, Utm::Hemisphere::N, 448251.795, 5411932.678).toLatLon();

    EXPECT_NEAR(point.lat(), 48.8582, 1e-8);
    EXPECT_NEAR(point.lon(), 2.2945, 1e-8);
    EXPECT_NEAR(point.getConvergence(), -0.531312211, 1e-9);
    EXPECT_NEAR(point.getScale(), 0.99963289743, 1e-12);
}

TEST(utm_unittest, utm_to_latlon_preserves_central_meridian_convergence_and_scale)
{
    const LatLonUtm point = Utm(31, Utm::Hemisphere::N, 500000.0, 0.0).toLatLon();

    EXPECT_NEAR(point.lat(), 0.0, 1e-12);
    EXPECT_NEAR(point.lon(), 3.0, 1e-12);
    EXPECT_NEAR(point.getConvergence(), 0.0, 1e-12);
    EXPECT_NEAR(point.getScale(), 0.9996, 1e-12);
}

TEST(utm_unittest, round_trip_preserves_northern_and_southern_coordinates)
{
    const LatLonUtm northern(51.47788, -0.00147, 46.0);
    const LatLonUtm northernRoundTrip = northern.toUtm().toLatLon();

    EXPECT_NEAR(northernRoundTrip.lat(), northern.lat(), 1e-9);
    EXPECT_NEAR(northernRoundTrip.lon(), northern.lon(), 1e-9);

    const LatLonUtm southern(-33.857, 151.215);
    const Utm southernUtm = southern.toUtm();
    const LatLonUtm southernRoundTrip = southernUtm.toLatLon();

    EXPECT_EQ(southernUtm.hemisphere(), Utm::Hemisphere::S);
    EXPECT_NEAR(southernRoundTrip.lat(), southern.lat(), 1e-9);
    EXPECT_NEAR(southernRoundTrip.lon(), southern.lon(), 1e-9);
}

TEST(utm_unittest, applies_norway_and_svalbard_zone_exceptions)
{
    EXPECT_EQ(LatLonUtm(60.0, 6.0).toUtm().zone(), 32);

    EXPECT_EQ(LatLonUtm(75.0, 5.0).toUtm().zone(), 31);
    EXPECT_EQ(LatLonUtm(75.0, 15.0).toUtm().zone(), 33);
    EXPECT_EQ(LatLonUtm(75.0, 25.0).toUtm().zone(), 35);
    EXPECT_EQ(LatLonUtm(75.0, 35.0).toUtm().zone(), 37);
}

TEST(utm_unittest, rejects_latitudes_outside_utm_limits_and_invalid_zone_override)
{
    EXPECT_THROW(LatLonUtm(84.1, 0.0).toUtm(), std::invalid_argument);
    EXPECT_THROW(LatLonUtm(-80.1, 0.0).toUtm(), std::invalid_argument);
    EXPECT_THROW(LatLonUtm(48.8582, 2.2945).toUtm(0), std::invalid_argument);
    EXPECT_THROW(LatLonUtm(48.8582, 2.2945).toUtm(61), std::invalid_argument);
}

TEST(utm_unittest, rejects_invalid_latlon_projection_metadata)
{
    LatLonUtm point(48.8582, 2.2945);

    EXPECT_THROW(point.setConvergence(std::numeric_limits<double>::infinity()), std::invalid_argument);
    EXPECT_THROW(point.setScale(std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
    EXPECT_THROW(point.setScale(0.0), std::invalid_argument);
}

}
