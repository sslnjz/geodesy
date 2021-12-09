/**********************************************************************************
*  MIT License                                                                    *
*                                                                                 *
*  Copyright (c) 2021 Binbin Song <ssln.jzs@gmail.com>                         *
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

#include "geodesy/latlon_ellipsoidal_datum.h"
#include "geodesy/cartesian_datum.h"

using LatLonE = geodesy::LatLonEllipsoidal;
using LEDatum = geodesy::LatLonEllipsoidalDatum;
using CDatum = geodesy::CartesianDatum;

TEST(latlon_ellipsoidal_datum_unittest, examples)
{
   geodesy::Dms::setSeparator("");
   EXPECT_EQ(LEDatum(53.3444, -6.2577, 17, LatLonE::datums().Irl1975)
      .toString(), "53.3444°N, 006.2577°W");
   EXPECT_EQ(LEDatum::parse("51.47736, 0.0000", 0, LEDatum::datums().OSGB36).toString(), "51.4774°N, 000.0000°E");
   EXPECT_EQ(LEDatum(51.47788, -0.00147).convertDatum(LEDatum::datums().OSGB36).toString(), "51.4774°N, 000.0001°E");
   EXPECT_EQ(CDatum(4027893.924, 307041.993, 4919474.294).toLatLon()
      .convertDatum(LEDatum::datums().OSGB36).toString(), "50.7971°N, 004.3612°E");
}

TEST(latlon_ellipsoidal_datum_unittest, getter)
{
   EXPECT_NO_THROW(LEDatum::ellipsoids().WGS84);
}

TEST(latlon_ellipsoidal_datum_unittest, convert_datum_Greenwich)
{
   const auto greenwichWGS84 = LEDatum(51.47788, -0.00147); // default WGS84
   const auto greenwichOSGB36 = greenwichWGS84.convertDatum(LEDatum::datums().OSGB36);
   // greenwichOSGB36.height = 0;
   EXPECT_EQ(greenwichOSGB36.toString(geodesy::Dms::D, 6), "51.477364°N, 000.000150°E"); // TODO: huh? should be 0°E? out by c. 10 metres / 0.5″! am I missing something?
   EXPECT_EQ(greenwichOSGB36.convertDatum(LEDatum::datums().WGS84).toString(geodesy::Dms::D, 5), "51.47788°N, 000.00147°W");
}

TEST(latlon_ellipsoidal_datum_unittest, convert_datum_Petroleum_Operations_Notices)
{
   // https://www.gov.uk/guidance/oil-and-gas-petroleum-operations-notices#test-point-using-osgb-petroleum-transformation-parameters
   EXPECT_EQ(LEDatum(53, 1, 50).convertDatum(LEDatum::datums().OSGB36).toString(geodesy::Dms::DMS, 3, 2), "52°59′58.719″N, 001°00′06.490″E +3.99m");
   // https://www.gov.uk/guidance/oil-and-gas-petroleum-operations-notices#test-point-using-common-offshore-transformation-parameters
   EXPECT_EQ(LEDatum(53, 1, 50).convertDatum(LEDatum::datums().ED50).toString(geodesy::Dms::DMS, 3, 2), "53°00′02.887″N, 001°00′05.101″E +2.72m");
   EXPECT_EQ(LEDatum(53, 1, 50).convertDatum(LEDatum::datums().OSGB36)
      .convertDatum(LEDatum::datums().ED50)
      .convertDatum(LEDatum::datums().WGS84).toString(geodesy::Dms::D, 4, 1), "53.0000°N, 001.0000°E +50.0m");
}

TEST(latlon_ellipsoidal_datum_unittest, equals)
{
   const auto p1 = LEDatum(51.47788, -0.00147, 1, LEDatum::datums().WGS84);
   const auto p2 = LEDatum(51.47788, -0.00147, 1, LEDatum::datums().WGS84);
   EXPECT_TRUE(p1 == p2);
   EXPECT_TRUE(p1.equals(p2));
   EXPECT_FALSE(p1.equals(LEDatum(0, -0.00147, 1, LEDatum::datums().WGS84)));
   EXPECT_FALSE(p1.equals(LEDatum(51.47788, 0, 1, LEDatum::datums().WGS84)));
   EXPECT_FALSE(p1.equals(LEDatum(51.47788, -0.00147, 99, LEDatum::datums().WGS84)));
   EXPECT_FALSE(p1.equals(LEDatum(51.47788, -0.00147, 1, LEDatum::datums().Irl1975)));
}

TEST(latlon_ellipsoidal_datum_unittest, cartesian)
{
   const auto p = LEDatum::parse("45N, 45E");
   EXPECT_EQ(p.toCartesian().toString(), "[3194419, 3194419, 4487348]");
   auto c = geodesy::CartesianDatum(3194419, 3194419, 4487348);
   EXPECT_EQ(c.toLatLon().toString(), "45.0000°N, 045.0000°E");
}

