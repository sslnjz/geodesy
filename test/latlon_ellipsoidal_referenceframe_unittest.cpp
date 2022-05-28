/**********************************************************************************
*  MIT License                                                                    *
*                                                                                 *
*  Copyright (c) 2022 Binbin Song <ssln.jzs@gmail.com>                         *
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
#include "geodesy/latlon_ellipsoidal_referenceframe.h"
#include "geodesy/cartesian.h"
#include "geodesy/ellipsoids.h"

using LLEREF = geodesy::LatLonEllipsoidalReferenceFrame;

TEST(ReferenceFrame, constructor)
{
   LLEREF TRF(0, 0, 0, geodesy::g_reference_frames.ITRF2014, "2000.0");
   EXPECT_EQ(TRF.toString(), "00.0000°N, 000.0000°E");
   EXPECT_THROW({
      try {
         LLEREF badTRF(0, 0, 0, std::nullopt);
      }
      catch (const std::exception& e)
      {
         EXPECT_STREQ("unrecognized reference frame", e.what());
         throw;
      }
   }, std::runtime_error);
   EXPECT_THROW({
      try
      {
         LLEREF badTRF(0, 0, 0, geodesy::g_reference_frames.ITRF2014, "xxx");
      }
      catch (const std::exception& e)
      {
         EXPECT_STREQ("invalid epoch xxx", e.what());
         throw;
      }
   }, std::runtime_error);
}

TEST(ReferenceFrame, examples)
{
   LLEREF ctor(51.47788, -0.00147, 0, geodesy::g_reference_frames.ITRF2000);
   EXPECT_STREQ("51.4779°N, 000.0015°W", ctor.toString().c_str());
   EXPECT_STREQ(LLEREF::parse(51.47788, -0.00147, 17, LLEREF::referenceFrames().ETRF2000).toString().c_str(),
                "51.4779°N, 000.0015°W");
   EXPECT_STREQ(LLEREF::parse("51.47788, -0.00147", 17, LLEREF::referenceFrames().ETRF2000).toString().c_str(),
                "51.4779°N, 000.0015°W");
}