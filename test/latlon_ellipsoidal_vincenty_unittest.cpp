#include <gtest/gtest.h>

#include "geodesy/latlon_ellipsoidal_vincenty.h"
#include "geodesy/latlon_ellipsoidal_datum.h"
#include "geodesy/dms.h"

using Vincenty = geodesy::LatLonEllipsoidalVincenty;

class latlon_ellipsoidal_vincenty_unittest : public testing::Test
{
protected:
   void SetUp() override
   {
      geodesy::Dms::setSeparator("");
   }
   void TearDown() override
   {
      
   }

public:
   const double m_circEquatorial = 40075016.686; // eslint-disable-line no-unused-vars
   const double m_circMeridional = 40007862.918;
};

TEST_F(latlon_ellipsoidal_vincenty_unittest, example)
{
   EXPECT_EQ("", geodesy::Dms::separator());
   EXPECT_NEAR(Vincenty(50.06632, -5.71475).distanceTo({58.64402, -3.07009}), 969954.166, 0.001);
   EXPECT_NEAR(Vincenty(50.06632, -5.71475).initialBearingTo({58.64402, -3.07009}), 9.1419, 0.0001);
   EXPECT_NEAR(Vincenty(50.06632, -5.71475).finalBearingTo({58.64402, -3.07009}), 11.2972, 0.0001);
   EXPECT_NEAR(Vincenty(-37.95103, 144.42487).finalBearingOn(54972.271, 306.86816), 307.1736, 0.0001);

   EXPECT_EQ(Vincenty(-37.95103, 144.42487).destinationPoint(54972.271, 306.86816).toString(), "37.6528°S, 143.9265°E");
   EXPECT_EQ(Vincenty(50.06632, -5.71475).intermediatePointTo({58.64402, -3.07009}, 0.5).toString(), "54.3639°N, 004.5304°W");
}

TEST_F(latlon_ellipsoidal_vincenty_unittest, Uk)
{
   const Vincenty le(50.06632, -5.71475);
   const Vincenty jog(58.64402, -3.07009);
   const auto dist = 969954.166, brngInit = 9.1418775, brngFinal = 11.2972204;
   EXPECT_NEAR(le.distanceTo(jog), dist, 0.001);
   EXPECT_NEAR(le.initialBearingTo(jog), brngInit, 0.0000001);
   EXPECT_NEAR(le.finalBearingTo(jog), brngFinal, 0.0000001);
   EXPECT_EQ(le.destinationPoint(dist, brngInit).toString(geodesy::Dms::D), jog.toString(geodesy::Dms::D));
   EXPECT_NEAR(le.finalBearingOn(dist, brngInit), brngFinal, 0.0000001);
   EXPECT_EQ(le.intermediatePointTo(jog, 0), le);
   EXPECT_TRUE(le.intermediatePointTo(jog, 1).equals(jog));
}