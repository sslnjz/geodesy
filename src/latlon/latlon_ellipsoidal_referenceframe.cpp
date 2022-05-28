/**********************************************************************************
*  MIT License                                                                    *
*                                                                                 *
*  Copyright (c) 2021 Binbin Song <ssln.jzs@gmail.com>                       *
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

#include "latlon_ellipsoidal_referenceframe.h"
#include "cartesian_referenceFrame.h"

using namespace geodesy;
using LatLonERF = LatLonEllipsoidalReferenceFrame;

LatLonERF::LatLonEllipsoidalReferenceFrame(
    double lat, double lon, double height,
    std::optional<ReferenceFrame> referenceFrame,
    std::optional<std::string> epoch)
    : LatLonEllipsoidal(lat, lon, height)
{
    if (!referenceFrame || !referenceFrame.value().epoch.has_value())
        throw std::runtime_error("unrecognized reference frame");
    if (epoch && 0.0f == std::strtof(epoch.value().c_str(), nullptr) )
        throw std::runtime_error("invalid epoch " + (epoch ? *epoch : ""));

    m_referenceFrame = referenceFrame;
    m_epoch = epoch;
}

std::optional<ReferenceFrame> LatLonERF::referenceFrame() const
{
    return m_referenceFrame;
}

std::optional<std::string> LatLonERF::epoch()
{
    return m_epoch ? *m_epoch : m_referenceFrame.value().epoch;
}

Ellipsoids LatLonERF::ellipsoids()
{
    return g_ellipsoids;
}

ReferenceFrames LatLonERF::referenceFrames()
{
    return g_reference_frames;
}

std::vector<HelmertTransforms> LatLonERF::transformParameters()
{
    return s_txParams;
}

LatLonEllipsoidalReferenceFrame LatLonERF::convertReferenceFrame(const ReferenceFrame &to) const
{
   if (to.epoch == std::nullopt)
      throw new std::invalid_argument("unrecognised reference frame");

   auto oldCartesian = toCartesian();                                   // convert geodetic to cartesian
   const auto newCartesian = oldCartesian.convertReferenceFrame(to); // convert TRF
   const auto newLatLon = newCartesian.toLatLon();                                 // convert cartesian back to to geodetic

   return newLatLon;
}

CartesianReferenceFrame LatLonERF::toCartesian() const
{
   const auto cartesian = LatLonEllipsoidal::toCartesian();
   const auto cartesianReferenceFrame = CartesianReferenceFrame(cartesian.x(), cartesian.y(), cartesian.z(),
                                 m_referenceFrame, m_epoch);
   return cartesianReferenceFrame;
}

LatLonERF LatLonERF::parse(double lat, double lon, double height,
                           std::optional<ReferenceFrame> referenceFrame,
                           std::optional<std::string> epoch)
{
   if (!referenceFrame || referenceFrame->epoch== std::nullopt)
      throw std::invalid_argument("unrecognised reference frame");

   const auto p = LatLonEllipsoidal::parse(lat, lon, height);
   return LatLonERF(p.lat(), p.lon(), p.height(), referenceFrame, epoch);
}

LatLonERF LatLonERF::parse(const std::string &dms, double height,
                           std::optional<ReferenceFrame> referenceFrame,
                           std::optional<std::string> epoch)
{
   if (!referenceFrame || referenceFrame->epoch== std::nullopt)
      throw std::invalid_argument("unrecognised reference frame");

   const auto p = LatLonEllipsoidal::parse(dms, height);
   return LatLonERF(p.lat(), p.lon(), p.height(), referenceFrame, epoch);
}

LatLonERF LatLonERF::parse(const std::string &lat, const std::string &lon, double height,
                           std::optional<ReferenceFrame> referenceFrame,
                           std::optional<std::string> epoch)
{
   if (!referenceFrame || referenceFrame->epoch== std::nullopt)
      throw std::invalid_argument("unrecognised reference frame");

   const auto p = LatLonEllipsoidal::parse(lat, lon, height);
   return LatLonERF(p.lat(), p.lon(), p.height(), referenceFrame, epoch);
}

LatLonERF LatLonERF::parse(const std::string &lat, const std::string &lon, std::string height,
                           std::optional<ReferenceFrame> referenceFrame,
                           std::optional<std::string> epoch)
{
   if (!referenceFrame || referenceFrame->epoch== std::nullopt)
      throw std::invalid_argument("unrecognised reference frame");

   const auto p = LatLonEllipsoidal::parse(lat, lon, height);
   return LatLonERF(p.lat(), p.lon(), p.height(), referenceFrame, epoch);
}
