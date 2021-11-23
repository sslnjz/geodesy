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

using namespace geodesy;

LatLonEllipsoidalReferenceFrame::LatLonEllipsoidalReferenceFrame(
    double lat, double lon, double height,
    std::optional<ReferenceFrame> referenceFrame,
    std::optional<std::string> epoch)
    : LatLonEllipsoidal(lat, lon, height)
{
    if (!referenceFrame || referenceFrame.value().epoch == std::nullopt)
        throw std::runtime_error("unrecognised reference frame");
    if (epoch.has_value() && std::isnan(std::stod(epoch.value())))
        throw std::runtime_error("invalid epoch");

    m_referenceFrame = referenceFrame;
    try
    {
        m_epoch = std::stof(epoch.value());
    }
    catch (const std::exception& e)
    {
        throw e;
    }
}

std::optional<ReferenceFrame> LatLonEllipsoidalReferenceFrame::referenceFrame()
{
    return m_referenceFrame;
}

std::optional<float> LatLonEllipsoidalReferenceFrame::epoch()
{
    return m_epoch || m_referenceFrame.value().epoch;
}

Ellipsoids LatLonEllipsoidalReferenceFrame::ellipsoids()
{
    return g_ellipsoids;
}

ReferenceFrames LatLonEllipsoidalReferenceFrame::referenceFrames()
{
    return g_reference_frames;
}

std::map<std::string, Helmert> LatLonEllipsoidalReferenceFrame::transformParameters()
{
    return s_txParams;
}
