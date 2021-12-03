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

#include "cartesian_referenceFrame.h"
#include "latlon_ellipsoidal_referenceframe.h"
#include "strutil.h"

using namespace geodesy;

CartesianReferenceFrame::CartesianReferenceFrame() = default;
CartesianReferenceFrame::CartesianReferenceFrame(double x, double y, double z,
                                                 std::optional<ReferenceFrame> referenceFrame,
                                                 std::optional<std::string> epoch)
      : Cartesian(x, y, z) {
   if (referenceFrame.has_value() && !referenceFrame.value().epoch.has_value())
      throw std::invalid_argument("recognized reference frame");

   if (!epoch.has_value())
      throw std::invalid_argument("invalid epoch");


   _epoch = epoch.value();
   if (referenceFrame.has_value())
      _referenceFrame = referenceFrame.value();

}

std::optional<ReferenceFrame> CartesianReferenceFrame::referenceFrame() const
{
   return _referenceFrame;
}

void CartesianReferenceFrame::setReferenceFrame(const ReferenceFrame& referenceFrame)
{
   if (!referenceFrame.epoch.has_value())
      throw std::invalid_argument("unrecognized reference frame");
   _referenceFrame = referenceFrame;
}

std::optional<std::string> CartesianReferenceFrame::epoch()
{
   return _epoch ? *_epoch : _referenceFrame ? _referenceFrame->epoch : std::nullopt;
}

void CartesianReferenceFrame::setEpoch(std::string epoch)
{
   if (_referenceFrame.has_value() &&
      _referenceFrame->epoch.has_value() &&
      _epoch != (*_referenceFrame).epoch.value())
   {
      _epoch = epoch;
   }
}

LatLonEllipsoidalReferenceFrame CartesianReferenceFrame::toLatLon() const
{
   if (_referenceFrame == std::nullopt)
      throw std::runtime_error("cartesian reference frame not defined");

   const auto latLon = Cartesian::toLatLon(_referenceFrame->ellipsoid);
   const auto point = LatLonEllipsoidalReferenceFrame(latLon.lat(), latLon.lon(),
                                                          latLon.height(),
                                                          _referenceFrame,
                                                          _epoch);
   return point;
}

CartesianReferenceFrame CartesianReferenceFrame::convertReferenceFrame(const ReferenceFrame &to) const
{
   if (to.epoch == std::nullopt)
      throw std::invalid_argument("unrecognized reference frame");
   if (!_referenceFrame)
      throw std::runtime_error("cartesian coordinate has no reference frame");

   if (_referenceFrame->name == to.name)
      return *this; // no-op!


   const auto oldTrf = *_referenceFrame;
   const auto newTrf = to;

   // WGS84(G730/G873/G1150) are coincident with ITRF at 10-centimetre level; WGS84(G1674) and
   // ITRF20014 / ITRF2008 ‘are likely to agree at the centimeter level’ (QPS)
   if (strutil::start_with(oldTrf.name, "ITRF") && strutil::start_with(newTrf.name, "WGS84")) return *this;
   if (strutil::start_with(oldTrf.name, "WGS84") && strutil::start_with(newTrf.name, "ITRF")) return *this;

   const CartesianReferenceFrame oldC = *this;
   CartesianReferenceFrame newC;

   const auto txFwd = std::find_if(s_txParams.begin(), s_txParams.end(), 
      [&](const auto& item) { return (item.from == oldTrf.name && item.to == newTrf.name);});
   const auto txRev = std::find_if(s_txParams.begin(), s_txParams.end(),
      [&](const auto& item){ return (item.from == newTrf.name && item.to == oldTrf.name);});

   auto reverseTransform = [](const Helmert& tx) -> Helmert
   {
      Helmert helmert;
      helmert.epoch = tx.epoch;
      for(int i = 0; i < 7; ++i) helmert.params[i] = -tx.params[i];
      for(int i = 0; i < 7; ++i) helmert.rates[i] = -tx.rates[i];

      return helmert;
   };

   if (txFwd != s_txParams.end() || txRev != s_txParams.end())
   {
      // yes, single step available (either forward or reverse)
      const auto tx = txFwd != s_txParams.end() ? txFwd->helmert : reverseTransform(txRev->helmert);
      const auto t = _epoch ? _epoch : _referenceFrame->epoch;
      const auto t0 = tx.epoch;//epoch || newTrf.epoch;
      newC = oldC.applyTransform(tx.params, tx.rates, std::stoi(*t) - std::stoi(t0)); // ...apply transform...
   }
   else
   {
      // find intermediate transform common to old & new to chain though; this is pretty yucky,
      // but since with current transform params we can transform in no more than 2 steps, it works!

      std::vector<HelmertTransforms> txAvailFromOld, txAvailToNew;
      std::for_each(s_txParams.begin(), s_txParams.end(), 
         [&](const HelmertTransforms& item)
         {
            if (item.from == oldTrf.name)
            {
               txAvailFromOld.emplace_back(item);
            }

            if (item.to == newTrf.name)
            {
               txAvailToNew.emplace_back(item);
            }
         });

      // search a suitable intermediate transform
      std::optional<std::pair<HelmertTransforms, HelmertTransforms>> intermediateTxFwd = std::nullopt, intermediateTxRev = std::nullopt;
      for (auto p = s_txParams.begin(); p != s_txParams.end(); ++p)
      {
         for (auto q = s_txParams.end() - 1; q != s_txParams.begin() && p != q; --q)
         {
            if(p->from == oldTrf.name && q->to == newTrf.name && p->to == q->from)
            {
               intermediateTxFwd->first = *p;
               intermediateTxFwd->second = *q;
            }

            if (p->from == newTrf.name && q->to == oldTrf.name && p->to == q->from)
            {
               intermediateTxRev->first = *p;
               intermediateTxRev->second = *q;
            }
         }
      }

      if(intermediateTxRev)
      {
         const auto tx1 = intermediateTxFwd ? intermediateTxFwd->first.helmert : reverseTransform(intermediateTxRev->second.helmert);
         const auto tx2 = intermediateTxFwd ? intermediateTxFwd->second.helmert : reverseTransform(intermediateTxRev->first.helmert);


         const auto t = _epoch ?  _epoch : _referenceFrame->epoch;
         newC = oldC.applyTransform(tx1.params, tx1.rates, std::stoi(*t) - std::stoi(tx1.epoch)); // ...apply transform 1...
         newC = newC.applyTransform(tx2.params, tx2.rates, std::stoi(*t) - std::stoi(tx2.epoch)); // ...app
      }
   }

   newC._referenceFrame = to;
   newC._epoch = oldC._epoch;

   return newC;
}

CartesianReferenceFrame CartesianReferenceFrame::applyTransform(const double params[7], const double rates[7], double deltat) const
{
   // this point
   const auto x1 = x(), y1 = y(), z1 = z();

   // base parameters
   const auto tx = params[0] / 1000; // x-shift: normalise millimetres to metres
   const auto ty = params[1] / 1000; // y-shift: normalise millimetres to metres
   const auto tz = params[2] / 1000; // z-shift: normalise millimetres to metres
   const auto s = params[3] / 1e9; // scale: normalise parts-per-billion
   const auto rx = toRadians((params[4] / 3600 / 1000)); // x-rotation: normalise milliarcseconds to radians
   const auto ry = toRadians((params[5] / 3600 / 1000)); // y-rotation: normalise milliarcseconds to radians
   const auto rz = toRadians((params[6] / 3600 / 1000)); // z-rotation: normalise milliarcseconds to radians

   // rate parameters
   const auto x_shift = rates[0] / 1000; // x-shift: normalise millimetres to metres
   const auto y_shift = rates[1] / 1000; // y-shift: normalise millimetres to metres
   const auto z_shift = rates[2] / 1000; // z-shift: normalise millimetres to metres
   const auto scale = rates[3] / 1e9; // scale: normalise parts-per-billion
   const auto x_rotation = toRadians((rates[4] / 3600 / 1000)); // x-rotation: normalise milliarcseconds to radians
   const auto y_rotation = toRadians((rates[5] / 3600 / 1000)); // y-rotation: normalise milliarcseconds to radians
   const auto z_rotation = toRadians((rates[6] / 3600 / 1000)); // z-rotation: normalise milliarcseconds to radians

   // combined (normalised) parameters
   const double T[] = {tx + x_shift * deltat, ty + y_shift * deltat, tz + z_shift * deltat};
   const double R[] = {rx + x_rotation * deltat, ry + y_rotation * deltat, rz + z_rotation * deltat};
   const auto S = 1 + s + scale * deltat;

   // apply transform (shift, scale, rotate)
   const auto x2 = T[0] + x1 * S - y1 * R[2] + z1 * R[1];
   const auto y2 = T[1] + x1 * R[2] + y1 * S - z1 * R[0];
   const auto z2 = T[2] - x1 * R[1] + y1 * R[0] + z1 * S;

   return CartesianReferenceFrame(x2, y2, z2);
}

std::string CartesianReferenceFrame::toString(int dp) const
{
   const auto epoch = _referenceFrame && _epoch != _referenceFrame->epoch ? *_epoch : "";
   const auto trf = _referenceFrame ? "(" + _referenceFrame->name + ((!epoch.empty()) ? ("@" + epoch) : epoch) : "";
   return Cartesian::toString(dp) + trf ;
}


