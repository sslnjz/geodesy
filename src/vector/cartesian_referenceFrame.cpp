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

#include "cartesian_referenceFrame.h"
#include "latlon_ellipsoidal_referenceframe.h"
#include "strutil.h"

using namespace geodesy;

CartesianReferenceFrame::CartesianReferenceFrame() {

}

CartesianReferenceFrame::CartesianReferenceFrame(double x, double y, double z,
                                                 std::optional<ReferenceFrame> referenceFrame,
                                                 std::optional<std::string> epoch)
      : Cartesian(x, y, z) {
   if (referenceFrame.has_value() && !referenceFrame.value().epoch.has_value())
      throw std::invalid_argument("unrecognised reference frame");

   if (!epoch.has_value())
      throw std::invalid_argument("invalid epoch");


   _epoch = epoch.value();
   if (referenceFrame.has_value())
      _referenceFrame = referenceFrame.value();

}

LatLonEllipsoidalReferenceFrame CartesianReferenceFrame::toLatLon() {
   if (_referenceFrame == std::nullopt)
      throw std::runtime_error("cartesian reference frame not defined");

   const auto latLon = Cartesian::toLatLon(_referenceFrame->ellipsoid);
   const auto point = LatLonEllipsoidalReferenceFrame(latLon.lat(), latLon.lon(),
                                                          latLon.height(),
                                                          _referenceFrame,
                                                          _epoch);
   return point;
}

CartesianReferenceFrame CartesianReferenceFrame::convertReferenceFrame(const ReferenceFrame &to)
{
   if (to.epoch == std::nullopt)
      throw std::invalid_argument("unrecognised reference frame");
   if (!_referenceFrame)
      throw std::runtime_error("cartesian coordinate has no reference frame");

   if (_referenceFrame->name == to.name)
      return *this; // no-op!


   const auto oldTrf = _referenceFrame;
   const auto newTrf = to;

   // WGS84(G730/G873/G1150) are coincident with ITRF at 10-centimetre level; WGS84(G1674) and
   // ITRF20014 / ITRF2008 ‘are likely to agree at the centimeter level’ (QPS)
   if (strutil::start_with(oldTrf->name, "ITRF") && strutil::start_with(newTrf.name, "WGS84")) return *this;
   if (strutil::start_with(oldTrf->name, "WGS84") && strutil::start_with(newTrf.name, "ITRF")) return *this;

   const CartesianReferenceFrame oldC = *this;
   CartesianReferenceFrame newC;

   // is requested transformation available in single step?
   const auto txFwd = s_txParams.find(oldTrf->name+"→"+newTrf.name);
   const auto txRev = s_txParams.find(newTrf.name+"→"+oldTrf->name);

   auto reverseTransform = [](const Helmert& tx) -> Helmert
   {
      Helmert helmert;
      helmert.epoch = tx.epoch;
      for(int i = 0; i < 7; ++i) helmert.params[i] = -tx.params[i];
      for(int i = 0; i < 7; ++i) helmert.rates[i] = -tx.rates[i];

      return helmert;
   };

//   if (txFwd != s_txParams.end() || txRev != s_txParams.end())
//   {
//      // yes, single step available (either forward or reverse)
//      const auto tx = txFwd != s_txParams.end() ? *txFwd : reverseTransform(*txRev);
//      const t = this.epoch || this.referenceFrame.epoch;
//      const t0 = tx.epoch;//epoch || newTrf.epoch;
//      newC = oldC.applyTransform(tx.params, tx.rates, t - t0); // ...apply transform...
//   }
//   else
//   {
//      // find intermediate transform common to old & new to chain though; this is pretty yucky,
//      // but since with current transform params we can transform in no more than 2 steps, it works!
//      // TODO: find cleaner method!
//      const auto txAvailFromOld = Object.keys(txParams).filter(tx = > tx.split('→')[0] == oldTrf.name).map(tx = >
//                                                                                                  tx.split(
//                                                                                                        '→')[1]);
//      const txAvailToNew = Object.keys(txParams).filter(tx = > tx.split('→')[1] == newTrf.name).map(tx = >
//                                                                                                tx.split(
//                                                                                                      '→')[0]);
//      const txIntermediateFwd = txAvailFromOld.filter(tx = > txAvailToNew.includes(tx))[0];
//      const txAvailFromNew = Object.keys(txParams).filter(tx = > tx.split('→')[0] == newTrf.name).map(tx = >
//                                                                                                  tx.split(
//                                                                                                        '→')[1]);
//      const txAvailToOld = Object.keys(txParams).filter(tx = > tx.split('→')[1] == oldTrf.name).map(tx = >
//                                                                                                tx.split(
//                                                                                                      '→')[0]);
//      const txIntermediateRev = txAvailFromNew.filter(tx = > txAvailToOld.includes(tx))[0];
//      const txFwd1 = txParams[oldTrf.name + '→' + txIntermediateFwd];
//      const txFwd2 = txParams[txIntermediateFwd + '→' + newTrf.name];
//      const txRev1 = txParams[newTrf.name + '→' + txIntermediateRev];
//      const txRev2 = txParams[txIntermediateRev + '→' + oldTrf.name];
//      const tx1 = txIntermediateFwd ? txFwd1 : reverseTransform(txRev2);
//      const tx2 = txIntermediateFwd ? txFwd2 : reverseTransform(txRev1);
//      const t = this.epoch || this.referenceFrame.epoch;
//      newC = oldC.applyTransform(tx1.params, tx1.rates, t - tx1.epoch); // ...apply transform 1...
//      newC = newC.applyTransform(tx2.params, tx2.rates, t - tx2.epoch); // ...apply transform 2...
//   }
//
//   newC.referenceFrame = toReferenceFrame;
//   newC.epoch = oldC.epoch;
//
//   return newC;
//
//   function reverseTransform(tx) {
//      return {epoch: tx.epoch, params: tx.params.map(p = > -p), rates: tx.rates.map(r => -r)};
//   }
}


