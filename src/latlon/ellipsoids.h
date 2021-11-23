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

#ifndef ELLIPSOIDS_H
#define ELLIPSOIDS_H

#include <numeric>
#include <cmath>
#include <string>
#include <optional>

namespace geodesy
{
   /*
   * Ellipsoid parameters; exposed through static getter below.
   */
   struct Ellipsoid
   {
      double a{ 0.0 };
      double b{ 0.0 };
      double f{ 0.0 };

      inline bool operator==(const Ellipsoid& rhs) const
      {
          return std::fabs(a- rhs.a) <= std::numeric_limits<double>::epsilon() &&
            std::fabs(b- rhs.b) <= std::numeric_limits<double>::epsilon() &&
            std::fabs(f- rhs.f) <= std::numeric_limits<double>::epsilon();
      }

       inline bool operator!=(const Ellipsoid& rhs) const
       {
          return !(*this == rhs);
       }
   };

   /*
   * A set of Ellipsoid
   */
   static const struct Ellipsoids
   {
      Ellipsoid WGS84         { 6378137,      6356752.314245, 1 / 298.257223563 };
      Ellipsoid Airy1830      { 6377563.396,  6356256.909,    1 / 299.3249646   };
      Ellipsoid AiryModified  { 6377340.189,  6356034.448,    1 / 299.3249646   };
      Ellipsoid Bessel1841    { 6377397.155,  6356078.962818, 1 / 299.1528128   };
      Ellipsoid Clarke1866    { 6378206.4,    6356583.8,      1 / 294.978698214 };
      Ellipsoid Clarke1880IGN { 6378249.2,    6356515.0,      1 / 293.466021294 };
      Ellipsoid GRS80         { 6378137,      6356752.314140, 1 / 298.257222101 };
      Ellipsoid Intl1924      { 6378388,      6356911.946,    1 / 297.000000000 };             // aka Hayford
      Ellipsoid WGS72         { 6378135,      6356750.5,      1 / 298.26        };
   } g_ellipsoids;


   struct Transform
   {
      double tx{ 0.0 };
      double ty{ 0.0 };
      double tz{ 0.0 };
      double s { 0.0 };
      double rx{ 0.0 };
      double ry{ 0.0 };
      double rz{ 0.0 };

      inline bool operator==(const Transform& rhs) const
      {
         return std::fabs(tx - rhs.tx) <= std::numeric_limits<double>::epsilon() &&
            std::fabs(ty - rhs.ty) <= std::numeric_limits<double>::epsilon() &&
            std::fabs(tz - rhs.tz) <= std::numeric_limits<double>::epsilon() &&
            std::fabs(s - rhs.s) <= std::numeric_limits<double>::epsilon() &&
            std::fabs(rx - rhs.rx) <= std::numeric_limits<double>::epsilon() &&
            std::fabs(ry - rhs.ry) <= std::numeric_limits<double>::epsilon() &&
            std::fabs(rz - rhs.rz) <= std::numeric_limits<double>::epsilon();
      }

      inline bool operator!=(const Transform& rhs) const
      {
         return !(*this == rhs);
      }

      Transform& inverse()
      {
         tx = -tx;
         ty = -ty;
         tz = -tz;
         s = -s;
         rx = -rx;
         ry = -ry;
         rz = -rz;

         return *this;
      }
   };

   /*
   * Datums; exposed through static getter below.
   */
   struct Datum
   {
      Ellipsoid ellipsoid;
      // transforms: t in metres, s in ppm, r in arcseconds
      Transform transforms;

      inline bool operator==(const Datum& d) const
      {
          return (ellipsoid == d.ellipsoid) && (transforms == d.transforms);
      }

      inline bool operator!=(const Datum& d) const
      {
         return !(*this == d);
      }
   };

   static const struct Datums
   {
      Datum ED50       { g_ellipsoids.Intl1924,      {   89.5,    93.8,    123.1,    -1.2,     0.0,      0.0,      0.156    } };
      Datum ETRS89     { g_ellipsoids.GRS80,         {    0,       0,        0,       0,       0,        0,        0        } };
      Datum rl1975     { g_ellipsoids.AiryModified,  { -482.530, 130.596, -564.557,  -8.150,   1.042,    0.214,    0.631    } };
      Datum NAD27      { g_ellipsoids.Clarke1866,    {    8,    -160,     -176,       0,       0,        0,        0        } };
      Datum NAD83      { g_ellipsoids.GRS80,         {    0.9956, -1.9103,  -0.5215, -0.00062, 0.025915, 0.009426, 0.011599 } };
      Datum NTF        { g_ellipsoids.Clarke1880IGN, {  168,      60,     -320,       0,       0,        0,        0        } };
      Datum OSGB36     { g_ellipsoids.Airy1830,      { -446.448, 125.157, -542.060,  20.4894, -0.1502,  -0.2470,  -0.8421   } };
      Datum Potsdam    { g_ellipsoids.Bessel1841,    { -582,    -105,     -414,      -8.3,     1.04,     0.35,    -3.08     } };
      Datum TokyoJapan { g_ellipsoids.Bessel1841,    {  148,    -507,     -685,       0,       0,        0,        0        } };
      Datum WGS72      { g_ellipsoids.WGS72,         {    0,       0,       -4.5,    -0.22,    0,        0,        0.554    } };
      Datum WGS84      { g_ellipsoids.WGS84,         {    0.0,     0.0,      0.0,     0.0,     0.0,      0.0,      0.0      } };
   } &g_datums = *new Datums;
   /* sources:
    * - ED50:       www.gov.uk/guidance/oil-and-gas-petroleum-operations-notices#pon-4
    * - Irl1975:    www.osi.ie/wp-content/uploads/2015/05/transformations_booklet.pdf
    * - NAD27:      en.wikipedia.org/wiki/Helmert_transformation
    * - NAD83:      www.uvm.edu/giv/resources/WGS84_NAD83.pdf [strictly, WGS84(G1150) -> NAD83(CORS96) @ epoch 1997.0]
    *               (note NAD83(1986) ≡ WGS84(Original); confluence.qps.nl/pages/viewpage.action?pageId=29855173)
    * - NTF:        Nouvelle Triangulation Francaise geodesie.ign.fr/contenu/fichiers/Changement_systeme_geodesique.pdf
    * - OSGB36:     www.ordnancesurvey.co.uk/docs/support/guide-coordinate-systems-great-britain.pdf
    * - Potsdam:    kartoweb.itc.nl/geometrics/Coordinate%20transformations/coordtrans.html
    * - TokyoJapan: www.geocachingtoolbox.com?page=datumEllipsoidDetails
    * - WGS72:      www.icao.int/safety/pbn/documentation/eurocontrol/eurocontrol wgs 84 implementation manual.pdf
    *
    * more transform parameters are available from earth-info.nga.mil/GandG/coordsys/datums/NATO_DT.pdf,
    * www.fieldenmaps.info/cconv/web/cconv_params.js
    */
   /* note:
    * - ETRS89 reference frames are coincident with WGS-84 at epoch 1989.0 (ie null transform) at the one metre level.
    */



   /*
    * Reference frames; exposed through static getter below.
    */
   struct ReferenceFrame
   {
      std::string name;
      std::optional<float> epoch;
      Ellipsoid ellipsoid;

       inline bool operator==(const ReferenceFrame& rhs) const
       {
           return name == rhs.name && epoch == rhs.epoch && ellipsoid == rhs.ellipsoid;
       }

       inline bool operator!=(const ReferenceFrame& rhs) const
       {
           return !(*this == rhs);
       }
   };

   /*
    * Reference frames; exposed through static getter below.
    */
   static const struct ReferenceFrames
   {
      ReferenceFrame ITRF2014  { "ITRF2014",   2010.0, g_ellipsoids.GRS80 };
      ReferenceFrame ITRF2008  { "ITRF2008",   2005.0, g_ellipsoids.GRS80 };
      ReferenceFrame ITRF2005  { "ITRF2005",   2000.0, g_ellipsoids.GRS80 };
      ReferenceFrame ITRF2000  { "ITRF2000",   1997.0, g_ellipsoids.GRS80 };
      ReferenceFrame ITRF93    { "ITRF93",     1988.0, g_ellipsoids.GRS80 };
      ReferenceFrame ITRF91    { "ITRF91",     1988.0, g_ellipsoids.GRS80 };
      ReferenceFrame WGS84g1762{ "WGS84g1762", 2005.0, g_ellipsoids.WGS84 };
      ReferenceFrame WGS84g1674{ "WGS84g1674", 2005.0, g_ellipsoids.WGS84 };
      ReferenceFrame WGS84g1150{ "WGS84g1150", 2001.0, g_ellipsoids.WGS84 };
      ReferenceFrame ETRF2000  { "ETRF2000",   2005.0, g_ellipsoids.GRS80 }; // ETRF2000(R08)
      ReferenceFrame NAD83     { "NAD83",      1997.0, g_ellipsoids.GRS80 }; // CORS96
      ReferenceFrame GDA94     { "GDA94",      1994.0, g_ellipsoids.GRS80 };
   } &g_reference_frames = *new ReferenceFrames;

}

#endif //ELLIPSOIDS_H
