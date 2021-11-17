//
// Created by Bin on 2021/11/15.
//

#ifndef ELLIPSOIDS_H
#define ELLIPSOIDS_H

namespace geodesy
{
   /*
   * Ellipsoid parameters; exposed through static getter below.
   */
   struct Ellipsoid
   {
      double a = 0.0;
      double b = 0.0;
      double f = 0.0;

      explicit operator bool() const
      {
        return std::fabs(a) <= std::numeric_limits<double>::epsilon() &&
           std::fabs(b) <= std::numeric_limits<double>::epsilon() &&
           std::fabs(f) <= std::numeric_limits<double>::epsilon();
      }
   };

   struct Transforms
   {
      double tx = 0.0;
      double ty = 0.0;
      double tz = 0.0;
      double s  = 0.0;
      double rx = 0.0;
      double ry = 0.0;
      double rz = 0.0;

      explicit operator bool() const
      {
        return std::fabs(tx) <= std::numeric_limits<double>::epsilon() &&
               std::fabs(ty) <= std::numeric_limits<double>::epsilon() &&
               std::fabs(tz) <= std::numeric_limits<double>::epsilon() &&
               std::fabs(s)  <= std::numeric_limits<double>::epsilon() &&
               std::fabs(rx) <= std::numeric_limits<double>::epsilon() &&
               std::fabs(ry) <= std::numeric_limits<double>::epsilon() &&
               std::fabs(rz) <= std::numeric_limits<double>::epsilon();
      }
   };

   /*
   * A set of Ellipsoid
   */
   struct Ellipsoids
   {
      Ellipsoid WGS84;
      Ellipsoid Airy1830;
      Ellipsoid AiryModified;
      Ellipsoid Bessel1841;
      Ellipsoid Clarke1866;
      Ellipsoid Clarke1880IGN;
      Ellipsoid GRS80;
      Ellipsoid Intl1924; // aka Hayford
      Ellipsoid WGS72;
   };


   static const Ellipsoids s_ellipsoids = {
      /* WGS84         */ { 6378137,      6356752.314245, 1 / 298.257223563 },
      /* Airy1830      */ { 6377563.396,  6356256.909,    1 / 299.3249646   },
      /* AiryModified  */ { 6377340.189,  6356034.448,    1 / 299.3249646   },
      /* Bessel1841    */ { 6377397.155,  6356078.962818, 1 / 299.1528128   },
      /* Clarke1866    */ { 6378206.4,    6356583.8,      1 / 294.978698214 },
      /* Clarke1880IGN */ { 6378249.2,    6356515.0,      1 / 293.466021294 },
      /* GRS80         */ { 6378137,      6356752.314140, 1 / 298.257222101 },
      /* Intl1924      */ { 6378388,      6356911.946,    1 / 297           }, // aka Hayford
      /* WGS72         */ { 6378135,      6356750.5,      1 / 298.26        }
   };

   /*
   * Datums; exposed through static getter below.
   */
   struct Datum
   {
      Ellipsoid ellipsoid;
      // transforms: t in metres, s in ppm, r in arcseconds
      Transforms transforms;

      operator bool() const
      {
         return ellipsoid && transforms;
      }
   };

   /*
   * A set of Datum
   */
   struct Datums
   {
      Datum ED50;
      Datum ETRS89;
      Datum rl1975;
      Datum NAD27;
      Datum NAD83;
      Datum NTF;
      Datum OSGB36;
      Datum Potsdam;
      Datum TokyoJapan;
      Datum WGS72;
      Datum WGS84;
   };

   static const Datums s_datums = {
      /* ED50       */ { s_ellipsoids.Intl1924,      {   89.5,    93.8,    123.1,    -1.2,     0.0,      0.0,      0.156    } }, // epsg.io/1311
      /* ETRS89     */ { s_ellipsoids.GRS80,         {    0,       0,        0,       0,       0,        0,        0        } }, // epsg.io/1149; @ 1-metre level
      /* rl1975     */ { s_ellipsoids.AiryModified,  { -482.530, 130.596, -564.557,  -8.150,   1.042,    0.214,    0.631    } }, // epsg.io/1954
      /* NAD27      */ { s_ellipsoids.Clarke1866,    {    8,    -160,     -176,       0,       0,        0,        0        } },
      /* NAD83      */ { s_ellipsoids.GRS80,         {    0.9956, -1.9103,  -0.5215, -0.00062, 0.025915, 0.009426, 0.011599 } },
      /* NTF        */ { s_ellipsoids.Clarke1880IGN, {  168,      60,     -320,       0,       0,        0,        0        } },
      /* OSGB36     */ { s_ellipsoids.Airy1830,      { -446.448, 125.157, -542.060,  20.4894, -0.1502,  -0.2470,  -0.8421   } }, // epsg.io/1314
      /* Potsdam    */ { s_ellipsoids.Bessel1841,    { -582,    -105,     -414,      -8.3,     1.04,     0.35,    -3.08     } },
      /* TokyoJapan */ { s_ellipsoids.Bessel1841,    {  148,    -507,     -685,       0,       0,        0,        0        } },
      /* WGS72      */ { s_ellipsoids.WGS72,         {    0,       0,       -4.5,    -0.22,    0,        0,        0.554    } },
      /* WGS84      */ { s_ellipsoids.WGS84,         {    0.0,     0.0,      0.0,     0.0,     0.0,      0.0,      0.0      } }
   };
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
      float epoch = std::numeric_limits<float>::epsilon();
      Ellipsoid ellipsoid;

      operator bool() const
      {
         return !name.empty() && std::fabs(epoch - 0) <= std::numeric_limits<float>::epsilon() && ellipsoid;
      }
   };

   struct ReferenceFrames
   {
      ReferenceFrame ITRF2014;
      ReferenceFrame ITRF2008;
      ReferenceFrame ITRF2005;
      ReferenceFrame ITRF2000;
      ReferenceFrame ITRF93;
      ReferenceFrame ITRF91;
      ReferenceFrame WGS84g1762;
      ReferenceFrame WGS84g1674;
      ReferenceFrame WGS84g1150;
      ReferenceFrame ETRF2000;
      ReferenceFrame NAD83;
      ReferenceFrame GDA94;
   };

   /*
    * Reference frames; exposed through static getter below.
    */
   static const ReferenceFrames s_reference_frames = {
      /* ITRF2014   */ { "ITRF2014",   2010.0, s_ellipsoids.GRS80 },
      /* ITRF2008   */ { "ITRF2008",   2005.0, s_ellipsoids.GRS80 },
      /* ITRF2005   */ { "ITRF2005",   2000.0, s_ellipsoids.GRS80 },
      /* ITRF2000   */ { "ITRF2000",   1997.0, s_ellipsoids.GRS80 },
      /* ITRF93     */ { "ITRF93",     1988.0, s_ellipsoids.GRS80 },
      /* ITRF91     */ { "ITRF91",     1988.0, s_ellipsoids.GRS80 },
      /* WGS84g1762 */ { "WGS84g1762", 2005.0, s_ellipsoids.WGS84 },
      /* WGS84g1674 */ { "WGS84g1674", 2005.0, s_ellipsoids.WGS84 },
      /* WGS84g1150 */ { "WGS84g1150", 2001.0, s_ellipsoids.WGS84 },
      /* ETRF2000   */ { "ETRF2000",   2005.0, s_ellipsoids.GRS80 }, // ETRF2000(R08)
      /* NAD83      */ { "NAD83",      1997.0, s_ellipsoids.GRS80 }, // CORS96
      /* GDA94      */ { "GDA94",      1994.0, s_ellipsoids.GRS80 }
   };

}

#endif //ELLIPSOIDS_H
