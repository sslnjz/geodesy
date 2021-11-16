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
    };

    struct Transforms
    {
        double tx = 0.0;
        double ty = 0.0;
        double tz = 0.0;
        double s = 0.0;
        double rx = 0.0;
        double ry = 0.0;
        double rz = 0.0;
    };

    /*
     * Datums; exposed through static getter below.
     */
    struct Datum
    {
        Ellipsoid ellipsoid;
        // transforms: t in metres, s in ppm, r in arcseconds
        Transforms transforms;
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

    static const Ellipsoids s_ellipsoids = {
        { 6378137,      6356752.314245, 1/298.257223563 },
        { 6377563.396,  6356256.909,    1/299.3249646   },
        { 6377340.189,  6356034.448,    1/299.3249646   },
        { 6377397.155,  6356078.962818, 1/299.1528128   },
        { 6378206.4,    6356583.8,      1/294.978698214 },
        { 6378249.2,    6356515.0,      1/293.466021294 },
        { 6378137,      6356752.314140, 1/298.257222101 },
        { 6378388,      6356911.946,    1/297           }, // aka Hayford
        { 6378135,      6356750.5,      1/298.26        }
    };

    static const Datums s_datums = {
        { s_ellipsoids.Intl1924,      {   89.5,    93.8,    123.1,    -1.2,     0.0,      0.0,      0.156    } }, // epsg.io/1311
        { s_ellipsoids.GRS80,         {    0,       0,        0,       0,       0,        0,        0        } }, // epsg.io/1149; @ 1-metre level
        { s_ellipsoids.AiryModified,  { -482.530, 130.596, -564.557,  -8.150,   1.042,    0.214,    0.631    } }, // epsg.io/1954
        { s_ellipsoids.Clarke1866,    {    8,    -160,     -176,       0,       0,        0,        0        } },
        { s_ellipsoids.GRS80,         {    0.9956, -1.9103,  -0.5215, -0.00062, 0.025915, 0.009426, 0.011599 } },
        { s_ellipsoids.Clarke1880IGN, {  168,      60,     -320,       0,       0,        0,        0        } },
        { s_ellipsoids.Airy1830,      { -446.448, 125.157, -542.060,  20.4894, -0.1502,  -0.2470,  -0.8421   } }, // epsg.io/1314
        { s_ellipsoids.Bessel1841,    { -582,    -105,     -414,      -8.3,     1.04,     0.35,    -3.08     } },
        { s_ellipsoids.Bessel1841,    {  148,    -507,     -685,       0,       0,        0,        0        } },
        { s_ellipsoids.WGS72,         {    0,       0,       -4.5,    -0.22,    0,        0,        0.554    } },
        { s_ellipsoids.WGS84,         {    0.0,     0.0,      0.0,     0.0,     0.0,      0.0,      0.0      } }
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
}

#endif //ELLIPSOIDS_H
