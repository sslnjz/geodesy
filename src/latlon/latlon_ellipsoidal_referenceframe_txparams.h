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
*  furnished to do so, subject to the following conditions,                       *
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
#ifndef LATLON_ELLIPSOIDAL_REFERENCEFRAME_TXPARAMS_H
#define LATLON_ELLIPSOIDAL_REFERENCEFRAME_TXPARAMS_H

#include <map>
#include <string>
#include <vector>

namespace geodesy
{
   /* Helmert transform parameters */
   struct Helmert
   {
      std::string epoch;
      double      params[7];
      double      rates[7];
   };

   struct HelmertTransforms
   {
      std::string from;
      std::string to;
      Helmert helmert;
   };

   static std::vector<HelmertTransforms> s_txParams = {
      {
         "ITRF2014",
         "ITRF2008",
         {
            "2010.0",
            /* tx(mm)   ty(mm)   tz(mm) s(ppb)  rx(mas)  ry(mas)  rz(mas) */
            {   1.6,    1.9,    2.4,  -0.02,   0.00,   0.00,   0.00   },
            {   0.0,    0.0,   -0.1,   0.03,   0.00,   0.00,   0.00   }
         }
      },
      {
         "ITRF2014",
         "ITRF2005",
         {
            "2010.0",
            {   2.6,    1.0,   -2.3,   0.92,   0.00,   0.00,   0.00   },
            {   0.3,    0.0,   -0.1,   0.03,   0.00,   0.00,   0.00   }
         }
      },
      {
         "ITRF2014",
         "ITRF2000",
         {
            "2010.0",
            {   0.7,    1.2,   -26.1,   2.12,   0.00,   0.00,   0.00   },
            {   0.1,    0.1,   -1.9,   0.11,   0.00,   0.00,   0.00   }
         }
      },
      {
         "ITRF2014",
         "ITRF97",
         {
            "2010.0",
            {   7.4,   -0.5,   -62.8,   3.80,   0.00,   0.00,   0.26   },
            {   0.1,   -0.5,   -3.3,   0.12,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2014",
         "ITRF96",
         {
            "2010.0",
            {   7.4,   -0.5,   -62.8,   3.80,   0.00,   0.00,   0.26   },
            {   0.1,   -0.5,   -3.3,   0.12,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2014",
         "ITRF94",
         {
            "2010.0",
            {   7.4,   -0.5,   -62.8,   3.80,   0.00,   0.00,   0.26   },
            {   0.1,   -0.5,   -3.3,   0.12,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2014",
         "ITRF93",
         {
            "2010.0",
            {  -50.4,    3.3,   -60.2,   4.29,  -2.81,   -3.38,   0.40   },
            {   -2.8,   -0.1,   -2.5,   0.12,  -0.11,   -0.19,   0.07   }
         }
      },
      {
         "ITRF2014",
         "ITRF92",
         {
            "2010.0",
            {   15.4,    1.5,   -70.8,   3.09,   0.00,   0.00,   0.26   },
            {   0.1,   -0.5,   -3.3,   0.12,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2014",
         "ITRF91",
         {
            "2010.0",
            {   27.4,   15.5,   -76.8,   4.49,   0.00,   0.00,   0.26   },
            {   0.1,   -0.5,   -3.3,   0.12,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2014",
         "ITRF90",
         {
            "2010.0",
            {   25.4,   11.5,   -92.8,   4.79,   0.00,   0.00,   0.26   },
            {   0.1,   -0.5,   -3.3,   0.12,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2014",
         "ITRF89",
         {
            "2010.0",
            {   30.4,   35.5,  -130.8,   8.19,   0.00,   0.00,   0.26   },
            {   0.1,   -0.5,   -3.3,   0.12,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2014",
         "ITRF88",
         {
            "2010.0",
            {   25.4,   -0.5,  -154.8,  11.29,   0.10,   0.00,   0.26   },
            {   0.1,   -0.5,   -3.3,   0.12,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2008",
         "ITRF2005",
         {
            "2000.0",
            {   -2.0,   -0.9,   -4.7,   0.94,   0.00,   0.00,   0.00   },
            {   0.3,    0.0,    0.0,   0.00,   0.00,   0.00,   0.00   }
         }
      },
      {
         "ITRF2008",
         "ITRF2000",
         {
            "2000.0",
            {   -1.9,   -1.7,   -10.5,   1.34,   0.00,   0.00,   0.00   },
            {   0.1,    0.1,   -1.8,   0.08,   0.00,   0.00,   0.00   }
         }
      },
      {
         "ITRF2008",
         "ITRF97",
         {
            "2000.0",
            {   4.8,    2.6,   -33.2,   2.92,   0.00,   0.00,   0.06   },
            {   0.1,   -0.5,   -3.2,   0.09,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2008",
         "ITRF96",
         {
            "2000.0",
            {   4.8,    2.6,   -33.2,   2.92,   0.00,   0.00,   0.06   },
            {   0.1,   -0.5,   -3.2,   0.09,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2008",
         "ITRF94",
         {
            "2000.0",
            {   4.8,    2.6,   -33.2,   2.92,   0.00,   0.00,   0.06   },
            {   0.1,   -0.5,   -3.2,   0.09,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2008",
         "ITRF93",
         {
            "2000.0",
            {  -24.0,    2.4,   -38.6,   3.41,  -1.71,   -1.48,   -0.30   },
            {   -2.8,   -0.1,   -2.4,   0.09,  -0.11,   -0.19,   0.07   }
         }
      },
      {
         "ITRF2008",
         "ITRF92",
         {
            "2000.0",
            {   12.8,    4.6,   -41.2,   2.21,   0.00,   0.00,   0.06   },
            {   0.1,   -0.5,   -3.2,   0.09,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2008",
         "ITRF91",
         {
            "2000.0",
            {   24.8,   18.6,   -47.2,   3.61,   0.00,   0.00,   0.06   },
            {   0.1,   -0.5,   -3.2,   0.09,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2008",
         "ITRF90",
         {
            "2000.0",
            {   22.8,   14.6,   -63.2,   3.91,   0.00,   0.00,   0.06   },
            {   0.1,   -0.5,   -3.2,   0.09,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2008",
         "ITRF89",
         {
            "2000.0",
            {   27.8,   38.6,  -101.2,   7.31,   0.00,   0.00,   0.06   },
            {   0.1,   -0.5,   -3.2,   0.09,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2008",
         "ITRF88",
         {
            "2000.0",
            {   22.8,    2.6,  -125.2,  10.41,   0.10,   0.00,   0.06   },
            {   0.1,   -0.5,   -3.2,   0.09,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2005",
         "ITRF2000",
         {
            "2000.0",
            {   0.1,   -0.8,   -5.8,   0.40,   0.000,   0.000,   0.000  },
            {   -0.2,    0.1,   -1.8,   0.08,   0.000,   0.000,   0.000  }
         }
      },
      {
         "ITRF2000",
         "ITRF97",
         {
            "1997.0",
            {   0.67,    0.61,   -1.85,  1.55,   0.00,   0.00,   0.00   },
            {   0.00,   -0.06,   -0.14,  0.01,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2000",
         "ITRF96",
         {
            "1997.0",
            {   0.67,    0.61,   -1.85,  1.55,   0.00,   0.00,   0.00   },
            {   0.00,   -0.06,   -0.14,  0.01,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2000",
         "ITRF94",
         {
            "1997.0",
            {   0.67,    0.61,   -1.85,  1.55,   0.00,   0.00,   0.00   },
            {   0.00,   -0.06,   -0.14,  0.01,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2000",
         "ITRF93",
         {
            "1988.0",
            {   12.7,    6.5,   -20.9,   1.95,  -0.39,   0.80,   -1.14   },
            {   -2.9,   -0.2,   -0.6,   0.01,  -0.11,   -0.19,   0.07   }
         }
      },
      {
         "ITRF2000",
         "ITRF92",
         {
            "1988.0",
            {   1.47,    1.35,   -1.39,  0.75,   0.00,   0.00,   -0.18   },
            {   0.00,   -0.06,   -0.14,  0.01,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2000",
         "ITRF91",
         {
            "1988.0",
            {   26.7,   27.5,   -19.9,   2.15,   0.00,   0.00,   -0.18   },
            {   0.0,   -0.6,   -1.4,   0.01,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2000",
         "ITRF90",
         {
            "1988.0",
            {   2.47,    2.35,   -3.59,  2.45,   0.00,   0.00,   -0.18   },
            {   0.00,   -0.06,   -0.14,  0.01,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2000",
         "ITRF89",
         {
            "1988.0",
            {   2.97,    4.75,   -7.39,  5.85,   0.00,   0.00,   -0.18   },
            {   0.00,   -0.06,   -0.14,  0.01,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2000",
         "ITRF88",
         {
            "1988.0",
            {   2.47,    1.15,   -9.79,  8.95,   0.10,   0.00,   -0.18   },
            {   0.00,   -0.06,   -0.14,  0.01,   0.00,   0.00,   0.02   }
         }
      },
      {
         "ITRF2000",
         "NAD83",
         {
            "1997.0", // note NAD83(CORS96)
            {  995.6, -1901.3,  -521.5,   0.62,  25.915,   9.426,  11.599  },
            {   0.7,   -0.7,    0.5,  -0.18,   0.067,  -0.757,  -0.051  }
         }
      },
      {
         "ITRF2014",
         "ETRF2000",
         {
            "2000.0",
            {   53.7,   51.2,   -55.1,   1.02,   0.891,   5.390,  -8.712  },
            {   0.1,    0.1,   -1.9,   0.11,   0.081,   0.490,  -0.792  }
         }
      },
      {
         "ITRF2008",
         "ETRF2000",
         {
            "2000.0",
            {   52.1,   49.3,   -58.5,   1.34,   0.891,   5.390,  -8.712  },
            {   0.1,    0.1,   -1.8,   0.08,   0.081,   0.490,  -0.792  }
         }
      },
      {
         "ITRF2005",
         "ETRF2000",
         {
            "2000.0",
            {   54.1,   50.2,   -53.8,   0.40,   0.891,   5.390,  -8.712  },
            {   -0.2,    0.1,   -1.8,   0.08,   0.081,   0.490,  -0.792  }
         }
      },
      {
         "ITRF2000",
         "ETRF2000",
         {
            "2000.0",
            {   54.0,   51.0,   -48.0,   0.00,   0.891,   5.390,  -8.712  },
            {   0.0,    0.0,    0.0,   0.00,   0.081,   0.490,  -0.792  }
         }
      },
      {
         "ITRF2008",
         "GDA94",
         {
            "1994.0",
            {  -84.68,  -19.42,   32.01,  9.710, -0.4254,  2.2578,  2.4015 },
            {   1.42,   1.34,   0.90,  0.109,  1.5461,  1.1820,  1.1551 }
         }
      },
      {
         "ITRF2005",
         "GDA94",
         {
            "1994.0",
            {  -79.73,   -6.86,   38.03,  6.636,  0.0351, -2.1211, -2.1411 },
            {   2.25,   -0.62,   -0.56,  0.294, -1.4707, -1.1443, -1.1701 }
         }
      },
      {
         "ITRF2000",
         "GDA94",
         {
            "1994.0",
            {  -45.91,  -29.85,  -20.37,  7.070, -1.6705,  0.4594,  1.9356 },
            {   -4.66,   3.55,   11.24,  0.249,  1.7454,  1.4868,  1.2240 }
         }
      }
   };

   /* Note WGS84(G730/G873/G1150) are coincident with ITRF at 10-centimetre level; WGS84(G1674) and
    * ITRF20014 / ITRF2008 ‘are likely to agree at the centimeter level’ (QPS).
    *
    * sources:
    * - ITRS: itrf.ensg.ign.fr/trans_para.php
    * - NAD83: Transforming Positions and Velocities between the International Terrestrial Reference
    *    Frame of 2000 and North American Datum of 1983, Soler & Snay, 2004;
    *    www.ngs.noaa.gov/CORS/Articles/SolerSnayASCE.pdf
    * - ETRS: etrs89.ensg.ign.fr/memo-V8.pdf / www.euref.eu/symposia/2016SanSebastian/01-02-Altamimi.pdf
    * - GDA:  ITRF to GDA94 coordinate transformations, Dawson & Woods, 2010
    *    (note sign of rotations for GDA94 reversed from Dawson & Woods 2010 as “Australia assumes rotation
    *    to be of coordinate axes” rather than the more conventional “position around the coordinate axes”)
    * more are available at:
    * confluence.qps.nl/qinsy/files/en/29856813/45482834/2/1453459502000/ITRF_Transformation_Parameters.xlsx
    */
}

#endif //LATLON_ELLIPSOIDAL_REFERENCEFRAME_TXPARAMS_H
