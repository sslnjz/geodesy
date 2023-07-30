#include <iostream>

#include "geodesy/dms.h"
#include "geodesy/latlon_utm_mgrs.h"

int main()
{
   std::cout << geodesy::Dms::toDms(000.12, geodesy::Dms::DMS, 2) <<std::endl;


   geodesy::LatLonUtmMgrs latlonutmmgrs(0.012, 0.014);

   auto utm = latlonutmmgrs.toUtm();

   std::cout << utm.northing() << " ++++++" << std::endl;

   return 0;
}