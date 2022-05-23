#include <iostream>

#include "geodesy/dms.h"

int main()
{
   std::cout << geodesy::Dms::toDms(000.12, geodesy::Dms::DMS, 2) <<std::endl;
   return 0;
}