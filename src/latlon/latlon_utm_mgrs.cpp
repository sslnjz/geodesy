#include "latlon_utm_mgrs.h"
#include "utm.h"

using geodesy::LatLonUtmMgrs;

LatLonUtmMgrs::LatLonUtmMgrs(double lat, double lon, double height, std::optional<Datum> datum, std::optional<ReferenceFrame> reference, std::optional<std::string> epoch)
	: LatLonUtm(lat, lon, height, datum, reference, epoch)
{
}

geodesy::UtmMgrs geodesy::LatLonUtmMgrs::toUtm(std::optional<int> zoneOverride)
{
	const auto utm = LatLonUtm::toUtm(zoneOverride);
	return UtmMgrs(utm.m_zone, utm.m_hemisphere,utm.easting(), utm.northing(), utm.m_datum);
}
