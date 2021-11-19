﻿/**********************************************************************************
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
#include "latlon_nvector_spherical.h"
#include "nvector_spherical.h"

using namespace geodesy;

LatLonNvectorSpherical::LatLonNvectorSpherical() : m_lat(0.0), m_lon(0.0)
{
}

LatLonNvectorSpherical::LatLonNvectorSpherical(double lat, double lon) : m_lat(lat), m_lon(lon)
{
}

NvectorSpherical LatLonNvectorSpherical::toNvector() const
{
   const auto φ = toRadians(m_lat);
   const auto λ = toRadians(m_lon);

   const auto sinφ = std::sin(φ), cosφ = std::cos(φ);
   const auto sinλ = std::sin(λ), cosλ = std::cos(λ);

   // right-handed vector: x -> 0°E,0°N; y -> 90°E,0°N, z -> 90°N
   const auto x = cosφ * cosλ;
   const auto y = cosφ * sinλ;
   const auto z = sinφ;

   return { x, y, z};
}

vector3d LatLonNvectorSpherical::greatCircle(double bearing) const
{
   const auto φ = toRadians(m_lat);
   const auto λ = toRadians(m_lon);
   const auto θ = toRadians(bearing);

   const auto x = std::sin(λ) * std::cos(θ) - std::sin(φ) * std::cos(λ) * std::sin(θ);
   const auto y = -std::cos(λ) * std::cos(θ) - std::sin(φ) * std::sin(λ) * std::sin(θ);
   const auto z = std::cos(φ) * std::sin(θ);

   return { x, y, z };
}

double LatLonNvectorSpherical::distanceTo(const LatLonNvectorSpherical& point, double radius) const
{
   const auto R = radius;

   const auto n1 = toNvector();
   const auto n2 = toNvector();

   const auto sinθ = n1.cross(n2).length();
   const auto cosθ = n1.dot(n2);
   const auto δ = std::atan2(sinθ, cosθ); // tanδ = |n₁×n₂| / n₁⋅n₂

   return δ * R;
}

double LatLonNvectorSpherical::initialBearingTo(const LatLonNvectorSpherical& point) const
{
   if (equals(point))
      return NAN; // coincident points

   const NvectorSpherical p1 = toNvector();
   const NvectorSpherical p2 = point.toNvector();

   const auto N = NvectorSpherical(0, 0, 1); // n-vector representing north pole

   const auto c1 = p1.cross(p2); // great circle through p1 & p2
   const auto c2 = p1.cross(N); // great circle through p1 & north pole

   const auto θ = c1.angleTo(c2, p1); // bearing is (signed) angle between c1 & c2

   return Dms::wrap360(toDegrees(θ)); // normalise to range 0..360°
}

double LatLonNvectorSpherical::finalBearingTo(const LatLonNvectorSpherical& point) const
{
   // get initial bearing from destination point to this point & reverse it by adding 180°
   return Dms::wrap360(point.initialBearingTo(*this) + 180);
}

LatLonNvectorSpherical LatLonNvectorSpherical::midpointTo(const LatLonNvectorSpherical& point) const
{
   const NvectorSpherical n1 = toNvector();
   const NvectorSpherical n2 = point.toNvector();

   const auto mid = n1 + n2;

   return NvectorSpherical(mid.x(), mid.y(), mid.z()).toLatLon();
}

LatLonNvectorSpherical LatLonNvectorSpherical::intermediatePointTo(const LatLonNvectorSpherical& point, double fraction) const
{
   // angular distance between points; tanδ = |n₁×n₂| / n₁⋅n₂
   const auto n1 = toNvector();
   const auto n2 = point.toNvector();
   const auto sinθ = n1.cross(n2).length();
   const auto cosθ = n1.dot(n2);
   const auto δ = std::atan2(sinθ, cosθ);

   // interpolated angular distance on straight line between points
   const auto δi = δ * fraction;
   const auto sinδi = std::sin(δi);
   const auto cosδi = std::cos(δi);

   // direction vector (perpendicular to n1 in plane of n2)
   const auto d = n1.cross(n2).unit().cross(n1); // unit(n₁×n₂) × n₁

   // interpolated position
   const auto ip = n1 * cosδi  + d * sinδi; // n₁⋅cosδᵢ + d⋅sinδᵢ

   return NvectorSpherical(ip.x(), ip.y(), ip.z()).toLatLon();
}

LatLonNvectorSpherical LatLonNvectorSpherical::intermediatePointOnChordTo(const LatLonNvectorSpherical& point, double fraction) const
{
   const auto n1 = toNvector();
   const auto n2 = point.toNvector();
   const auto inter = (n1 + (n2 -n1)) * fraction; // n₁ + (n₂−n₁)·f ≡ n₁·(1-f) + n₂·f

   return NvectorSpherical(inter.x(), inter.y(), inter.z()).toLatLon();
}

LatLonNvectorSpherical LatLonNvectorSpherical::destinationPoint(double distance, double bearing, double radius) const
{
   const auto n1 = toNvector();        // Gade's n_EA_E
   const auto δ = distance / radius;            // angular distance in radians
   const auto θ = toRadians(bearing);           // initial bearing in radians

   const auto N = NvectorSpherical(0, 0, 1);          // north pole

   const auto de = N.cross(n1).unit();               // east direction vector @ n1 (Gade's k_e_E)
   const auto dn = n1.cross(de);                     // north direction vector @ n1 (Gade's (k_n_E)

   const auto deSinθ = de * (std::sin(θ));
   const auto dnCosθ = dn * (std::cos(θ));

   const auto d = dnCosθ + deSinθ;               // direction vector @ n1 (≡ C×n1; C = great circle)

   const auto x = n1 * (std::cos(δ));             // component of n2 parallel to n1
   const auto y = d *  (std::sin(δ));              // component of n2 perpendicular to n1

   const auto n2 = x + y;                        // Gade's n_EB_E

   return NvectorSpherical(n2.x(), n2.y(), n2.z()).toLatLon();
}

LatLonNvectorSpherical LatLonNvectorSpherical::intersection(
   const LatLonNvectorSpherical& path1start,
   const LatLonNvectorSpherical& path1brngEnd,
   const LatLonNvectorSpherical& path2start,
   const LatLonNvectorSpherical& path2brngEnd)
{

   // if c1 & c2 are great circles through start and end points (or defined by start point + bearing),
   // then candidate intersections are simply c1 × c2 & c2 × c1; most of the work is deciding correct
   // intersection point to select! if bearing is given, that determines which intersection, if both
   // paths are defined by start/end points, take closer intersection

   const auto p1 = path1start.toNvector();
   const auto p2 = path2start.toNvector();

   const vector3d c1 = p1.cross(path1brngEnd.toNvector());
   const vector3d c2 = p2.cross(path2brngEnd.toNvector());

   // c1 & c2 are vectors defining great circles through start & end points; p × c gives initial bearing vector
   const auto i1 = c1.cross(c2);
   const auto i2 = c2.cross(c1);

   const auto mid = p1 + p2 + path1brngEnd.toNvector() + path2brngEnd.toNvector(); // eslint-disable-line no-case-declarations
   const auto intersec = mid.dot(i1) > 0 ? i1 : i2;
   return NvectorSpherical(intersec.x(), intersec.y(), intersec.z()).toLatLon();
}

LatLonNvectorSpherical LatLonNvectorSpherical::intersection(const LatLonNvectorSpherical& path1start,
   double path1brngEnd, const LatLonNvectorSpherical& path2start, double path2brngEnd)
{
   // if c1 & c2 are great circles through start and end points (or defined by start point + bearing),
   // then candidate intersections are simply c1 × c2 & c2 × c1; most of the work is deciding correct
   // intersection point to select! if bearing is given, that determines which intersection, if both
   // paths are defined by start/end points, take closer intersection

   const auto p1 = path1start.toNvector();
   const auto p2 = path2start.toNvector();

   const vector3d c1 = path1start.greatCircle(path1brngEnd);
   const vector3d c2 = path2start.greatCircle(path2brngEnd);

   // c1 & c2 are vectors defining great circles through start & end points; p × c gives initial bearing vector
   const auto i1 = c1.cross(c2);
   const auto i2 = c2.cross(c1);

   // if c×p⋅i1 is +ve, the initial bearing is towards i1, otherwise towards antipodal i2
   const auto dir1 = sign(c1.cross(p1).dot(i1)); // c1×p1⋅i1 +ve means p1 bearing points to i1
   const auto dir2 = sign(c2.cross(p2).dot(i1)); // c2×p2⋅i1 +ve means p2 bearing points to i1

   vector3d intersec{};
   switch (dir1 + dir2)
   {
   case  2: // dir1, dir2 both +ve, 1 & 2 both pointing to i1
      intersec = i1;
      break;
   case -2: // dir1, dir2 both -ve, 1 & 2 both pointing to i2
      intersec = i2;
      break;
   case  0: // dir1, dir2 opposite; intersection is at further-away intersection point
       // take opposite intersection from mid-point of p1 & p2 [is this always true?]
      intersec = (p1 + p2).dot(i1) > 0 ? i2 : i1;
      break;
   default: 
      break;
   }

   return NvectorSpherical(intersec.x(), intersec.y(), intersec.z()).toLatLon();
}

double LatLonNvectorSpherical::crossTrackDistanceTo(const LatLonNvectorSpherical& pathStart,
                                                    const LatLonNvectorSpherical& pathBrngEnd, double radius) const
{
   if (equals(pathStart)) return 0;

   const auto p = toNvector();
   const auto gc = pathStart.toNvector().cross(pathBrngEnd.toNvector()); // great circle defined by two points
   const auto α = gc.angleTo(p) - π / 2; // angle between point & great-circle

   return α * radius;
}

double LatLonNvectorSpherical::crossTrackDistanceTo(const LatLonNvectorSpherical& pathStart, double pathBrngEnd,
   double radius) const
{
   if (equals(pathStart)) return 0;

   const auto p = toNvector();
   const auto gc = pathStart.greatCircle(pathBrngEnd);  // great circle defined by two points
   const auto α = gc.angleTo(p) - π / 2; // angle between point & great-circle

   return α * radius;
}

double LatLonNvectorSpherical::alongTrackDistanceTo(const LatLonNvectorSpherical& pathStart,
                                                    const LatLonNvectorSpherical& pathBrngEnd, double radius) const
{
   const auto p = toNvector();
   const auto gc = pathStart.toNvector().cross(pathBrngEnd.toNvector()); // great circle defined by two points
   const auto pat = gc.cross(p).cross(gc); // along-track point c × p × c
   const auto α = pathStart.toNvector().angleTo(pat, gc); // angle between start point and along-track point

   return α * radius;
}

double LatLonNvectorSpherical::alongTrackDistanceTo(const LatLonNvectorSpherical& pathStart, double pathBrngEnd,
   double radius) const
{
   const auto p = toNvector();
   const auto gc = pathStart.greatCircle(pathBrngEnd); // great circle defined by two points
   const auto pat = gc.cross(p).cross(gc); // along-track point c × p × c
   const auto α = pathStart.toNvector().angleTo(pat, gc); // angle between start point and along-track point

   return α * radius;
}

LatLonNvectorSpherical LatLonNvectorSpherical::nearestPointOnSegment(const LatLonNvectorSpherical& point1,
                                                                     const LatLonNvectorSpherical& point2) const
{
   if (isWithinExtent(point1, point2) && !point1.equals(point2))
   {
      // closer to segment than to its endpoints, find closest point on segment
      const auto n0 = toNvector(), n1 = point1.toNvector(), n2 = point2.toNvector();
      const auto c1 = n1.cross(n2); // n1×n2 = vector representing great circle through p1, p2
      const auto c2 = n0.cross(c1); // n0×c1 = vector representing great circle through p0 normal to c1
      const auto n = c1.cross(c2); // c2×c1 = nearest point on c1 to n0
      return NvectorSpherical(n.x(), n.y(), n.z()).toLatLon();
   }
   else
   {
      // beyond segment extent, take closer endpoint
      const auto d1 = distanceTo(point1);
      const auto d2 = distanceTo(point2);
      const auto pCloser = d1 < d2 ? point1 : point2;
      return LatLonNvectorSpherical(pCloser.lat(), pCloser.lon());
   }
}

LatLonNvectorSpherical LatLonNvectorSpherical::triangulate(const LatLonNvectorSpherical& point1, double bearing1,
                                                           const LatLonNvectorSpherical& point2, double bearing2)
{
   const auto n1 = point1.toNvector();
   const auto n2 = point2.toNvector();

   const auto θ1 = toRadians(bearing1);
   const auto θ2 = toRadians(bearing2);

   const auto N = NvectorSpherical(0, 0, 1); // north pole

   const auto de1 = N.cross(n1).unit(); // east vector @ n1
   const auto dn1 = n1.cross(de1); // north vector @ n1
   const auto de1Sinθ = de1 * (std::sin(θ1));
   const auto dn1Cosθ = dn1 * (std::cos(θ1));
   const auto d1 = dn1Cosθ + de1Sinθ; // direction vector @ n1

   const auto c1 = n1.cross(d1); // great circle p1 + bearing1

   const auto de2 = N.cross(n2).unit(); // east vector @ n2
   const auto dn2 = n2.cross(de2); // north vector @ n2
   const auto de2Sinθ = de2 * sin(θ2);
   const auto dn2Cosθ = dn2 * cos(θ2);
   const auto d2 = dn2Cosθ + de2Sinθ; // direction vector @ n2

   const auto c2 = n2.cross(d2); // great circle p2 + bearing2
   const auto ni = c1.cross(c2); // n-vector of intersection point

   return NvectorSpherical(ni.x(), ni.y(), ni.z()).toLatLon();
}

LatLonNvectorSpherical LatLonNvectorSpherical::trilaterate(const LatLonNvectorSpherical& point1, double distance1,
                                                           const LatLonNvectorSpherical& point2, double distance2,
                                                           const LatLonNvectorSpherical& point3, double distance3,
                                                           double radius)
{
   // from en.wikipedia.org/wiki/Trilateration

   const auto n1 = point1.toNvector();
   const auto n2 = point2.toNvector();
   const auto n3 = point3.toNvector();
   const auto δ1 = distance1 / radius;
   const auto δ2 = distance2 / radius;
   const auto δ3 = distance3 / radius;

   // the following uses x,y coordinate system with origin at n1, x axis n1->n2
   const auto eX = (n2 - n1).unit(); // unit vector in x direction n1->n2
   const auto i = eX.dot(n3 - n1); // signed magnitude of x component of n1->n3
   const auto eY = ((n3 - n1) - (eX * i )).unit(); // unit vector in y direction
   const auto d = (n2 - n1).length(); // distance n1->n2
   const auto j = eY.dot(n3 - n1); // signed magnitude of y component of n1->n3
   const auto x = (δ1 * δ1 - δ2 * δ2 + d * d) / (2 * d); // x component of n1 -> intersection
   const auto y = (δ1 * δ1 - δ3 * δ3 + i * i + j * j) / (2 * j) - x * i / j; // y component of n1 -> intersection
   // const eZ = eX.cross(eY);                            // unit vector perpendicular to plane
   // const z = Math.sqrt(δ1*δ1 - x*x - y*y);             // z will be NaN for no intersections

   if (!std::isfinite(x) || !std::isfinite(y))
      return {}; // coincident points?

   const auto n = (n1 + (eX * x)) + (eY * y); // note don't use z component; assume points at same height

   return NvectorSpherical(n.x(), n.y(), n.z()).toLatLon();
}

bool LatLonNvectorSpherical::isEnclosedBy(std::vector<LatLonNvectorSpherical>& polygon) const
{
   // this method uses angle summation test; on a plane, angles for an enclosed point will sum
   // to 360°, angles for an exterior point will sum to 0°. On a sphere, enclosed point angles
   // will sum to less than 360° (due to spherical excess), exterior point angles will be small
   // but non-zero. TODO: are any winding number optimisations applicable to spherical surface?

   // close the polygon so that the last point equals the first point
   const auto closed = polygon[0].equals(polygon[polygon.size() - 1]);
   if (!closed) polygon.push_back(polygon[0]);

   const size_t nVertices = polygon.size() - 1;

   const auto p = toNvector();

   // get vectors from p to each vertex
   std::vector<vector3d> vectorToVertex;
   for (size_t v = 0; v < nVertices; ++v) 
   {
      vectorToVertex[v] = (p - polygon[v].toNvector());
   }
   vectorToVertex.push_back(vectorToVertex[0]);

   // sum subtended angles of each edge (using vector p to determine sign)
   double Σθ = 0.0;
   for (size_t v = 0; v < nVertices; ++v) 
   {
      Σθ += vectorToVertex[v].angleTo(vectorToVertex[v + 1], p);
   }

   if (!closed)
   {
      polygon.pop_back(); // restore polygon to pristine condition
   }

   return std::abs(Σθ) > π;
}

bool LatLonNvectorSpherical::equals(const LatLonNvectorSpherical& point) const
{
   if (std::abs(m_lat - point.lat()) > std::numeric_limits<double>::epsilon()) return false;
   if (std::abs(m_lon - point.lon()) > std::numeric_limits<double>::epsilon()) return false;

   return true;
}

bool LatLonNvectorSpherical::isWithinExtent(const LatLonNvectorSpherical& point1,
                                            const LatLonNvectorSpherical& point2) const
{
   if (point1.equals(point2))
      return equals(point1); // null segment

   const auto n0 = toNvector();
   const auto n1 = point1.toNvector();
   const auto n2 = point2.toNvector(); // n-vectors

   // get vectors representing p0->p1, p0->p2, p1->p2, p2->p1
   const auto δ10 = n0 - n1, δ12 = n2 - n1;
   const auto δ20 = n0 - n2, δ21 = n1 - n2;

   // dot product δ10⋅δ12 tells us if p0 is on p2 side of p1, similarly for δ20⋅δ21
   const auto extent1 = δ10.dot(δ12);
   const auto extent2 = δ20.dot(δ21);

   const auto isSameHemisphere = n0.dot(n1) >= 0 && n0.dot(n2) >= 0;

   return extent1 >= 0 && extent2 >= 0 && isSameHemisphere;
}

double LatLonNvectorSpherical::areaOf(std::vector<LatLonNvectorSpherical> &polygon, double radius)
{
    // close the polygon so that the last point equals the first point
    const auto closed = polygon[0].equals(polygon[polygon.size()-1]);
    if (!closed) polygon.push_back(polygon[0]);

    const auto n = polygon.size() - 1; // number of vertices

    // get great-circle vector for each segment
    std::vector<vector3d> c{};
    for (size_t v=0; v < n; ++v)
    {
        const auto i = polygon[v].toNvector();
        const auto j = polygon[v+1].toNvector();
        c[v] = i.cross(j); // great circle for segment v..v+1
    }
    c.push_back(c[0]);

    // sum interior angles; depending on whether polygon is cw or ccw, angle between edges is
    // π−α or π+α, where α is angle between great-circle vectors; so sum α, then take n·π − |Σα|
    // (cannot use Σ(π−|α|) as concave polygons would fail); use vector to 1st point as plane
    // normal for sign of α
    const auto n1 = polygon[0].toNvector();
    double Σα = 0.0;
    for (size_t v=0; v<n; v++) {
        Σα += c[v].angleTo(c[v+1], n1);
    }
    const auto Σθ = n*π - std::abs(Σα);

    // note: angle between two sides of a spherical triangle is acos(c₁·c₂) where cₙ is the
    // plane normal vector to the great circle representing the triangle side - use this instead
    // of angleTo()?

    const auto E = (Σθ - (n-2)*π); // spherical excess (in steradians)
    const auto A = E * radius*radius;        // area in units of R²

    if (!closed)
        polygon.pop_back(); // restore polygon to pristine condition

    return A;
}

LatLonNvectorSpherical LatLonNvectorSpherical::meanOf(const std::vector<LatLonNvectorSpherical> &points)
{
    auto m = NvectorSpherical(0, 0, 0); // null vector

    // add all vectors
    for (size_t p = 0; p < points.size(); p++)
    {
        m += points[p].toNvector();
    }
    // m is now geographic mean
    return NvectorSpherical(m.x(), m.y(), m.z()).toLatLon();
}

std::string LatLonNvectorSpherical::toString(Dms::eFormat format)
{
    std::stringstream ss;
    if (format == Dms::N)
    { // signed numeric degrees
        ss << std::fixed << std::setprecision(4) << m_lat << "," << m_lon;
        return ss.str();
    }
    const auto lat = Dms::toLatitude(m_lat, format);
    const auto lon = Dms::toLongitude(m_lon, format);

    ss << Dms::toLatitude(m_lat, format) << ", " << Dms::toLongitude(m_lon, format);
    return ss.str();
}
