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

#include "osgridref.h"

#include "algorithm.h"
#include "latlon_osgridref.h"
#include "strutil.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <iomanip>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

constexpr double maxEastingMetres = 700000.0;
constexpr double maxNorthingMetres = 1300000.0;
constexpr double meridionalArcToleranceMetres = 0.00001;
constexpr int maxMeridionalArcIterations = 20;

[[nodiscard]] bool isFiniteInRange(double value, double min, double max)
{
    return std::isfinite(value) && value >= min && value <= max;
}

[[nodiscard]] std::string formatNumericMetres(double value)
{
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(3) << value;

    std::string result = stream.str();
    while (!result.empty() && result.back() == '0')
        result.pop_back();
    if (!result.empty() && result.back() == '.')
        result.pop_back();

    const std::size_t decimal = result.find('.');
    const std::size_t integerWidth = decimal == std::string::npos ? result.length() : decimal;
    if (integerWidth < 6)
        result.insert(0, 6 - integerWidth, '0');

    return result;
}

[[nodiscard]] int letterIndex(char letter)
{
    int index = std::toupper(static_cast<unsigned char>(letter)) - 'A';
    if (index > 7)
        --index; // The National Grid skips I.
    return index;
}

[[nodiscard]] int64_t truncatedGridDigits(double metres, int digits)
{
    const double scale = std::pow(10.0, 5 - digits / 2.0);
    return static_cast<int64_t>(std::floor(std::fmod(metres, 100000.0) / scale));
}

}

namespace geodesy
{

OsGridRef::OsGridRef(double easting, double northing)
   : m_easting(easting)
   , m_northing(northing)
{
    if (!isFiniteInRange(m_easting, 0.0, maxEastingMetres))
        throw std::invalid_argument("invalid easting");
    if (!isFiniteInRange(m_northing, 0.0, maxNorthingMetres))
        throw std::invalid_argument("invalid northing");
}

double OsGridRef::easting() const
{
    return m_easting;
}

double OsGridRef::northing() const
{
    return m_northing;
}

LatLonOsGridRef OsGridRef::toLatLon(Datum datum) const
{
    const double E = m_easting;
    const double N = m_northing;
    const auto [ a, b, f ] = nationalGrid.ellipsoid;
    const double phi0 = toRadians(nationalGrid.trueOrigin.lat);
    const double lambda0 = toRadians(nationalGrid.trueOrigin.lon);
    const double E0 = -nationalGrid.falseOrigin.easting;
    const double N0 = -nationalGrid.falseOrigin.northing;
    const double F0 = nationalGrid.scaleFactor;

    const double e2 = 1.0 - (b * b) / (a * a);
    const double n = (a - b) / (a + b);
    const double n2 = n * n;
    const double n3 = n2 * n;

    double phi = phi0;
    double meridionalArc = 0.0;
    int iterations = 0;
    do {
        phi = (N - N0 - meridionalArc) / (a * F0) + phi;

        // OSGB36 inverse projection iterates latitude until the meridional arc matches northing.
        const double Ma = (1.0 + n + (5.0 / 4.0) * n2 + (5.0 / 4.0) * n3) * (phi - phi0);
        const double Mb = (3.0 * n + 3.0 * n2 + (21.0 / 8.0) * n3)
            * std::sin(phi - phi0) * std::cos(phi + phi0);
        const double Mc = ((15.0 / 8.0) * n2 + (15.0 / 8.0) * n3)
            * std::sin(2.0 * (phi - phi0)) * std::cos(2.0 * (phi + phi0));
        const double Md = (35.0 / 24.0) * n3
            * std::sin(3.0 * (phi - phi0)) * std::cos(3.0 * (phi + phi0));
        meridionalArc = b * F0 * (Ma - Mb + Mc - Md);
        ++iterations;
    } while (std::abs(N - N0 - meridionalArc) >= meridionalArcToleranceMetres
        && iterations < maxMeridionalArcIterations);

    if (iterations == maxMeridionalArcIterations)
        throw std::runtime_error("OS grid inverse projection failed to converge");

    const double cosPhi = std::cos(phi);
    const double sinPhi = std::sin(phi);
    const double nu = a * F0 / std::sqrt(1.0 - e2 * sinPhi * sinPhi);
    const double rho = a * F0 * (1.0 - e2) / std::pow(1.0 - e2 * sinPhi * sinPhi, 1.5);
    const double eta2 = nu / rho - 1.0;

    const double tanPhi = std::tan(phi);
    const double tan2Phi = tanPhi * tanPhi;
    const double tan4Phi = tan2Phi * tan2Phi;
    const double tan6Phi = tan4Phi * tan2Phi;
    const double secPhi = 1.0 / cosPhi;
    const double nu3 = nu * nu * nu;
    const double nu5 = nu3 * nu * nu;
    const double nu7 = nu5 * nu * nu;
    const double VII = tanPhi / (2.0 * rho * nu);
    const double VIII = tanPhi / (24.0 * rho * nu3) * (5.0 + 3.0 * tan2Phi + eta2 - 9.0 * tan2Phi * eta2);
    const double IX = tanPhi / (720.0 * rho * nu5) * (61.0 + 90.0 * tan2Phi + 45.0 * tan4Phi);
    const double X = secPhi / nu;
    const double XI = secPhi / (6.0 * nu3) * (nu / rho + 2.0 * tan2Phi);
    const double XII = secPhi / (120.0 * nu5) * (5.0 + 28.0 * tan2Phi + 24.0 * tan4Phi);
    const double XIIA = secPhi / (5040.0 * nu7)
        * (61.0 + 662.0 * tan2Phi + 1320.0 * tan4Phi + 720.0 * tan6Phi);

    const double dE = E - E0;
    const double dE2 = dE * dE;
    const double dE3 = dE2 * dE;
    const double dE4 = dE2 * dE2;
    const double dE5 = dE3 * dE2;
    const double dE6 = dE4 * dE2;
    const double dE7 = dE5 * dE2;
    phi = phi - VII * dE2 + VIII * dE4 - IX * dE6;
    const double lambda = lambda0 + X * dE - XI * dE3 + XII * dE5 - XIIA * dE7;

    LatLonOsGridRef point(toDegrees(phi), toDegrees(lambda), 0.0, LatLonEllipsoidalDatum::datums().OSGB36);

    if (datum != LatLonEllipsoidalDatum::datums().OSGB36) {
        // Grid formulae are OSGB36; convert through the datum bridge only after the projection step.
        const LatLonOsGridRef converted = point.convertDatum(datum);
        point = LatLonOsGridRef(converted.lat(), converted.lon(), converted.height(), converted.datum());
    }

    return point;
}

OsGridRef OsGridRef::parse(const std::string& gridref)
{
    const std::string ref = strutil::strip(gridref);

    std::smatch numericMatch;
    if (std::regex_match(ref, numericMatch, std::regex(R"(^(\d+(?:\.\d+)?),\s*(\d+(?:\.\d+)?)$)"))) {
        return OsGridRef(std::stod(numericMatch[1].str()), std::stod(numericMatch[2].str()));
    }

    std::smatch gridMatch;
    if (!std::regex_match(ref, gridMatch,
            std::regex(R"(^([HNSThnst])([ABCDEFGHJKLMNOPQRSTUVWXYZabcdefghjklmnopqrstuvwxyz])\s*([0-9]+)\s*([0-9]*)$)"))) {
        throw std::invalid_argument("invalid grid reference");
    }

    const int l1 = letterIndex(ref[0]);
    const int l2 = letterIndex(ref[1]);

    const int e100km = ((l1 - 2) % 5) * 5 + (l2 % 5);
    const int n100km = (19 - static_cast<int>(std::floor(l1 / 5.0)) * 5)
        - static_cast<int>(std::floor(l2 / 5.0));

    std::vector<std::string> en = strutil::split_regex(strutil::strip(ref.substr(2)), "\\s+");
    if (en.size() == 1) {
        if (en[0].length() % 2 != 0)
            throw std::invalid_argument("invalid grid reference");
        en = { en[0].substr(0, en[0].length() / 2), en[0].substr(en[0].length() / 2) };
    }

    if (en.size() != 2 || en[0].empty() || en[0].length() != en[1].length() || en[0].length() > 5)
        throw std::invalid_argument("invalid grid reference");

    en[0] = strutil::padRight(en[0], 5, '0');
    en[1] = strutil::padRight(en[1], 5, '0');

    const double easting = e100km * 100000.0 + std::stod(en[0]);
    const double northing = n100km * 100000.0 + std::stod(en[1]);

    return OsGridRef(easting, northing);
}

std::string OsGridRef::toString(int digits) const
{
    constexpr std::array<int, 9> validPrecisions { 0, 2, 4, 6, 8, 10, 12, 14, 16 };
    if (std::find(validPrecisions.begin(), validPrecisions.end(), digits) == validPrecisions.end())
        throw std::invalid_argument("invalid precision");

    if (digits == 0)
        return formatNumericMetres(m_easting) + "," + formatNumericMetres(m_northing);

    const int e100km = static_cast<int>(std::floor(m_easting / 100000.0));
    const int n100km = static_cast<int>(std::floor(m_northing / 100000.0));

    int l1 = (19 - n100km) - (19 - n100km) % 5 + static_cast<int>(std::floor((e100km + 10) / 5.0));
    int l2 = (19 - n100km) * 5 % 25 + e100km % 5;
    if (l1 > 7)
        ++l1;
    if (l2 > 7)
        ++l2;

    std::string letterPair;
    letterPair += static_cast<char>(l1 + 'A');
    letterPair += static_cast<char>(l2 + 'A');

    const std::string easting = strutil::padLeft(std::to_string(truncatedGridDigits(m_easting, digits)),
        static_cast<std::size_t>(digits / 2), '0');
    const std::string northing = strutil::padLeft(std::to_string(truncatedGridDigits(m_northing, digits)),
        static_cast<std::size_t>(digits / 2), '0');

    return letterPair + " " + easting + " " + northing;
}

}
