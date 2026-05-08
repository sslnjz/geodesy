# Geodesy C++ Initial Port Analysis

Date: 2026-05-08

This document preserves the local static analysis of `D:\workspace\github\geodesy` for follow-up sessions. The local repository is a C++ initial port/implementation of Chris Veness' JavaScript geodesy library, not the original JS package layout.

## Source Facts

- Local repository remote: `git@github.com:sslnjz/geodesy.git`.
- Local `README.md`: "Libraries of geodesy functions implemented in c++".
- Local source layout is CMake/C++: `src/`, `include/`, `test/`, `ui/`.
- The original JS project is `chrisveness/geodesy`, with flat ES module files such as `dms.js`, `vector3d.js`, `latlon-spherical.js`, `latlon-ellipsoidal*.js`, `utm.js`, `mgrs.js`, and `osgridref.js`.

## Directory And Module Structure

| Module / Directory | Main Files | Purpose |
|---|---|---|
| `src/dms` | `dms.h`, `dms.cpp` | Degree/minute/second parsing, formatting, compass points, angle wrapping. |
| `src/utils` | `algorithm.h`, `strutil.h/.cpp` | Math constants/helpers, angle conversion, string splitting/padding/trim. |
| `src/latlon` | `latlon*.h/.cpp`, `ellipsoids.h`, `latlon_ellipsoidal_referenceframe_txparams.h` | Base lat/lon, spherical, ellipsoidal, datum, reference-frame, Vincenty, n-vector, UTM/OSGrid lat/lon adapters. |
| `src/vector` | `vector3d`, `cartesian`, `cartesian_datum`, `cartesian_referenceFrame`, `nvector_*`, `ned` | 3D vector, ECEF cartesian, datum/reference-frame cartesian transforms, n-vector and NED types. |
| `src/utm` | `utm`, `utm_mgrs`, `mgrs`, `osgridref` | UTM, MGRS, and UK OS grid references. |
| `include/geodesy` | generated/copied headers | Public headers copied from `src` by the build. |
| `test` | `*_unittest.cpp` | GTest tests for DMS, LatLon, spherical, ellipsoidal, n-vector, datum/reference-frame, Vincenty, vector3d. UTM/MGRS/OSGrid are not currently in `TEST_LIST`. |
| `ui` | `main.cpp` | Small local executable entrypoint. |

## Implemented Feature List

| Module | Implemented Features | Current Assessment |
|---|---|---|
| `dms` | `Dms::parse`, `toDms`, `toLat`, `toLon`, `toBearing`, `compassPoint`, `wrap90/180/360`, separator setter/getter. | Basic port exists, but global separator design is risky. |
| `latlon` | Base `LatLon`, normalized numeric/string construction, comma-separated string parse, `toString`, `toGeoJSON`, equality. | Basic usable subset; JS object/GeoJSON parse is missing. |
| `latlon_spherical` | Distance, bearings, midpoint, interpolation, destination, intersections, max latitude, crossing parallels, rhumb line, cross/along-track, polygon area. | Broad feature coverage; boundary/TODO risks remain. |
| `latlon_ellipsoidal` | Lat/lon/height, datum/reference-frame fields, ECEF conversion, parse overloads, formatting, equality. | Basic path exists; optional handling is unsafe. |
| `latlon_ellipsoidal_datum` | Datum-aware point, datum conversion, datum cartesian conversion. | Main chain exists, but datum identity and error semantics need hardening. |
| `latlon_ellipsoidal_referenceframe` | Reference frame point, epoch, transform parameter access, conversion entrypoint. | High risk; contains P0 defects. |
| `latlon_ellipsoidal_vincenty` | Vincenty inverse/direct, distance, initial/final bearing, destination, final bearing on path, intermediate point. | Algorithm port exists; convergence and exception semantics need test lock. |
| `latlon_nvector_spherical` | Spherical n-vector conversion, bearings, intersections, cross/along-track, nearest point, triangulate, trilaterate, polygon enclosure, area, mean. | Feature-rich but boundary behavior is not fully settled. |
| `latlon_nvector_ellipsoidal` | Ellipsoidal n-vector, NED delta, destination from delta, cartesian conversion. | Narrow but useful; needs datum/height verification. |
| `vector` | `vector3d`, `Cartesian`, `CartesianDatum`, `CartesianReferenceFrame`, `NvectorSpherical`, `NvectorEllipsoidal`, `NvectorCartesian`, `Ned`. | Core infrastructure exists; `vector3d` and reference-frame cartesian contain hard bugs. |
| `utm/mgrs/osgridref` | UTM parse/convert/format, MGRS parse/convert/format, OSGrid parse/convert/format. | Present in source, but UTM has obvious logic bugs and test coverage gap. |
| Build/test | CMake static library, generated include headers, GTest test targets. | Build framework exists; test dependency and packaging are fragile. |

## Core Classes And Dependencies

| Class / Type | Core Functions | Depends On |
|---|---|---|
| `Dms` | `parse`, `toDms`, `toLat`, `toLon`, `toBearing`, `compassPoint`, wrapping functions. | `strutil`, regex, global separator. |
| `vector3d` | vector arithmetic, dot/cross, unit, angle, rotation, formatting. | `algorithm.h` math helpers. |
| `LatLon` | parse, getters/setters, string/GeoJSON output, equality. | `Dms`, `strutil`. |
| `LatLonSpherical` | spherical trig geodesy and rhumb-line operations. | `LatLon`, `Dms`, `vector3d`, `algorithm`. |
| `LatLonEllipsoidal` | ECEF conversion, datum/reference-frame-aware storage. | `LatLon`, `Ellipsoid`, `Datum`, `ReferenceFrame`, `Cartesian`. |
| `LatLonEllipsoidalDatum` | datum conversion. | `LatLonEllipsoidal`, `CartesianDatum`, `g_datums`. |
| `LatLonEllipsoidalReferenceFrame` | reference-frame conversion. | `LatLonEllipsoidal`, `CartesianReferenceFrame`, `s_txParams`. |
| `LatLonEllipsoidalVincenty` | Vincenty geodesics. | `LatLonEllipsoidal`, datum ellipsoid, math helpers. |
| `LatLonNvectorSpherical` | spherical n-vector operations. | `LatLon`, `NvectorSpherical`, `vector3d`. |
| `LatLonNvectorEllipsoidal` | ellipsoidal n-vector and NED delta. | `LatLonEllipsoidal`, `Cartesian`, `NvectorCartesian`, `Ned`. |
| `Cartesian` | ECEF to geodetic. | `vector3d`, `Ellipsoid`, `LatLonEllipsoidal`. |
| `CartesianDatum` | Helmert 7-parameter datum transform. | `Cartesian`, `Datum`, `LatLonEllipsoidalDatum`. |
| `CartesianReferenceFrame` | Helmert 14-parameter reference-frame transform. | `Cartesian`, `ReferenceFrame`, `s_txParams`, epoch. |
| `Utm` / `LatLonUtm` | UTM parse and forward/reverse projection. | `Datum`, `LatLonEllipsoidal`, Karney/Kruger series. |
| `Mgrs` / `UtmMgrs` | MGRS parse/format and UTM conversion. | `Utm`, MGRS letter tables. |
| `OsGridRef` / `LatLonOsGridRef` | UK National Grid conversion. | `Datum`, OSGB36/WGS84 conversion, national grid constants. |

## Data Types, Constants, And Global State

| Item | Description | Risk |
|---|---|---|
| `Ellipsoid` | `a`, `b`, `f`. | Precision-sensitive but appropriate as `double`. |
| `Transform` | Helmert 7 parameters: translation, scale, rotations. | Equality uses very tight epsilon; datum identity may be brittle. |
| `Datum` | `Ellipsoid + Transform`. | Value type is acceptable; optional usage is unsafe in places. |
| `ReferenceFrame` | `name`, optional `epoch`, `ellipsoid`. | Epoch stored as string; numeric epoch math uses `stoi`, losing decimals. |
| `HelmertTransforms` | `from`, `to`, `Helmert{epoch, params[7], rates[7]}`. | Header-level mutable vector and raw arrays are not ideal C++17 API. |
| `Dms::_separator` | Static mutable separator string. | Global mutable state, non-thread-safe, intentionally leaked with `*new`. |
| `g_ellipsoids`, `g_datums`, `g_reference_frames` | Header-defined global constant tables via `*new` references. | Multiple TU copies, leaks, static initialization/lifetime smell. |
| `latBands`, `e100kLetters`, `n100kLetters` | MGRS letter tables. | Reasonable constants, but should be `inline constexpr std::string_view`/arrays. |
| `nationalGrid` | OS grid constants via nested static refs and `*new`. | Same global design smell as datum constants. |

## Asynchronous And Event Logic

The project is a pure calculation library. No async/event loop, callbacks, Promises, threads, I/O watchers, or observer/event bus patterns were found in the core library. Module interaction is synchronous value-in/value-out calls:

- parse/construct value objects;
- convert between coordinate representations;
- return new objects for transformed coordinates;
- throw exceptions for invalid input or unsupported transforms.

## Differences From Original JS Library

| Difference | Impact |
|---|---|
| Original JS uses flat ES module files; C++ uses CMake object libraries plus generated public headers. | Build artifacts and source headers can drift. |
| JS `parse` accepts number/string/object/GeoJSON shapes; C++ only implements numeric, two-string, and comma-string subsets. | API compatibility is incomplete. |
| JS can prototype-mixin methods between LatLon variants; C++ uses fixed classes/inheritance. | Dynamic extension behavior must be redesigned explicitly. |
| JS module constants are singletons; C++ constants are header-defined `static` refs to heap objects. | ABI/ODR/lifetime and test isolation risks. |
| JS exceptions are `TypeError`, `RangeError`, `EvalError`; C++ uses mixed standard exception types. | Consumers cannot rely on a stable error taxonomy. |
| JS can dynamically add fields like convergence/scale; C++ uses dedicated setters/fields. | Better typed, but API must be intentionally exposed. |
| JS numerical values use `Number`; C++ uses `double`. | Appropriate, but tolerance policy must be consistent. |

## Potential Migration Difficulties

| Topic | Detail | Recommended Direction |
|---|---|---|
| Dynamic JS input types | JS accepts numbers, strings, objects, GeoJSON. | Use overloads plus a separate JSON/GeoJSON adapter, not implicit ad-hoc parsing everywhere. |
| Prototype mixins | JS README explicitly documents runtime prototype copying. | Replace with explicit composition, policy types, or named adapter classes. |
| Floating precision | Geodesy algorithms are precision-sensitive. | Standardize tolerance helpers and test all known JS examples. |
| Vincenty convergence | Inverse/direct may fail for hard cases. | Use explicit `std::optional`/exception policy and regression tests. |
| Units | Degrees/radians, metres, mm, ppm, ppb, mas, decimal epoch. | Make unit conversions explicit and avoid string epoch arithmetic. |
| Unicode DMS | Degree, prime, double-prime, narrow spaces, NSEW suffixes. | Keep encoding tests and locale boundaries explicit. |
| Global state | `Dms.separator` and constant tables are global. | Prefer immutable constants and parameterized formatting state. |
| Header/API hygiene | Generated headers mirror source headers. | Treat `include/geodesy` as API and validate with install-tree tests. |

## Bugs And Design Defects

| Severity | Location | Problem | Must Change |
|---|---|---|---|
| P0 | `src/latlon/latlon_ellipsoidal_referenceframe.cpp:79` | `throw new std::invalid_argument(...)` throws a pointer and leaks. | `LatLonEllipsoidalReferenceFrame::convertReferenceFrame`. |
| P0 | `src/vector/cartesian_referenceFrame.cpp:230` | `applyTransform()` returns `CartesianReferenceFrame(x2,y2,z2)` with no epoch/reference frame; constructor rejects empty epoch. | `CartesianReferenceFrame::applyTransform` and/or constructor semantics. |
| P0 | `src/vector/cartesian_referenceFrame.cpp:167` | Empty `std::optional` is dereferenced/assigned via `intermediateTxFwd->first`. | `CartesianReferenceFrame::convertReferenceFrame`. |
| P0 | `src/utm/utm.cpp:61` and `:68` | UTM northing range checks are inverted, throwing for valid ranges. | `Utm::Utm`. |
| P0 | `src/utm/utm.cpp:87` | Regex parses one non-space token but later expects 4 capture groups. | `Utm::parse`. |
| P0 | `src/utm/utm.cpp:139` and `src/latlon/latlon_utm.cpp:75` | Karney/Kruger series loops start at `j=0`; formulas require 1-based terms or corrected indexing. | `Utm::toLatLon`, `LatLonUtm::toUtm`. |
| P0 | `src/vector/vector3d.h:366` | Unary minus returns the original vector instead of negating. | `operator-(const vector3d&)`. |
| P0 | `src/vector/vector3d.h:70` and related getters | Extra qualification inside class definition is non-standard and compiler-sensitive. | `vector3d` header. |
| P1 | `src/latlon/latlon_ellipsoidal.cpp:74` | `datum()` dereferences empty `m_datum`. | `LatLonEllipsoidal::datum`. |
| P1 | `src/dms/dms.cpp:40` | Separator uses leaked heap global mutable state. | `Dms` separator implementation. |
| P1 | `src/latlon/ellipsoids.h:83`, `:173`, `:229` | Header-defined global `*new` constant tables. | Ellipsoid/datum/reference-frame constants. |
| P1 | `src/latlon/latlon.cpp:57` | String constructor parses twice and supports fewer JS input forms. | `LatLon` string parse path. |
| P1 | `src/latlon/latlon_spherical.cpp:370` | Area/pole handling has TODO and platform-specific test notes. | `LatLonSpherical::areaOf`. |
| P2 | `src/CMakeLists.txt` | Non-Apple path creates `libgeodesy.so` linker script pointing at static archive. | Packaging/install layout. |
| P2 | `test/CMakeLists.txt` | Tests fetch googletest from network and omit UTM/MGRS/OSGrid targets. | Test infrastructure. |

## C++17 Type And Class Design Evaluation

| Area | Assessment |
|---|---|
| Value semantics | Mostly appropriate: coordinate objects are small value types. |
| `std::optional` | Overused and sometimes unsafe. Empty optional behavior is not consistently modeled. |
| Constants | Not idiomatic C++17. Prefer `inline constexpr`, `inline const`, `std::array`, `std::string_view`, or function-local immutable statics. |
| Exceptions | Needs a consistent taxonomy. Never throw pointers. |
| Headers | Public headers should compile cleanly standalone and avoid generated-copy drift. |
| Naming | Mixed class naming (`vector3d` vs `LatLon...`) should be cleaned only if API compatibility allows. |
| Floating comparison | Tolerances are inconsistent: `epsilon`, `1e-9`, and custom helpers are mixed. |
| Raw arrays | `double params[7]` / `rates[7]` should become `std::array<double, 7>` for safer C++17 APIs. |
| Epoch | String epoch plus `stoi` is wrong for decimal years such as `2010.0`; use `double` or a typed epoch wrapper. |

## Coupling And Maintainability

| Module | Coupling / Maintainability | Migration Risk |
|---|---|---|
| `utils/algorithm` | Simple, shared math helpers. | Low |
| `utils/strutil` | Parse behavior affects DMS/UTM/MGRS/OSGrid. | Medium |
| `dms` | Foundational formatting/parser with mutable global separator. | Medium |
| `vector3d` | Core dependency of Cartesian/n-vector/geometric algorithms; current hard bugs make blast radius large. | High |
| `latlon` | Basic coordinate base; parse compatibility incomplete. | Medium |
| `latlon_spherical` | Mostly independent algorithm module; edge cases complex. | Medium |
| `latlon_ellipsoidal` | Shared base for datum/reference-frame/Vincenty/UTM. | High |
| `datum/cartesian_datum` | Clear conversion chain but depends on brittle constants/equality. | High |
| `referenceFrame/cartesian_referenceFrame` | Highest defect density; transform chaining and epoch handling are fragile. | High |
| `vincenty` | Algorithmically complex but relatively isolated. | Medium |
| `nvector_spherical` | Many geometric features and edge cases. | High |
| `nvector_ellipsoidal/NED` | Narrower surface but depends on datum/cartesian correctness. | Medium |
| `utm` | Current parse/range/series defects make it high priority. | High |
| `mgrs` | Depends on UTM correctness and precision/truncation rules. | High |
| `osgridref` | Smaller domain but datum/grid constants and boundaries need tests. | Medium |
| `CMake/install` | Works as a first pass, but generated headers and pseudo `.so` are fragile. | Medium |

## Must-Change Priority List

1. Fix `src/vector/vector3d.h`: remove non-standard extra qualification and correct unary minus.
2. Fix `src/utm/utm.cpp` and `src/latlon/latlon_utm.cpp`: UTM parse, northing validation, and series indexing.
3. Fix `src/vector/cartesian_referenceFrame.cpp` and `src/latlon/latlon_ellipsoidal_referenceframe.cpp`: optional misuse, default epoch semantics, transform chaining, and pointer exception.
4. Fix `src/latlon/latlon_ellipsoidal.cpp`: safe datum access/default behavior.
5. Replace `*new` global constants in `src/dms/dms.cpp`, `src/latlon/ellipsoids.h`, and `src/latlon/latlon_osgridref.h`.
6. Expand tests to include UTM, MGRS, OSGrid, install-tree header compilation, and JS original behavior examples.

## Suggested Next Verification Plan

| Step | Goal |
|---|---|
| 1 | Build with MSVC and one strict compiler mode if available. |
| 2 | Run existing GTest suite and record failures. |
| 3 | Add/enable UTM, MGRS, OSGrid tests before changing projection logic. |
| 4 | Port JS original examples as golden tests module by module. |
| 5 | Fix P0 defects first, then P1 design debt. |
| 6 | Re-run tests from source tree and installed package headers. |

## Bottom Line

The C++ initial port already covers much of the JS library's functional surface, but it should not be treated as a stable geodesy library yet. The highest risks are not missing files; they are hard correctness defects in `vector3d`, UTM, and reference-frame conversion, plus non-idiomatic global constants and unsafe optional handling. The migration should proceed by locking JS behavior with tests and fixing P0/P1 issues in the narrow modules listed above.
