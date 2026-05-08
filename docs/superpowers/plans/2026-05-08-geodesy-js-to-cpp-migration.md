# Geodesy JS To C++ Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Migrate the Chris Veness geodesy JavaScript library concepts into the local C++17 codebase in controlled, testable phases until the full library surface is implemented and hardened.

**Architecture:** Build from foundational numeric utilities, DMS parsing, vector math, and base coordinate value types before migrating higher-risk geodesy algorithms. Preserve JavaScript behaviour through golden tests while using C++17 value types, explicit units, immutable constants, typed result structs, and predictable exceptions.

**Tech Stack:** C++17, CMake, GoogleTest, `std::vector`, `std::array`, `std::optional`, `std::future`, `std::async`, standard exceptions.

---

## Module Map

| JS module | C++ module | Priority | Migration target | Main risks |
|---|---|---|---|---|
| `dms.js` | `Dms`, `utils/algorithm`, `strutil` | Core | Angle parsing, formatting, wrapping, compass points | Unicode DMS input, signed zero, global separator |
| `vector3d.js` | `vector3d` | Core | 3D vector operations used by Cartesian and n-vector calculations | Zero vectors, cross-product direction, operator semantics |
| `latlon-spherical.js` | `LatLon`, `LatLonSpherical` | Core | Spherical distances, bearings, paths, rhumb lines, polygon operations | Poles, anti-meridian paths, coincident points |
| `latlon-ellipsoidal.js` | `LatLonEllipsoidal`, `Cartesian`, `Ellipsoid` | Core | Ellipsoidal geodetic coordinates and ECEF conversions | Degree/radian boundaries, height semantics, default WGS84 |
| `latlon-ellipsoidal-datum.js` | `LatLonEllipsoidalDatum`, `CartesianDatum`, `Datum`, `Transform` | Core | Historical datums and Helmert 7-parameter transforms | Datum constant precision, transform direction and units |
| `latlon-ellipsoidal-vincenty.js` | `LatLonEllipsoidalVincenty` | Core | Vincenty direct and inverse geodesic calculations | Near-antipodal non-convergence, iteration tolerances |
| `utm.js` | `LatLonUtm`, `Utm` | Core | Latitude/longitude and UTM conversions | Norway/Svalbard zones, Kruger series indexing |
| `mgrs.js` | `LatLonUtmMgrs`, `UtmMgrs`, `Mgrs` | Secondary | MGRS parsing, formatting, and UTM bridge | Precision truncation, invalid band and grid letters |
| `osgridref.js` | `LatLonOsGridRef`, `OsGridRef` | Secondary | British National Grid conversion and formatting | OSGB36/WGS84 datum conversion, meridional arc iteration |
| `latlon-ellipsoidal-referenceframe.js` | `LatLonEllipsoidalReferenceFrame`, `CartesianReferenceFrame`, `ReferenceFrame`, `HelmertTransforms` | Secondary, high risk | Dynamic reference frames and epoch-aware transforms | Decimal epoch handling, 14-parameter transforms, unsupported frame paths |
| `latlon-nvector-spherical.js` | `LatLonNvectorSpherical`, `NvectorSpherical` | Secondary | Spherical n-vector operations | Great-circle direction and segment edge cases |
| `latlon-nvector-ellipsoidal.js` | `LatLonNvectorEllipsoidal`, `NvectorEllipsoidal`, `NvectorCartesian`, `Ned` | Secondary | Ellipsoidal n-vector and NED operations | Datum and height inheritance semantics |

## Common Migration Steps

| Step | Applies to | Required work | Completion evidence |
|---|---|---|---|
| Interface refactor | Every migrated module | Preserve JavaScript-style method names such as `distanceTo`, `toLatLon`, `convertDatum`, and `destinationPoint`; map dynamic JavaScript inputs to explicit overloads or named factories. | Public headers compile standalone and tests call the intended API directly. |
| Data structure implementation | Every migrated module | Map `Object` to `struct` or `class`, `Array` to `std::vector` or `std::array`, JavaScript object registries to immutable constant tables or lookup maps. | No raw owning pointers or mutable global constants in newly migrated code. |
| Function implementation | Every migrated module | Port behaviour from the JavaScript reference examples first, then fill boundary cases. Keep public angles in degrees and internal trigonometry in radians. | Golden tests pass with meaningful tolerances. |
| Unit testing | Every migrated module | Add direct GoogleTest coverage for ordinary examples, invalid inputs, and relevant edge cases. | Focused module test passes and all available tests pass. |
| Validation | Every phase | Run CMake configure, build, focused tests, and full `ctest`. | Commands and outcomes are recorded in the phase summary. |

## Phase Plan

### Phase 0: Build And API Baseline

**Modules:** Build system, public headers, current tests.

**Priority:** Core.

**Files:**
- Inspect: `CMakeLists.txt`
- Inspect: `src/CMakeLists.txt`
- Inspect: `include/geodesy/*.h`
- Inspect: `test/CMakeLists.txt`
- Create or modify as needed: `test/js_reference_examples_*_unittest.cpp`

- [ ] Confirm the current Debug build command:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
```

- [ ] Confirm the current build command:

```bash
cmake --build build --parallel
```

- [ ] Confirm the current test command:

```bash
ctest --test-dir build --output-on-failure
```

- [ ] Add a test organization for JavaScript reference examples before broad algorithm changes.

**Completion standard:** Current build and tests are reproducible, and every later phase has a clear place for JavaScript golden examples.

**Pitfalls:** Do not begin by reshaping the whole source tree. The first objective is stable verification, not cosmetic layout.

### Phase 1: DMS And Angle Utilities

**Modules:** `dms.js`, shared angle helpers.

**Priority:** Core.

**Files:**
- Modify: `src/dms/dms.h`
- Modify: `src/dms/dms.cpp`
- Modify: `src/utils/algorithm.h`
- Test: `test/dms_unittest.cpp`

- [ ] Refactor the `Dms` public interface around explicit numeric and string inputs.
- [ ] Implement or verify `parse`, `toDms`, `toLat`, `toLon`, `toBearing`, `compassPoint`, `wrap90`, `wrap180`, and `wrap360`.
- [ ] Replace new mutable constants with `constexpr` or controlled state where practical.
- [ ] Add tests for decimal degrees, DMS with hemisphere suffixes, negative values, wrapping, and compass precision.
- [ ] Run the focused DMS test target.

**Completion standard:** JavaScript DMS examples and C++ boundary tests pass.

**Pitfalls:** Unicode degree, prime, and double-prime characters must be tested. Signed zero and separator behaviour can affect exact output strings.

### Phase 2: Vector3d

**Modules:** `vector3d.js`.

**Priority:** Core.

**Files:**
- Modify: `src/vector/vector3d.h`
- Modify: `src/vector/vector3d.cpp`
- Test: `test/vector3d_unittest.cpp`

- [ ] Correct public getters and operator declarations so the header is portable C++17.
- [ ] Implement or verify plus, minus, scalar multiplication, scalar division, unary minus, dot product, cross product, unit vector, angle, and rotation.
- [ ] Add tests for unary minus, cross-product sign, zero-length normalization, angle sign, and rotation.
- [ ] Run the focused vector test target.

**Completion standard:** Vector operations are numerically safe enough for Cartesian and n-vector modules.

**Pitfalls:** Cross-product operand order changes signs. Zero-vector handling must be explicit rather than accidental.

### Phase 3: Base LatLon

**Modules:** Shared `LatLon` base used by spherical and ellipsoidal points.

**Priority:** Core.

**Files:**
- Modify: `src/latlon/latlon.h`
- Modify: `src/latlon/latlon.cpp`
- Test: `test/latlon_unittest.cpp`

- [ ] Define the base coordinate as a value type with latitude and longitude stored in degrees.
- [ ] Preserve JavaScript-style constructor or parse behaviour through explicit overloads.
- [ ] Implement or verify `toString`, `toGeoJSON`, equality, and longitude normalization.
- [ ] Add tests for numeric input, string input, invalid coordinates, anti-meridian longitudes, and formatting.

**Completion standard:** Base coordinate parsing and formatting are stable before algorithm modules depend on it.

**Pitfalls:** JavaScript accepts dynamic object shapes. In C++17, prefer overloads and adapters rather than a large implicit variant-based constructor.

### Phase 4: Spherical LatLon

**Modules:** `latlon-spherical.js`.

**Priority:** Core.

**Files:**
- Modify: `src/latlon/latlon_spherical.h`
- Modify: `src/latlon/latlon_spherical.cpp`
- Test: `test/latlon_spherical_unittest.cpp`

- [ ] Implement or verify distance, initial bearing, final bearing, midpoint, intermediate point, destination point, and path intersection.
- [ ] Implement or verify rhumb-line distance, bearing, destination, midpoint, and cross/along-track operations.
- [ ] Implement or verify polygon area with pole and anti-meridian handling.
- [ ] Add JavaScript reference examples plus tests for poles, anti-meridian crossing, coincident points, and zero distance.

**Completion standard:** Spherical module examples and boundary cases pass against meaningful distance and angular tolerances.

**Pitfalls:** Area and intersection edge cases can appear correct for ordinary examples while failing near poles or the anti-meridian.

### Phase 5: Ellipsoidal Coordinates And Cartesian

**Modules:** `latlon-ellipsoidal.js`.

**Priority:** Core.

**Files:**
- Modify: `src/latlon/ellipsoids.h`
- Modify: `src/latlon/latlon_ellipsoidal.h`
- Modify: `src/latlon/latlon_ellipsoidal.cpp`
- Modify: `src/vector/cartesian.h`
- Modify: `src/vector/cartesian.cpp`
- Test: `test/latlon_ellipsoidal_unittest.cpp`

- [ ] Define `Ellipsoid` constants with sufficient reference precision.
- [ ] Implement or verify geodetic latitude/longitude/height to ECEF Cartesian conversion.
- [ ] Implement or verify Cartesian to geodetic conversion.
- [ ] Add tests for WGS84, non-zero height, poles, equator, and round-trip conversion.

**Completion standard:** Ellipsoidal and Cartesian conversion is reliable before datum and projection modules build on it.

**Pitfalls:** Public degrees and internal radians must not leak across the boundary. Height is metres above ellipsoid, not orthometric height.

### Phase 6: Datum And Helmert 7-Parameter Transforms

**Modules:** `latlon-ellipsoidal-datum.js`.

**Priority:** Core.

**Files:**
- Modify: `src/latlon/ellipsoids.h`
- Modify: `src/latlon/latlon_ellipsoidal_datum.h`
- Modify: `src/latlon/latlon_ellipsoidal_datum.cpp`
- Modify: `src/vector/cartesian_datum.h`
- Modify: `src/vector/cartesian_datum.cpp`
- Test: `test/latlon_ellipsoidal_datum_unittest.cpp`

- [ ] Define `Datum` and `Transform` as immutable value-oriented data structures.
- [ ] Implement or verify WGS84, OSGB36, NAD83, and other JavaScript reference datums.
- [ ] Implement or verify `convertDatum` using Cartesian conversion and Helmert 7-parameter transform.
- [ ] Add tests for WGS84 to OSGB36 and back, datum identity, and invalid datum paths.

**Completion standard:** Datum conversion examples match the JavaScript reference within geodesy-appropriate tolerances.

**Pitfalls:** Helmert rotations are in arc-seconds and scale is in parts per million. These conversions should happen exactly once at the calculation boundary.

### Phase 7: Vincenty Geodesics

**Modules:** `latlon-ellipsoidal-vincenty.js`.

**Priority:** Core.

**Files:**
- Modify: `src/latlon/latlon_ellipsoidal_vincenty.h`
- Modify: `src/latlon/latlon_ellipsoidal_vincenty.cpp`
- Test: `test/latlon_ellipsoidal_vincenty_unittest.cpp`

- [ ] Define typed inverse and direct result structs for distance, bearings, and destination data.
- [ ] Implement or verify Vincenty inverse solution.
- [ ] Implement or verify Vincenty direct solution.
- [ ] Add tests for distance, initial bearing, final bearing, destination, coincident points, and near-antipodal non-convergence.

**Completion standard:** Vincenty examples and failure cases are deterministic.

**Pitfalls:** Do not return plausible-looking values on non-convergence. Use explicit iteration limits and named tolerances.

### Phase 8: UTM

**Modules:** `utm.js`.

**Priority:** Core.

**Files:**
- Modify: `src/utm/utm.h`
- Modify: `src/utm/utm.cpp`
- Modify: `src/latlon/latlon_utm.h`
- Modify: `src/latlon/latlon_utm.cpp`
- Test: add or enable UTM tests in `test/CMakeLists.txt`

- [ ] Define `Utm` as a value type with `zone`, `Hemisphere`, `easting`, `northing`, `datum`, `convergence`, and `scale` where required.
- [ ] Implement or verify UTM parsing and formatting.
- [ ] Implement or verify latitude/longitude to UTM conversion.
- [ ] Implement or verify UTM to latitude/longitude conversion.
- [ ] Add tests for ordinary zones, Norway and Svalbard exceptions, invalid zones, invalid easting, invalid northing, and round-trip conversion.

**Completion standard:** UTM forward and reverse projection matches JavaScript reference examples and boundary expectations.

**Pitfalls:** The Kruger series is index-sensitive. Zone exceptions and northing range checks are common sources of hard failures.

### Phase 9: MGRS

**Modules:** `mgrs.js`.

**Priority:** Secondary.

**Files:**
- Modify: `src/utm/mgrs.h`
- Modify: `src/utm/mgrs.cpp`
- Modify: `src/utm/utm_mgrs.h`
- Modify: `src/utm/utm_mgrs.cpp`
- Modify: `src/latlon/latlon_utm_mgrs.h`
- Modify: `src/latlon/latlon_utm_mgrs.cpp`
- Test: add or enable MGRS tests in `test/CMakeLists.txt`

- [ ] Define MGRS grid zone, band, 100 km square, easting, northing, and precision fields.
- [ ] Implement or verify MGRS parsing and formatting.
- [ ] Implement or verify MGRS to UTM and UTM to MGRS conversion.
- [ ] Add tests for precision levels, illegal grid letters, illegal bands, and round-trip conversion.

**Completion standard:** MGRS examples and UTM round-trips pass.

**Pitfalls:** Precision-dependent digits are truncation-sensitive. Invalid letters must be rejected rather than normalized.

### Phase 10: OS Grid Reference

**Modules:** `osgridref.js`.

**Priority:** Secondary.

**Files:**
- Modify: `src/utm/osgridref.h`
- Modify: `src/utm/osgridref.cpp`
- Modify: `src/latlon/latlon_osgridref.h`
- Modify: `src/latlon/latlon_osgridref.cpp`
- Test: add or enable OS Grid tests in `test/CMakeLists.txt`

- [ ] Define national grid constants for Airy 1830, true origin, false origin, and scale.
- [ ] Implement or verify OS grid reference parsing.
- [ ] Implement or verify OS grid reference formatting.
- [ ] Implement or verify WGS84/OSGB36 conversion bridge through datum conversion.
- [ ] Add tests for multiple formatting precisions and malformed grid references.

**Completion standard:** OS Grid examples pass and datum bridge behaviour is explicit.

**Pitfalls:** Meridional arc iteration must preserve precision. OSGB36 and WGS84 conversion cannot be treated as a simple spherical transform.

### Phase 11: Reference Frames And Epoch Transforms

**Modules:** `latlon-ellipsoidal-referenceframe.js`.

**Priority:** Secondary, high risk.

**Files:**
- Modify: `src/latlon/latlon_ellipsoidal_referenceframe_txparams.h`
- Modify: `src/latlon/latlon_ellipsoidal_referenceframe.h`
- Modify: `src/latlon/latlon_ellipsoidal_referenceframe.cpp`
- Modify: `src/vector/cartesian_referenceFrame.h`
- Modify: `src/vector/cartesian_referenceFrame.cpp`
- Test: `test/latlon_ellipsoidal_referenceframe_unittest.cpp`

- [ ] Replace string-only epoch arithmetic with decimal-year semantics.
- [ ] Define 14-parameter Helmert transform data with named fields or `std::array<double, 7>` pairs.
- [ ] Implement or verify reference-frame conversion chains.
- [ ] Add tests for direct transform, chained transform, epoch rate adjustment, unsupported transforms, and metadata preservation.

**Completion standard:** Reference-frame conversion preserves frame and epoch metadata and fails predictably for unsupported paths.

**Pitfalls:** Existing optional misuse and pointer-throw behaviour must be fixed under test. Decimal epochs such as `2010.0` must not be truncated.

### Phase 12: Spherical N-Vector

**Modules:** `latlon-nvector-spherical.js`.

**Priority:** Secondary.

**Files:**
- Modify: `src/vector/nvector_spherical.h`
- Modify: `src/vector/nvector_spherical.cpp`
- Modify: `src/latlon/latlon_nvector_spherical.h`
- Modify: `src/latlon/latlon_nvector_spherical.cpp`
- Test: `test/latlon_nvector_spherical_unittest.cpp`

- [ ] Implement or verify latitude/longitude to n-vector conversion.
- [ ] Implement or verify n-vector to latitude/longitude conversion.
- [ ] Implement or verify great-circle, path intersection, nearest point, triangulation, trilateration, polygon enclosure, area, and mean point operations.
- [ ] Add tests for poles, anti-meridian paths, segment endpoints, and sign conventions.

**Completion standard:** Spherical n-vector behaviour matches JavaScript examples and does not regress spherical trig tests.

**Pitfalls:** Great-circle direction and cross-product sign are easy to invert. Segment boundary logic must be tested separately from infinite great-circle logic.

### Phase 13: Ellipsoidal N-Vector And NED

**Modules:** `latlon-nvector-ellipsoidal.js`.

**Priority:** Secondary.

**Files:**
- Modify: `src/vector/nvector_ellipsoidal.h`
- Modify: `src/vector/nvector_ellipsoidal.cpp`
- Modify: `src/vector/nvector_cartesian.h`
- Modify: `src/vector/nvector_cartesian.cpp`
- Modify: `src/vector/ned.h`
- Modify: `src/vector/ned.cpp`
- Modify: `src/latlon/latlon_nvector_ellipsoidal.h`
- Modify: `src/latlon/latlon_nvector_ellipsoidal.cpp`
- Test: `test/latlon_nvector_ellipsoidal_unittest.cpp`

- [ ] Implement or verify ellipsoidal latitude/longitude to n-vector conversion.
- [ ] Implement or verify n-vector Cartesian conversion.
- [ ] Implement or verify NED delta and destination from delta.
- [ ] Add tests for height, datum, round-trip conversion, and local NED movements.

**Completion standard:** Ellipsoidal n-vector and NED operations are compatible with ellipsoidal Cartesian semantics.

**Pitfalls:** NED is local north/east/down in metres. It must not be confused with ECEF x/y/z metres.

### Phase 14: API Compatibility Layer

**Modules:** Public API surface across all modules.

**Priority:** Core.

**Files:**
- Modify: `include/geodesy/*.h`
- Modify: corresponding `src/**/*.h`
- Test: public header compile tests or install-tree tests

- [ ] Preserve JavaScript-style method names where they are already part of the C++ public API.
- [ ] Add explicit named factories for parse and conversion entrypoints where constructors become ambiguous.
- [ ] Normalize exception taxonomy to `std::invalid_argument`, `std::domain_error`, and `std::runtime_error`.
- [ ] Ensure every public header compiles independently.

**Completion standard:** Public C++17 API is stable, documented through names and tests, and close to JavaScript concepts without dynamic runtime tricks.

**Pitfalls:** Avoid a broad `std::variant` constructor that tries to emulate every JavaScript input shape. Use focused overloads and adapters.

### Phase 15: Full Library Convergence

**Modules:** All modules.

**Priority:** Core.

**Files:**
- Modify only files required by failing tests or documented compatibility gaps.

- [ ] Run Debug configure, build, and all tests.
- [ ] Run Release configure, build, and all tests.
- [ ] Compare module outputs against JavaScript reference examples.
- [ ] Remove or quarantine fragile global mutable state in touched modules.
- [ ] Record remaining limitations in project documentation.

**Completion standard:** Full C++17 library functionality is implemented, tested, and aligned with the JavaScript reference behaviour as closely as practical.

**Pitfalls:** Do not mix new features into convergence work. This phase is for behaviour alignment, stability, and documented limitations.

## Async Strategy

| JavaScript concept | C++17 design | Type | Description | Notes |
|---|---|---|---|---|
| Synchronous geodesy methods | Synchronous C++ member functions | Core library | All parse, convert, and format operations return values synchronously. | Keep the core library deterministic and thread-agnostic. |
| Promise-like batch flow | `std::future<std::vector<T>>` | Application adapter | Batch point conversion or parsing can run outside the core library. | Do not change core type semantics for async convenience. |
| Callback completion | `std::function<void(Result)>` | UI or batch layer | Optional completion notification for application code. | In C++17, represent result boundaries explicitly because `std::expected` is unavailable. |
| Cancellation and progress | Token plus progress callback | Optional adapter | Useful for large batch conversion or import workflows. | Keep this outside geodesy value types. |

## Validation Commands

Use these commands as the default baseline unless a phase identifies a more specific target:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

For release checks:

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure
```

## Notes

- Core library APIs should remain synchronous and value-oriented.
- Public angular inputs and outputs use degrees unless explicitly named otherwise.
- Internal trigonometric calculations use radians.
- Distances and heights use metres unless a specific API documents another unit.
- Every behavioural change must be covered by tests before being considered complete.
