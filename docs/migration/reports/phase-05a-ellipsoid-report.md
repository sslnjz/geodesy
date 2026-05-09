# Phase Report

## Phase Name

Phase 5A: Ellipsoid constants and LatLon -> Cartesian

## Objective

- Define stable C++17 ellipsoid/datum/reference-frame constant registries without heap-backed global references.
- Harden geodetic latitude/longitude/height to ECEF Cartesian conversion coverage for WGS84, height, equator, poles, and round-trip behaviour.
- Fix Cartesian-to-geodetic conversion at the poles, where longitude is undefined and Bowring's general formula divides by the distance from the minor axis.

## Modified Files

- `src/latlon/ellipsoids.h`
  - Replaced heap-backed global references with C++17 `inline const` registries for `g_ellipsoids`, `g_datums`, and `g_reference_frames`.
  - Preserved existing public registry names and WGS84 constants.
- `src/vector/cartesian.cpp`
  - Added invalid ellipsoid validation.
  - Added explicit geocentric-origin rejection with `std::domain_error`.
  - Added polar-axis handling for `Cartesian::toLatLon()`, returning +/-90 degrees latitude, 0 degrees longitude, and polar height `abs(z)-b`.
- `test/latlon_ellipsoidal_unittest.cpp`
  - Added tests for WGS84 constants, equator conversion, height conversion, Greenwich ECEF coordinates, Cartesian equator inverse conversion, pole conversion, and LatLon/Cartesian round-trip.
- `docs/migration/phase_status.md`
  - Recorded Phase 5A completion and validation status.

## Added Or Modified Tests

- `latlon_ellipsoidal_unittest.wgs84_ellipsoid_constants`
- `latlon_ellipsoidal_unittest.converts_equator_origin_to_cartesian`
- `latlon_ellipsoidal_unittest.converts_height_above_equator_to_cartesian`
- `latlon_ellipsoidal_unittest.converts_greenwich_reference_point_to_cartesian`
- `latlon_ellipsoidal_unittest.converts_cartesian_equator_to_latlon`
- `latlon_ellipsoidal_unittest.converts_poles_between_latlon_and_cartesian`
- `latlon_ellipsoidal_unittest.round_trips_latlon_height_through_cartesian`

## Commands Executed

- `cmake --build build --target latlon_ellipsoidal_unittest --parallel`
- `ctest --test-dir build -C Debug -R latlon_ellipsoidal_unittest --output-on-failure`
- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`
- `cmake --build build --parallel`
- `ctest --test-dir build --output-on-failure`
- `ctest --test-dir build -C Debug --output-on-failure`

## Test Results

- Focused red check: `ctest --test-dir build -C Debug -R latlon_ellipsoidal_unittest --output-on-failure` initially failed only on `converts_poles_between_latlon_and_cartesian`, confirming the missing polar inverse conversion branch.
- Focused final check: `ctest --test-dir build -C Debug -R latlon_ellipsoidal_unittest --output-on-failure` passed, 1/1 tests.
- Debug configure: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug` passed after allowing Windows SDK metadata access. The first sandboxed attempt failed because MSBuild could not read `C:\Users\jbzs_\AppData\Local\Microsoft SDKs`.
- Debug build: `cmake --build build --parallel` passed after allowing Windows SDK metadata access.
- Required raw CTest: `ctest --test-dir build --output-on-failure` failed because the Visual Studio multi-configuration build requires `-C Debug`; all tests were reported `Not Run` with `Missing "-C <config>"`.
- Actual Debug CTest: `ctest --test-dir build -C Debug --output-on-failure` ran 10 tests; 9 passed and `latlon_ellipsoidal_vincenty_unittest.example` failed on the existing degree-symbol encoding expectation mismatch recorded by earlier phase status.

## Unfinished Items

- None in Phase 5A scope.
- Existing baseline issue remains: Vincenty example string expectations use non-UTF-8 degree-symbol bytes and fail against current UTF-8 formatted output.
- Raw CTest without `-C Debug` remains incompatible with the current Visual Studio multi-configuration build layout.

## Next Phase Handoff

Suggested next phase start summary:

Execute Phase 5B: Cartesian -> LatLon edge hardening and ellipsoidal API cleanup. Start from `docs/migration/phase_status.md`, keep scope inside Phase 5 files unless the build requires otherwise, preserve public degrees/metres semantics, add tests before code for invalid Cartesian origin, invalid ellipsoid input, non-WGS84 ellipsoid inverse conversion, and any remaining round-trip boundary cases, then run focused ellipsoidal tests, Debug configure/build, and full Debug CTest.
