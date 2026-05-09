# Phase Report

## Phase Name

Phase 8B: LatLon to UTM forward and reverse projection

## Objective

- Harden latitude/longitude to UTM conversion in `LatLonUtm::toUtm()`.
- Harden UTM to latitude/longitude reverse projection in `Utm::toLatLon()`.
- Preserve JavaScript reference behaviour for ordinary WGS84 examples, Norway/Svalbard zone exceptions, and round trips.
- Keep public angular values in degrees, internal trigonometry in radians, and projected coordinates in metres.

## Modified Files

- `src/latlon/latlon_utm.cpp`
  - Added the required source header.
  - Corrected the longitude radians boundary by subtracting the central meridian in radians.
  - Corrected Kruger-series indexing to use one-based series terms over zero-based C++ coefficient arrays.
  - Added WGS84 datum fallback for default `LatLonUtm` points.
  - Added zone override validation and longitude 180 degree clamping to zone 60.
  - Applied Norway and Svalbard zone exceptions only for automatic zone selection.
  - Added named constants and comments for UTM false origins, central meridian scale, and Karney projection terms.
- `src/utm/utm.cpp`
  - Corrected the required source header.
  - Corrected inverse Kruger-series indexing and floating-point meridional-radius coefficients.
  - Added explicit inverse-projection iteration limit with `std::runtime_error` on non-convergence.
  - Normalized recovered longitude with `wrap180()`.
  - Kept Phase 8A parsing, formatting, and value-type validation behaviour.
- `test/utm_unittest.cpp`
  - Added focused forward projection, reverse projection, round-trip, Norway/Svalbard exception, and invalid-zone tests.

## Added Or Modified Tests

- `utm_unittest.latlon_to_utm_matches_reference_example`
- `utm_unittest.utm_to_latlon_matches_reference_example`
- `utm_unittest.round_trip_preserves_northern_and_southern_coordinates`
- `utm_unittest.applies_norway_and_svalbard_zone_exceptions`
- `utm_unittest.rejects_latitudes_outside_utm_limits_and_invalid_zone_override`

## Commands Executed

- `git -c safe.directory=D:/workspace/github/geodesy status --short`
- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`
- `cmake --build build --target utm_unittest --parallel`
- `ctest --test-dir build -C Debug -R utm_unittest --output-on-failure`
- `cmake --build build --target utm_unittest --config Debug --clean-first --parallel`
- `cmake --build build --target utm_unittest --config Debug --parallel`
- `cmake --build build --parallel`
- `ctest --test-dir build --output-on-failure`
- `ctest --test-dir build -C Debug --output-on-failure`
- `git -c safe.directory=D:/workspace/github/geodesy restore -- build`

## Test Results

- Initial `git status` without command-scoped `safe.directory` failed due to this checkout's dubious-ownership guard; command-scoped `safe.directory` status succeeded.
- Initial configure and build attempts inside the sandbox failed on local Windows SDK metadata access under `C:\Users\jbzs_\AppData\Local\Microsoft SDKs`; reruns with local SDK metadata access succeeded.
- Initial focused UTM test run failed after adding Phase 8B tests, exposing stale build outputs and then the projection defects fixed in this phase.
- Focused UTM clean build: `cmake --build build --target utm_unittest --config Debug --clean-first --parallel` passed.
- Focused UTM test: `ctest --test-dir build -C Debug -R utm_unittest --output-on-failure` passed, 1/1 tests.
- Debug configure: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug` passed.
- Debug build: `cmake --build build --parallel` passed.
- Required raw CTest: `ctest --test-dir build --output-on-failure` failed because the Visual Studio multi-configuration build requires `-C Debug`; all 11 tests were reported `Not Run`.
- Actual Debug CTest: `ctest --test-dir build -C Debug --output-on-failure` passed, 11/11 tests.

## Unfinished Items

- The required raw CTest command remains incompatible with the current Visual Studio multi-configuration build layout without `-C Debug`.
- MGRS conversion remains open for Phase 9.

## Next Phase Handoff

Continue with Phase 9: MGRS. Start from the migration plan, current phase status, and this Phase 8B report. Keep scope to `src/utm/mgrs.h`, `src/utm/mgrs.cpp`, `src/utm/utm_mgrs.h`, `src/utm/utm_mgrs.cpp`, `src/latlon/latlon_utm_mgrs.h`, `src/latlon/latlon_utm_mgrs.cpp`, and required MGRS tests unless compilation requires a documented exception. Add focused tests first for MGRS parsing, formatting precision, illegal grid letters, illegal bands, UTM bridge conversion, and round trips, then implement narrowly and run focused MGRS tests plus full Debug validation.
