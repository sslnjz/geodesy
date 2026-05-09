# Phase Report

## Phase Name

Phase 4A: Spherical distance, bearing, midpoint, destination, and intersection

## Objective

- Harden the great-circle subset of `LatLonSpherical`.
- Preserve Chris Veness spherical JavaScript behaviour for distance, initial bearing, final bearing, midpoint, intermediate point, destination point, and path intersection.
- Add focused boundary coverage for zero distance, coincident points, anti-meridian crossings, poles, invalid radius input, and polar intersections.

## Modified Files

- `src/latlon/latlon_spherical.h`
  - Added direct standard-library includes used by the public declarations.
- `src/latlon/latlon_spherical.cpp`
  - Added a local positive finite radius validator for great-circle distance and destination calculations.
  - Centralized haversine central-angle calculation and clamped the intermediate value before square-root evaluation.
  - Reused the central-angle helper in `distanceTo`, `intermediatePointTo`, and `intersection`.
  - Replaced an exact zero comparison in intersection ambiguity handling with an explicit angular tolerance.
  - Replaced a stale final-bearing comment with a concise description of the reverse-bearing calculation.
- `test/latlon_spherical_unittest.cpp`
  - Added great-circle tests for zero distance, coincident points, anti-meridian crossing, and pole handling.
  - Added invalid-radius coverage for `distanceTo` and `destinationPoint`.
  - Replaced a polar intersection string check with numeric latitude and finite-longitude assertions, because longitude is undefined at the exact pole.
- `docs/migration/phase_status.md`
  - Recorded Phase 4A completion and left remaining Phase 4 slices open.
- `docs/migration/reports/phase-04a-spherical-dist-report.md`
  - Added this report.

## Added Or Modified Tests

- Added `great_circle_zero_distance_and_coincident_points`.
- Added `great_circle_crosses_anti_meridian`.
- Added `great_circle_handles_poles`.
- Extended `dist_brng_dest_fails` for invalid great-circle radius input.
- Adjusted the polar intersection check to assert meaningful numerical behaviour instead of an arbitrary pole longitude string.

## Commands Executed

- `git -c safe.directory=D:/workspace/github/geodesy status --short --branch`
- `git -c safe.directory=D:/workspace/github/geodesy switch -c <phase-4a-branch>`
- `ctest --test-dir build -C Debug -R latlon_spherical_unittest --output-on-failure`
- `cmake --build build --parallel --target latlon_spherical_unittest`
- `ctest --test-dir build -C Debug -R latlon_spherical_unittest --output-on-failure`
- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`
- `cmake --build build --parallel`
- `ctest --test-dir build --output-on-failure`
- `ctest --test-dir build -C Debug --output-on-failure`
- `cmake --build build --parallel --target latlon_spherical_unittest --clean-first`
- `ctest --test-dir build -C Debug -R latlon_spherical_unittest --output-on-failure`

## Test Results

- Initial focused spherical CTest failed on the existing polar intersection string expectation with a trailing apostrophe.
- Focused spherical build target passed after allowing the Windows toolchain to access local SDK metadata.
- Focused spherical CTest with `-C Debug` passed, 1/1.
- Final clean rebuild of the focused spherical target passed, followed by focused spherical CTest passing, 1/1.
- Required Debug configure command passed after allowing the Windows toolchain to access local SDK metadata. CMake still reports existing deprecation and developer warnings outside Phase 4A scope.
- Required Debug build command passed.
- Required raw CTest command `ctest --test-dir build --output-on-failure` failed under the Visual Studio multi-config generator because no configuration was selected.
- Full Debug CTest with `-C Debug` ran 10 tests, 9 passed and 1 failed:
  - `latlon_ellipsoidal_vincenty_unittest.example` still fails due to existing non-UTF-8 degree-byte expectations. This is outside Phase 4A scope and was already recorded in prior phase reports.

## Unfinished Items

- Remaining Phase 4 work is not completed in this slice: rhumb-line operations, cross/along-track review beyond existing coverage, and polygon area hardening.
- Full raw CTest still needs `-C Debug` with the current Visual Studio multi-config build.
- Existing Vincenty encoding-related test failure remains outside Phase 4A scope.

## Next Phase Handoff

- Continue with the next Phase 4 slice rather than Phase 5.
- Recommended next slice: Phase 4B for spherical rhumb-line and track operations.
- Keep edits limited to `src/latlon/latlon_spherical.h`, `src/latlon/latlon_spherical.cpp`, and `test/latlon_spherical_unittest.cpp` unless compilation requires a documented exception.
- Preserve Phase 4A behaviour: finite positive radius validation, clamped haversine central angle, NaN bearings for coincident points, normalized anti-meridian output, and pole assertions that do not depend on arbitrary longitude.
