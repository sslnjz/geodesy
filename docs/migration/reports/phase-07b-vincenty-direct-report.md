# Phase Report

## Phase Name

Phase 7B: Vincenty direct and convergence

## Objective

- Implement and harden the Vincenty direct solution for destination point and final bearing.
- Preserve public distance units in metres and angular inputs/outputs in degrees, with internal trigonometry in radians.
- Make direct zero-distance, non-finite inputs, height rejection, anti-meridian wrapping, and near-antipodal distance behaviour deterministic under test.
- Keep work scoped to the Phase 7 Vincenty files and migration bookkeeping.

## Modified Files

- `src/latlon/latlon_ellipsoidal_vincenty.h`
  - Retained the Phase 7 typed private result names used by the direct and inverse helpers.
  - No public Phase 7B API expansion was required.
- `src/latlon/latlon_ellipsoidal_vincenty.cpp`
  - Added named direct convergence tolerance and iteration limit.
  - Changed direct non-finite distance and bearing validation to `std::invalid_argument`.
  - Changed direct convergence failure to `std::runtime_error`.
  - Documented reduced-latitude and longitude wrapping boundaries in the direct calculation.
  - Explicitly wraps direct destination longitude before constructing the public destination point.
- `test/latlon_ellipsoidal_vincenty_unittest.cpp`
  - Added direct-solution tests for zero distance, non-finite distance/bearing rejection, height rejection, anti-meridian normalization, and near-antipodal distance.
- `docs/migration/phase_status.md`
  - Recorded Phase 7 completion and validation status.
- `docs/migration/reports/phase-07b-vincenty-direct-report.md`
  - Added this Phase 7B report.

## Added Or Modified Tests

- `latlon_ellipsoidal_vincenty_unittest.direct_returns_same_point_and_nan_bearing_for_zero_distance`
- `latlon_ellipsoidal_vincenty_unittest.direct_rejects_non_finite_distance_and_bearing`
- `latlon_ellipsoidal_vincenty_unittest.direct_rejects_starting_point_above_ellipsoid`
- `latlon_ellipsoidal_vincenty_unittest.direct_normalizes_destination_across_antimeridian`
- `latlon_ellipsoidal_vincenty_unittest.direct_handles_near_antipodal_distance`

## Commands Executed

- `git -c safe.directory=D:/workspace/github/geodesy status --short --branch`
- `cmake --build build --target latlon_ellipsoidal_vincenty_unittest --parallel`
- `cmake --build build --target latlon_ellipsoidal_vincenty_unittest --clean-first --parallel`
- `ctest --test-dir build -C Debug -R latlon_ellipsoidal_vincenty_unittest --output-on-failure`
- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`
- `cmake --build build --parallel`
- `ctest --test-dir build --output-on-failure`
- `ctest --test-dir build -C Debug --output-on-failure`
- `git -c safe.directory=D:/workspace/github/geodesy restore -- build`

## Test Results

- Initial focused build without elevated local SDK metadata access failed because MSBuild could not read `C:\Users\jbzs_\AppData\Local\Microsoft SDKs`.
- Initial focused Vincenty test after adding direct tests failed because non-finite direct inputs threw `std::runtime_error` instead of `std::invalid_argument`; the same run also showed stale inverse output from an older linked library.
- Clean focused build after allowing Windows SDK metadata access passed.
- Focused final check: `ctest --test-dir build -C Debug -R latlon_ellipsoidal_vincenty_unittest --output-on-failure` passed, 1/1 tests.
- Debug configure: sandboxed `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug` failed because the googletest subbuild could not read local Windows SDK metadata; rerun with local SDK metadata access passed.
- Debug build: sandboxed `cmake --build build --parallel` failed for the same Windows SDK metadata access reason; rerun with local SDK metadata access passed.
- Required raw CTest: `ctest --test-dir build --output-on-failure` failed because the Visual Studio multi-configuration build requires `-C Debug`; all 10 tests were reported `Not Run` with `Missing "-C <config>"`.
- Actual Debug CTest: `ctest --test-dir build -C Debug --output-on-failure` passed, 10/10 tests.

## Unfinished Items

- Raw CTest without `-C Debug` remains incompatible with the current Visual Studio multi-configuration build layout.
- Phase 8 UTM remains open.

## Next Phase Handoff

Execute Phase 8: UTM. Start from `docs/migration/phase_status.md`, this Phase 7B report, and the migration plan. Keep scope to `src/utm/utm.h`, `src/utm/utm.cpp`, `src/latlon/latlon_utm.h`, `src/latlon/latlon_utm.cpp`, and required UTM test/CMake enablement unless compilation requires a documented exception. Add or enable UTM tests first for ordinary conversions, Norway/Svalbard zone exceptions, invalid zones/easting/northing, formatting/parsing, and round-trip conversion. Then harden UTM forward/reverse projection with named constants and explicit standard exceptions, run focused UTM tests, Debug configure/build, raw CTest, and full Debug CTest with `-C Debug`.
