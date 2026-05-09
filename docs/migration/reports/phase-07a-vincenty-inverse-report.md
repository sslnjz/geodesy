# Phase Report

## Phase Name

Phase 7A: Vincenty inverse

## Objective

- Implement and harden the Vincenty inverse solution for distance, initial bearing, and final bearing.
- Preserve public angular inputs and outputs in degrees, with internal trigonometry in radians.
- Make coincident points, exact antipodal points, near-antipodal non-convergence, and iteration convergence deterministic under test.
- Keep Phase 7B direct-solution work out of scope except for type-name consistency required by the shared private result declarations.

## Modified Files

- `src/latlon/latlon_ellipsoidal_vincenty.h`
  - Renamed private forward-declared result structs from short names to `VincentyInverseResult` and `VincentyDirectResult`.
  - Corrected Vincenty documentation typo in the module comment.
- `src/latlon/latlon_ellipsoidal_vincenty.cpp`
  - Added named inverse convergence constants.
  - Preserved the JavaScript reference inverse iteration structure while using `std::runtime_error` for inverse non-convergence.
  - Made `distanceTo()`, `initialBearingTo()`, and `finalBearingTo()` return `NaN` for inverse convergence failure.
  - Used the reference `1e-24` `sinSqSigma` threshold for coincident and exactly antipodal handling.
  - Restored the required file header formatting and removed the old BOM from the file.
- `test/latlon_ellipsoidal_vincenty_unittest.cpp`
  - Added the required file header.
  - Corrected UTF-8 degree-symbol expectations.
  - Added inverse-focused tests for coincident points, exact antipodal meridional output, and near-antipodal non-convergence.
- `docs/migration/phase_status.md`
  - Recorded Phase 7A completion and validation status.

## Added Or Modified Tests

- `latlon_ellipsoidal_vincenty_unittest.example`
  - Corrected expected degree-symbol strings to UTF-8.
- `latlon_ellipsoidal_vincenty_unittest.inverse_returns_nan_bearings_for_coincident_points`
- `latlon_ellipsoidal_vincenty_unittest.inverse_handles_exact_antipodal_meridional_path`
- `latlon_ellipsoidal_vincenty_unittest.inverse_reports_nan_for_near_antipodal_non_convergence`

## Commands Executed

- `git -c safe.directory=D:/workspace/github/geodesy status --short --branch`
- `cmake --build build --target latlon_ellipsoidal_vincenty_unittest --parallel`
- `ctest --test-dir build -C Debug -R latlon_ellipsoidal_vincenty_unittest --output-on-failure`
- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`
- `cmake --build build --parallel`
- `ctest --test-dir build --output-on-failure`
- `ctest --test-dir build -C Debug --output-on-failure`

## Test Results

- Initial focused Vincenty test after adding the near-antipodal regression failed because `initialBearingTo()` and `finalBearingTo()` allowed the inverse non-convergence exception to escape.
- Focused final check: `ctest --test-dir build -C Debug -R latlon_ellipsoidal_vincenty_unittest --output-on-failure` passed, 1/1 tests.
- Debug configure: sandboxed `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug` failed because MSBuild could not read `C:\Users\jbzs_\AppData\Local\Microsoft SDKs`; rerun with local SDK metadata access passed.
- Debug build: `cmake --build build --parallel` passed after allowing Windows SDK metadata access.
- Required raw CTest: `ctest --test-dir build --output-on-failure` failed because the Visual Studio multi-configuration build requires `-C Debug`; all 10 tests were reported `Not Run` with `Missing "-C <config>"`.
- Actual Debug CTest: `ctest --test-dir build -C Debug --output-on-failure` passed, 10/10 tests.

## Unfinished Items

- Phase 7B remains open for Vincenty direct-solution hardening, including destination and final-bearing-on-path convergence checks.
- Raw CTest without `-C Debug` remains incompatible with the current Visual Studio multi-configuration build layout.

## Next Phase Handoff

Suggested next phase start summary:

Execute Phase 7B: Vincenty direct solution and convergence. Start from `docs/migration/phase_status.md` and `docs/migration/reports/phase-07a-vincenty-inverse-report.md`, keep scope inside `src/latlon/latlon_ellipsoidal_vincenty.h`, `src/latlon/latlon_ellipsoidal_vincenty.cpp`, and `test/latlon_ellipsoidal_vincenty_unittest.cpp` unless compilation requires a documented exception. Add or repair tests first for `destinationPoint()`, `finalBearingOn()`, zero distance, invalid distance/bearing, height rejection, longitude normalization, and direct-iteration convergence failure if a reachable case exists. Then harden `direct()` with named tolerances and standard exception semantics, run focused Vincenty tests, Debug configure/build, raw CTest, and full Debug CTest with `-C Debug`.
