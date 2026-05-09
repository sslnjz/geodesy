# Phase Report

## Phase Name

Phase 6: Datum And Helmert 7-Parameter Transforms

## Objective

- Harden historical datum conversion against the JavaScript reference flow.
- Preserve public latitude/longitude in degrees and height/ECEF coordinates in metres.
- Apply Helmert translation, scale, and rotation parameters exactly once at the Cartesian transform boundary.
- Lock WGS84 identity, WGS84<->OSGB36 conversion, chained datum conversion, and invalid datum paths with focused tests.

## Modified Files

- `src/latlon/latlon_ellipsoidal_datum.h`
  - Made `datum()` const and value-returning.
  - Moved inline datum/ellipsoid getters into the implementation file.
- `src/latlon/latlon_ellipsoidal_datum.cpp`
  - Added explicit datum ellipsoid validation.
  - Preserved WGS84 as the default datum when no datum is supplied.
  - Validated source and target datums before geodetic-to-Cartesian conversion.
- `src/vector/cartesian_datum.h`
  - Made `datum()` const and made `convertDatum()` const.
  - Marked conversion return values `[[nodiscard]]`.
- `src/vector/cartesian_datum.cpp`
  - Added explicit datum ellipsoid validation.
  - Removed library `stdout` output from deprecated `toLatLon()` datum parameter handling.
  - Implemented WGS84 identity as a no-op Cartesian conversion.
  - Replaced mutating transform inversion with a local inverse transform.
  - Kept non-WGS84 datum pairs chained through WGS84.
- `test/latlon_ellipsoidal_datum_unittest.cpp`
  - Added Phase 6 tests for WGS84 identity, WGS84<->OSGB36 numerical conversion, Cartesian identity, and invalid datum paths.
- `docs/migration/phase_status.md`
  - Recorded Phase 6 completion and validation status.

## Added Or Modified Tests

- `latlon_ellipsoidal_datum_unittest.convert_datum_preserves_identity_for_wgs84`
- `latlon_ellipsoidal_datum_unittest.converts_wgs84_to_osgb36_and_back_numerically`
- `latlon_ellipsoidal_datum_unittest.cartesian_identity_datum_does_not_apply_false_transform`
- `latlon_ellipsoidal_datum_unittest.rejects_unrecognised_datum_paths`

## Commands Executed

- `git -c safe.directory=D:/workspace/github/geodesy status --short --branch`
- `cmake --build build --target latlon_ellipsoidal_datum_unittest --parallel`
- `ctest --test-dir build -C Debug -R latlon_ellipsoidal_datum_unittest --output-on-failure`
- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`
- `cmake --build build --parallel`
- `ctest --test-dir build --output-on-failure`
- `ctest --test-dir build -C Debug --output-on-failure`

## Test Results

- Initial focused build after adding tests failed before implementation because `LatLonEllipsoidalDatum::datum()` and `CartesianDatum::convertDatum()` were not const-friendly and because WGS84 identity still depended on datum/transform truthiness.
- Focused final check: `ctest --test-dir build -C Debug -R latlon_ellipsoidal_datum_unittest --output-on-failure` passed, 1/1 tests.
- Debug configure: sandboxed `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug` failed because MSBuild could not read `C:\Users\jbzs_\AppData\Local\Microsoft SDKs`; rerun with local SDK metadata access passed.
- Debug build: `cmake --build build --parallel` passed after allowing Windows SDK metadata access. Linker warnings for missing GoogleTest PDB files remain non-fatal.
- Required raw CTest: `ctest --test-dir build --output-on-failure` failed because the Visual Studio multi-configuration build requires `-C Debug`; all 10 tests were reported `Not Run` with `Missing "-C <config>"`.
- Actual Debug CTest: `ctest --test-dir build -C Debug --output-on-failure` ran 10 tests; 9 passed and `latlon_ellipsoidal_vincenty_unittest.example` failed on the existing degree-symbol encoding expectation mismatch recorded by earlier phase status.

## Unfinished Items

- None in Phase 6 scope.
- Existing baseline issue remains: Vincenty example string expectations use replacement-character degree-symbol bytes and fail against current UTF-8 formatted output.
- Raw CTest without `-C Debug` remains incompatible with the current Visual Studio multi-configuration build layout.

## Next Phase Handoff

Suggested next phase start summary:

Execute Phase 7: Vincenty Geodesics. Start from `docs/migration/phase_status.md` and `docs/migration/reports/phase-06-datum-report.md`, keep scope inside `src/latlon/latlon_ellipsoidal_vincenty.h`, `src/latlon/latlon_ellipsoidal_vincenty.cpp`, and `test/latlon_ellipsoidal_vincenty_unittest.cpp` unless compilation requires a documented exception. Add or repair tests first for inverse distance, initial/final bearings, direct destination, coincident points, near-antipodal non-convergence, and the existing degree-symbol encoding mismatch, then harden Vincenty direct/inverse iteration and run focused Vincenty tests, Debug configure/build, raw CTest, and full Debug CTest.
