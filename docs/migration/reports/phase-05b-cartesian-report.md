# Phase Report

## Phase Name

Phase 5B: Cartesian -> LatLon and round-trip

## Objective

- Harden ECEF Cartesian to geodetic latitude/longitude/height conversion boundaries.
- Preserve public degrees/metres semantics and internal radians/metres calculations.
- Lock reverse conversion behaviour with tests for invalid inputs, polar-axis handling, non-WGS84 ellipsoids, and round-trip conversion.

## Modified Files

- `src/vector/cartesian.cpp`
  - Added explicit rejection of non-finite ECEF x/y/z coordinates.
  - Tightened ellipsoid validation before Bowring reverse conversion.
  - Preserved geocentric-origin rejection and polar-axis handling.
- `src/latlon/latlon_ellipsoidal.cpp`
  - Added explicit non-finite height rejection before forward ECEF conversion.
  - Added ellipsoid validation before using semi-major axis and flattening.
- `test/latlon_ellipsoidal_unittest.cpp`
  - Added Phase 5B tests for non-finite Cartesian coordinates, geocentric origin, invalid ellipsoid input, non-finite height, and Airy1830/OSGB36 round-trip conversion.
  - Preserved Phase 5A WGS84, height, pole, and round-trip tests.
- `docs/migration/phase_status.md`
  - Recorded Phase 5B completion and validation status.

## Added Or Modified Tests

- `latlon_ellipsoidal_unittest.rejects_non_finite_cartesian_coordinates`
- `latlon_ellipsoidal_unittest.rejects_invalid_cartesian_origin_and_ellipsoid`
- `latlon_ellipsoidal_unittest.rejects_non_finite_geodetic_height`
- `latlon_ellipsoidal_unittest.round_trips_non_wgs84_ellipsoid_through_cartesian`

## Commands Executed

- `git -c safe.directory=D:/workspace/github/geodesy status --short`
- `cmake --build build --target latlon_ellipsoidal_unittest --parallel`
- `ctest --test-dir build -C Debug -R latlon_ellipsoidal_unittest --output-on-failure`
- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`
- `cmake --build build --parallel`
- `ctest --test-dir build --output-on-failure`
- `ctest --test-dir build -C Debug --output-on-failure`

## Test Results

- Initial focused check: `ctest --test-dir build -C Debug -R latlon_ellipsoidal_unittest --output-on-failure` failed before the Phase 5B fix. It exposed missing polar-axis handling in the rebuilt target and missing rejection for geocentric origin / invalid ellipsoid input.
- Focused final check: `ctest --test-dir build -C Debug -R latlon_ellipsoidal_unittest --output-on-failure` passed, 1/1 tests.
- Debug configure: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug` passed after allowing Windows SDK metadata access. The sandboxed attempt failed because MSBuild could not read `C:\Users\jbzs_\AppData\Local\Microsoft SDKs`.
- Debug build: `cmake --build build --parallel` passed after allowing Windows SDK metadata access. A transient focused-link retry failed on `latlon_ellipsoidal_unittest.pdb`; deleting the generated stale PDB and rebuilding resolved it.
- Required raw CTest: `ctest --test-dir build --output-on-failure` failed because the Visual Studio multi-configuration build requires `-C Debug`; all 10 tests were reported `Not Run` with `Missing "-C <config>"`.
- Actual Debug CTest: `ctest --test-dir build -C Debug --output-on-failure` ran 10 tests; 9 passed and `latlon_ellipsoidal_vincenty_unittest.example` failed on the existing degree-symbol encoding expectation mismatch recorded by earlier phase status.

## Unfinished Items

- None in Phase 5B scope.
- Existing baseline issue remains: Vincenty example string expectations use non-UTF-8 degree-symbol bytes and fail against current UTF-8 formatted output.
- Raw CTest without `-C Debug` remains incompatible with the current Visual Studio multi-configuration build layout.

## Next Phase Handoff

Suggested next phase start summary:

Execute Phase 6: Datum And Helmert 7-Parameter Transforms. Start from `docs/migration/phase_status.md` and `docs/migration/reports/phase-05b-cartesian-report.md`, keep scope inside `src/latlon/ellipsoids.h`, `src/latlon/latlon_ellipsoidal_datum.h`, `src/latlon/latlon_ellipsoidal_datum.cpp`, `src/vector/cartesian_datum.h`, `src/vector/cartesian_datum.cpp`, and `test/latlon_ellipsoidal_datum_unittest.cpp` unless compilation requires a documented exception. Add tests first for WGS84-to-OSGB36 and reverse datum conversion, datum identity, invalid datum paths, transform direction/units, and preservation of height/metres semantics, then implement or harden Helmert 7-parameter conversion and run focused datum tests, Debug configure/build, raw CTest, and full Debug CTest.
