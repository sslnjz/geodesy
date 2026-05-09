# Phase Report

## Phase Name

Phase 10: OS Grid Reference

## Objective

- Harden British National Grid parsing and precision formatting.
- Preserve OSGB36 projection formulae and WGS84/OSGB36 datum bridge behaviour.
- Cover malformed grid references, invalid easting/northing values, and reference examples with focused tests.

## Modified Files

- `src/utm/osgridref.h`
  - Added the required source header.
  - Made formatting const-correct and added easting/northing accessors.
  - Documented the value type as metres from the OS false origin.
- `src/utm/osgridref.cpp`
  - Reworked constructor validation to reject non-finite and out-of-range easting/northing with `std::invalid_argument`.
  - Fixed standard and compact grid parsing, including numeric comma-separated references and 100 km square scaling.
  - Fixed precision formatting and numeric formatting.
  - Corrected OSGB36 inverse projection coefficients by avoiding integer division, added explicit meridional-arc iteration tolerance and convergence failure handling.
- `src/latlon/latlon_osgridref.h`
  - Replaced dynamically allocated national-grid constants with value-oriented public structs.
  - Made `toOsGrid()` const.
- `src/latlon/latlon_osgridref.cpp`
  - Added the required source header.
  - Corrected forward projection coefficients by avoiding integer division.
  - Preserved the OSGB36-first datum bridge and millimetre rounding before `OsGridRef` construction.
- `test/CMakeLists.txt`
  - Added the focused `osgridref_unittest` target. This file is outside the listed source files but is required by the phase plan to enable OS Grid tests.
- `test/osgridref_unittest.cpp`
  - Added focused Phase 10 coverage.
- `docs/migration/phase_status.md`
  - Advanced Phase 10 to completed.

## Added Or Modified Tests

- `osgridref_unittest.formats_reference_example_at_multiple_precisions`
- `osgridref_unittest.parse_accepts_standard_compact_and_numeric_references`
- `osgridref_unittest.os_grid_to_latlon_matches_reference_in_wgs84_and_osgb36`
- `osgridref_unittest.latlon_to_os_grid_matches_reference_example`
- `osgridref_unittest.osgb36_latlon_to_os_grid_does_not_apply_second_datum_conversion`
- `osgridref_unittest.rejects_invalid_constructor_parse_and_format_inputs`

## Commands Executed

- `git status --short --branch`
- `git -c safe.directory=D:/workspace/github/geodesy status --short --branch`
- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`
- `cmake --build build --target osgridref_unittest --config Debug --parallel`
- `ctest --test-dir build -C Debug -R osgridref_unittest --output-on-failure`
- `cmake --build build --target utm mgrs_unittest --config Debug --clean-first --parallel`
- `cmake --build build --parallel`
- `ctest --test-dir build --output-on-failure`
- `ctest --test-dir build -C Debug --output-on-failure`
- `git -c safe.directory=D:/workspace/github/geodesy restore -- build`

## Test Results

- Initial plain `git status --short --branch` failed due to this checkout's dubious-ownership guard; command-scoped `safe.directory` status succeeded.
- Initial configure and build attempts inside the sandbox needed local Windows SDK metadata access under `C:\Users\jbzs_\AppData\Local\Microsoft SDKs`; reruns with that access succeeded.
- Initial focused OS Grid build failed before implementation due to non-const `toString()` and `toOsGrid()` APIs. Later focused CTest exposed overly strict DMS-derived reference tolerances; tests were adjusted to stay below the source example's 0.001-second precision.
- Focused OS Grid build: `cmake --build build --target osgridref_unittest --config Debug --parallel` passed.
- Focused OS Grid test: `ctest --test-dir build -C Debug -R osgridref_unittest --output-on-failure` passed, 1/1 tests.
- Debug configure: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug` passed.
- First full Debug build failed in `mgrs_unittest` with unresolved MGRS symbols from stale UTM/MGRS object files in the existing Phase 9 working tree. Clean rebuilding `utm` and `mgrs_unittest` fixed the stale-object issue without source changes.
- Debug build: `cmake --build build --parallel` passed.
- Required raw CTest: `ctest --test-dir build --output-on-failure` failed because the Visual Studio multi-configuration build requires `-C Debug`; all 13 tests were reported `Not Run`.
- Actual Debug CTest: `ctest --test-dir build -C Debug --output-on-failure` passed, 13/13 tests.

## Unfinished Items

- The required raw CTest command remains incompatible with the current Visual Studio multi-configuration build layout without `-C Debug`.
- Reference Frames and Epoch Transforms remain open for Phase 11.

## Next Phase Handoff

Continue with Phase 11: Reference Frames And Epoch Transforms. Start from `docs/migration/phase_status.md`, this Phase 10 report, and the migration plan. Keep scope to `src/latlon/latlon_ellipsoidal_referenceframe_txparams.h`, `src/latlon/latlon_ellipsoidal_referenceframe.h`, `src/latlon/latlon_ellipsoidal_referenceframe.cpp`, `src/vector/cartesian_referenceFrame.h`, `src/vector/cartesian_referenceFrame.cpp`, and `test/latlon_ellipsoidal_referenceframe_unittest.cpp` unless compilation requires a documented exception. Add focused tests first for decimal epoch handling, direct and chained transforms, rate-adjusted transform parameters, unsupported transform paths, and frame/epoch metadata preservation, then implement narrowly and run focused reference-frame tests plus full Debug validation.
