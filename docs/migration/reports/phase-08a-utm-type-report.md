# Phase Report

## Phase Name

Phase 8A: UTM value type, parse, format

## Objective

- Harden the `Utm` value type around zone, hemisphere, easting, northing, datum, convergence, and scale.
- Preserve JavaScript-compatible UTM parse and format behaviour for four-field coordinates such as `31 N 448251 5411932`.
- Use explicit standard exceptions for invalid UTM fields.
- Add focused parsing and formatting tests before implementation changes.

## Modified Files

- `src/utm/utm.h`
  - Added value-type accessors for zone, hemisphere, datum, convergence, and scale.
  - Kept existing easting and northing accessors.
- `src/utm/utm.cpp`
  - Replaced the broken single-token parser with four whitespace-separated field parsing.
  - Added finite-value validation for easting, northing, optional convergence, and optional scale.
  - Corrected UTM easting and northing range checks.
  - Normalized invalid UTM inputs to `std::invalid_argument`.
  - Fixed `toString()` to pad the zone and include spaces between hemisphere, easting, and northing.
- `test/utm_unittest.cpp`
  - Added focused Phase 8A tests for field preservation, parse, lowercase hemisphere, formatting, invalid inputs, and `verifyEN=false`.
- `test/CMakeLists.txt`
  - Added `utm_unittest` to the CTest target list.
- `src/latlon/latlon_spherical.h`
  - Removed a duplicated include guard block that prevented any test target from compiling.
- `src/latlon/latlon_spherical.cpp`
  - Moved helper functions out of an accidental duplicate constructor body.
  - Removed a stale duplicate `areaOf()` block that attempted to mutate a `const` polygon parameter.
- `test/latlon_spherical_unittest.cpp`
  - Removed stray statements outside a test body that prevented full test compilation.

## Added Or Modified Tests

- `utm_unittest.value_type_preserves_fields_and_optional_projection_metadata`
- `utm_unittest.parse_accepts_space_separated_coordinate`
- `utm_unittest.parse_accepts_lowercase_hemisphere`
- `utm_unittest.to_string_formats_zone_hemisphere_easting_and_northing`
- `utm_unittest.to_string_pads_single_digit_zone`
- `utm_unittest.rejects_invalid_value_type_inputs`
- `utm_unittest.rejects_invalid_parse_inputs`
- `utm_unittest.verify_flag_allows_extended_coordinate_values`

## Commands Executed

- `git -c safe.directory=D:/workspace/github/geodesy status --short --branch`
- `cmake --build build --target utm_unittest --parallel`
- `ctest --test-dir build -C Debug -R utm_unittest --output-on-failure`
- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`
- `cmake --build build --parallel`
- `ctest --test-dir build --output-on-failure`
- `ctest --test-dir build -C Debug --output-on-failure`
- `git -c safe.directory=D:/workspace/github/geodesy restore -- build`

## Test Results

- Initial `git status` without command-scoped `safe.directory` failed due to this checkout's dubious-ownership guard; command-scoped `safe.directory` status succeeded.
- Initial focused build failed because CMake/MSBuild could not read local Windows SDK metadata from `C:\Users\jbzs_\AppData\Local\Microsoft SDKs`; rerun with local SDK metadata access passed after source fixes.
- Focused UTM build: `cmake --build build --target utm_unittest --parallel` passed.
- Focused UTM test: `ctest --test-dir build -C Debug -R utm_unittest --output-on-failure` passed, 1/1 tests.
- Debug configure: sandboxed `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug` failed on local Windows SDK metadata access; rerun with local SDK metadata access passed.
- Debug build: `cmake --build build --parallel` passed after two pre-existing spherical syntax defects were repaired as compile unblockers.
- Required raw CTest: `ctest --test-dir build --output-on-failure` failed because the Visual Studio multi-configuration build requires `-C Debug`; all 11 tests were reported `Not Run`.
- Actual Debug CTest: `ctest --test-dir build -C Debug --output-on-failure` passed, 11/11 tests.

## Unfinished Items

- Phase 8B should harden latitude/longitude to UTM conversion and UTM to latitude/longitude conversion.
- Norway and Svalbard zone exception tests remain for the projection-focused Phase 8 slice.
- The required raw CTest command remains incompatible with the current Visual Studio multi-configuration build layout without `-C Debug`.

## Next Phase Handoff

Execute the next Phase 8 slice for UTM forward and reverse projection. Start from `docs/migration/phase_status.md`, this Phase 8A report, and the migration plan. Keep scope to `src/utm/utm.h`, `src/utm/utm.cpp`, `src/latlon/latlon_utm.h`, `src/latlon/latlon_utm.cpp`, and required UTM tests unless compilation requires a documented exception. Add or extend focused tests first for ordinary latitude/longitude conversion, UTM reverse conversion, round trips, Norway and Svalbard zone exceptions, latitude limits, invalid zone overrides, convergence, scale, and anti-meridian-adjacent longitudes. Then fix the Kruger-series loops and projection validation narrowly, run focused UTM tests, Debug configure/build, raw CTest, and full Debug CTest with `-C Debug`.
