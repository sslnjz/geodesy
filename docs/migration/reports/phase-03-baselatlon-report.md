# Phase Report

## Phase Name

Phase 3: Base LatLon

## Objective

- Stabilize the shared `LatLon` value type before spherical and ellipsoidal algorithm phases depend on it.
- Keep public latitude and longitude values in decimal degrees.
- Preserve JavaScript-style numeric and string construction through existing constructors and parse overloads.
- Verify formatting, GeoJSON output, equality, invalid inputs, and anti-meridian normalization.

## Modified Files

- `src/latlon/latlon.h`
  - Added direct standard-library includes required by the public header.
  - Documented `LatLon` as a decimal-degree value type.
  - Routed latitude and longitude setters through central validation and normalization helpers.
- `src/latlon/latlon.cpp`
  - Added finite-value validation for numeric and parsed coordinates.
  - Normalized stored longitude values to the half-open interval `[-180, 180)`.
  - Preserved established formatted anti-meridian output as `180E` while keeping stored longitude normalized.
  - Replaced machine-epsilon equality with an explicit coordinate tolerance for normalized floating-point values.
- `test/latlon_unittest.cpp`
  - Added coverage for DMS string construction, non-finite numeric input, malformed comma-separated parse input, setter failures, anti-meridian normalization, and normalized equality.
- `docs/migration/phase_status.md`
  - Marked Phase 3 completed and pointed to this report.

## Added Or Modified Tests

- Extended `constructor_with_strings` with DMS latitude and longitude input.
- Extended `constructor_fail` with NaN and infinity cases.
- Extended `parse_fail` with malformed comma-separated point input.
- Extended `setters_fail` with infinity cases.
- Added `normalizes_anti_meridian_longitudes`.
- Extended equality coverage for normalized longitudes and inequality.

## Commands Executed

- `git -c safe.directory=D:/workspace/github/geodesy status --short --branch`
- `cmake --build build --parallel --target latlon_unittest`
- `ctest --test-dir build -C Debug -R latlon_unittest --output-on-failure`
- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`
- `cmake --build build --parallel`
- `ctest --test-dir build --output-on-failure`
- `ctest --test-dir build -C Debug --output-on-failure`
- `ctest --test-dir build -C Debug -R "latlon_unittest|latlon_nvector_spherical_unittest" --output-on-failure`

## Test Results

- Focused LatLon build target: passed after allowing the Windows toolchain to read local SDK metadata.
- Focused LatLon CTest with `-C Debug`: passed, 1/1.
- Required Debug configure command: passed after allowing the Windows toolchain to read local SDK metadata.
- Required Debug build command: passed. One intermediate full-build attempt failed with `LNK1201` while writing `dms_unittest.pdb`; retrying the same command succeeded.
- Required raw CTest command `ctest --test-dir build --output-on-failure`: failed under the Visual Studio multi-config generator because no configuration was selected.
- Full Debug CTest with `-C Debug`: ran 10 tests, 8 passed and 2 failed.
  - `latlon_spherical_unittest.intersection`: existing expected string contains a trailing apostrophe not present in the actual result.
  - `latlon_ellipsoidal_vincenty_unittest.example`: existing expected strings contain non-UTF-8 degree-byte text while actual output uses UTF-8 degree text.
- Focused regression check for `latlon_unittest|latlon_nvector_spherical_unittest`: passed, 2/2.

## Unfinished Items

- Full raw CTest still needs `-C Debug` with the current Visual Studio multi-config build.
- Existing spherical and Vincenty baseline failures remain outside Phase 3 scope.
- The repository still contains earlier Phase 0/1/2 and documentation changes in the working tree; this phase did not revert or modify them.

## Next Phase Handoff

- Start Phase 4: Spherical LatLon.
- Read `AGENTS.md`, `docs/CODING_RULES.md`, the migration plan, `docs/migration/phase_status.md`, and this report before editing.
- Keep Phase 4 focused on `src/latlon/latlon_spherical.h`, `src/latlon/latlon_spherical.cpp`, and `test/latlon_spherical_unittest.cpp`.
- Preserve the completed base `LatLon` numeric behavior: public degrees, finite input rejection, normalized stored longitude, and established anti-meridian formatting.
- Address spherical distance, bearings, destinations, intersections, rhumb-line operations, polygon area, and boundary cases under focused tests.
