# Phase Report

## Phase Name

Phase 1: DMS And Angle Utilities

## Objective

- Refactor the DMS public interface around explicit numeric and string parsing entry points.
- Implement and harden DMS parsing, formatting, hemisphere suffixes, compass points, and angle wrapping.
- Move shared angle wrapping helpers into `src/utils/algorithm.h`.
- Add DMS boundary coverage for Unicode input, hemisphere suffixes, negative values, rounding, and compass precision.

## Modified Files

- `src/dms/dms.h`
  - Replaced the broad parse template with explicit `double`, `std::string`, and `const char*` overloads plus a narrow arithmetic overload.
  - Added `[[nodiscard]]` to value-returning DMS helpers.
  - Replaced the leaked static separator reference with an internal mutable separator accessor.
  - Retained prior transitive includes used by existing modules to avoid widening this phase into unrelated include repair.
- `src/dms/dms.cpp`
  - Moved string parsing implementation out of the public header.
  - Added complete-token numeric validation for DMS fields.
  - Fixed hemisphere suffix handling with trailing whitespace and lowercase `n/s/e/w`.
  - Reworked degree/minute/second rounding so carries are handled before formatting.
  - Kept the current C++ default separator as empty string to preserve existing test behaviour.
  - Delegated wrap helpers to `src/utils/algorithm.h`.
- `src/utils/algorithm.h`
  - Added `positiveModulo`, `wrap90`, `wrap180`, and `wrap360` helpers with negative-safe modulo behaviour.
  - Increased `pi` precision.
- `test/dms_unittest.cpp`
  - Added tests for Unicode DMS parsing, hemisphere suffixes, negative values, bearing rounding, and compass precision boundaries.
- `docs/migration/phase_status.md`
  - Marked Phase 1 completed and pointed to this report.

## Added Or Modified Tests

- Added `parse_unicode_hemisphere_and_negative_values`.
- Added `compass_precision_boundaries`.
- Added `format_boundary_values`.
- Updated DMS exception assertions to consume `[[nodiscard]]` return values cleanly.

## Commands Executed

- `git -c safe.directory=D:/workspace/github/geodesy status --short --branch`
- `cmake --build build --parallel --target dms_unittest`
- `ctest --test-dir build -C Debug -R dms_unittest --output-on-failure`
- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`
- `cmake --build build --parallel`
- `ctest --test-dir build --output-on-failure`
- `ctest --test-dir build -C Debug --output-on-failure`

## Test Results

- Focused DMS build target: passed.
- Focused DMS CTest with `-C Debug`: passed, 1/1.
- Required Debug configure command: passed after rerunning with filesystem access to local Windows SDK metadata.
- Required Debug build command: passed after rerunning with filesystem access to local Windows SDK metadata.
- Required raw CTest command `ctest --test-dir build --output-on-failure`: failed under the Visual Studio multi-config generator because no configuration was selected.
- Full Debug CTest with `-C Debug`: ran 10 tests, 8 passed and 2 failed.
  - `latlon_spherical_unittest.intersection`: existing expected string contains a trailing apostrophe not present in the actual result.
  - `latlon_ellipsoidal_vincenty_unittest.example`: existing expected strings contain non-UTF-8 degree-byte text while actual output uses UTF-8 degree text.

## Unfinished Items

- Root-level `Coding_rules.md` is absent; the applicable rules were read from `docs/CODING_RULES.md`.
- Full raw CTest still needs `-C Debug` with the current Visual Studio multi-config build.
- Existing spherical and Vincenty baseline failures remain outside Phase 1 scope.
- Some non-DMS headers still rely on transitive includes from `dms.h`; Phase 1 preserved that behaviour rather than changing unrelated module headers.

## Next Phase Handoff

- Start Phase 2: Vector3d.
- Read `AGENTS.md`, `docs/CODING_RULES.md`, the migration plan, `docs/migration/phase_status.md`, and this report before editing.
- Keep Phase 2 focused on `src/vector/vector3d.h`, `src/vector/vector3d.cpp`, and `test/vector3d_unittest.cpp`.
- Add tests first for unary minus, cross-product sign, zero-length normalization, signed angle behaviour, and rotation.
- Preserve the current known full-CTest baseline failures unless the next task explicitly widens scope.
