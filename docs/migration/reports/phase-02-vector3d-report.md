# Phase Report

## Phase Name

Phase 2: Vector3d

## Objective

- Correct the `vector3d` public getters and operator declarations so the header is portable C++17.
- Verify and harden vector arithmetic, scalar multiplication and division, unary negation, dot product, cross product, unit vector, signed angle, and rotation.
- Add focused tests for zero-vector normalization, cross-product sign, signed angle direction, and degree-based rotation.

## Modified Files

- `src/vector/vector3d.h`
  - Removed non-standard extra class qualification from inline member declarations.
  - Fixed unary `operator-` to return the negated vector.
  - Routed scalar operators through the existing checked arithmetic helpers.
  - Added explicit invalid-scalar checks for multiplication and division.
  - Reused the computed cross product inside `angleTo` and documented the rotation degree-to-radian boundary.
- `src/vector/vector3d.cpp`
  - Added the required project license header.
  - Added direct standard-library includes used by the implementation.
  - Changed invalid vector construction to throw `std::invalid_argument`.
- `test/vector3d_unittest.cpp`
  - Added direct tests for zero-vector `unit()`, cross-product operand order, signed angle normal direction, degree-based rotation, and unary negation.
- `docs/migration/phase_status.md`
  - Marked Phase 2 completed and pointed to this report.

## Added Or Modified Tests

- Added `zero_vector_unit_is_noop`.
- Added `cross_product_preserves_operand_order`.
- Added `angle_sign_uses_supplied_plane_normal`.
- Added `rotation_uses_degree_input`.
- Extended `operators` coverage for unary negation.

## Commands Executed

- `git -c safe.directory=D:/workspace/github/geodesy status --short --branch`
- `cmake --build build --parallel --target vector3d_unittest`
- `ctest --test-dir build -C Debug -R vector3d_unittest --output-on-failure`
- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`
- `cmake --build build --parallel`
- `ctest --test-dir build --output-on-failure`
- `ctest --test-dir build -C Debug --output-on-failure`

## Test Results

- Initial focused Vector3d build required filesystem access to local Windows SDK metadata. The first retry reached link and failed with `LNK1168` because `build/test/Debug/vector3d_unittest.exe` was temporarily unavailable for writing; a later retry succeeded.
- Focused Vector3d build target: passed.
- Focused Vector3d CTest with `-C Debug`: passed, 1/1.
- Required Debug configure command: passed.
- Required Debug build command: passed.
- Required raw CTest command `ctest --test-dir build --output-on-failure`: failed under the Visual Studio multi-config generator because no configuration was selected.
- Full Debug CTest with `-C Debug`: ran 10 tests, 8 passed and 2 failed.
  - `latlon_spherical_unittest.intersection`: existing expected string contains a trailing apostrophe not present in the actual result.
  - `latlon_ellipsoidal_vincenty_unittest.example`: existing expected strings contain non-UTF-8 degree-byte text while actual output uses UTF-8 degree text.

## Unfinished Items

- Full raw CTest still needs `-C Debug` with the current Visual Studio multi-config build.
- Existing spherical and Vincenty baseline failures remain outside Phase 2 scope.
- The repository still contains earlier Phase 0/1 and documentation changes in the working tree; this phase did not revert or modify them.

## Next Phase Handoff

- Start Phase 3: Base LatLon.
- Read `AGENTS.md`, `docs/CODING_RULES.md`, the migration plan, `docs/migration/phase_status.md`, and this report before editing.
- Keep Phase 3 focused on `src/latlon/latlon.h`, `src/latlon/latlon.cpp`, and `test/latlon_unittest.cpp`.
- Add tests first for numeric input, string input, invalid coordinates, anti-meridian longitudes, and formatting.
- Preserve the current known full-CTest baseline failures unless the next task explicitly widens scope.
