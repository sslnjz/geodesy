# Phase Report

## Phase Name

Phase 0: Build And API Baseline

## Objective

- Confirm the current Debug build path.
- Make top-level CTest discovery usable for this repository.
- Add an initial JavaScript reference example test target for later migration phases.
- Record baseline build and test status before changing geodesy behavior.

## Modified Files

- `CMakeLists.txt`
  - Enabled top-level testing so `ctest --test-dir build` discovers tests from the root build directory.
- `test/CMakeLists.txt`
  - Removed stale installed-package lookup for in-tree tests.
  - Switched the GoogleTest tag from `release-1.13.0` to the available `v1.13.0` tag.
  - Added `js_reference_examples_baseline_unittest` to the test list.
  - Registered tests with target-aware CMake test commands.
  - Built tests as part of the default build.
  - Registered memory-check tests only when a memory-check tool is available.
- `test/js_reference_examples_baseline_unittest.cpp`
  - Added a minimal reference-example baseline target.
- `ui/CMakeLists.txt`
  - Removed stale installed-package lookup for the in-tree UI target.
  - This file is outside the Phase 0 inspect list, but the top-level build includes `ui`; without this change the required Debug configure could not complete.
- `docs/migration/phase_status.md`
  - Marked Phase 0 status and report path.

## Added Or Modified Tests

- Added `test/js_reference_examples_baseline_unittest.cpp`.
- Added the `js_reference_examples_baseline_unittest` CTest target.

## Commands Executed

- `git -c safe.directory=D:/workspace/github/geodesy status --short --branch`
- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`
- `cmake --build build --parallel --target js_reference_examples_baseline_unittest`
- `ctest --test-dir build -C Debug -R js_reference_examples_baseline_unittest --output-on-failure`
- `cmake --build build --parallel`
- `ctest --test-dir build --output-on-failure`
- `ctest --test-dir build -C Debug --output-on-failure`

## Test Results

- Debug configure: passed after in-tree CMake package lookup fixes.
- Focused build target: passed.
- Focused CTest with `-C Debug`: passed, 1/1.
- Full Debug build: passed.
- Required raw CTest command `ctest --test-dir build --output-on-failure`: failed under the Visual Studio multi-config generator because no configuration was selected.
- Full CTest with `-C Debug`: ran 10 tests, 8 passed and 2 failed.
  - `latlon_spherical_unittest.intersection`: expected string contains a trailing apostrophe not present in the actual result.
  - `latlon_ellipsoidal_vincenty_unittest.example`: expected strings contain non-UTF-8 degree-byte text while the actual output uses UTF-8 degree text.

## Unfinished Items

- `Coding_rules.md` is missing from the repository root and was not discoverable by filename search.
- Top-level CMake still emits policy/deprecation warnings around `CMakeLists_Headers.txt` and `add_custom_command(TARGET ... DEPENDS ...)`.
- Full CTest requires `-C Debug` with the Visual Studio multi-config generator.
- Existing spherical and Vincenty tests fail as described above. They were not fixed in Phase 0 because that would expand into module-specific behavior/test alignment work.

## Next Phase Handoff

- Start Phase 1: DMS And Angle Utilities.
- Read `AGENTS.md`, the migration plan, this status file, and this report before editing.
- Keep Phase 1 focused on `src/dms/dms.h`, `src/dms/dms.cpp`, `src/utils/algorithm.h`, and `test/dms_unittest.cpp` unless a build-only fix is required.
- Preserve the new JavaScript reference example test organization by adding DMS reference examples to a dedicated `test/js_reference_examples_*_unittest.cpp` target or extending the baseline organization in a narrow way.
