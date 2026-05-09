# Phase Report

## Phase Name

Phase 8C: UTM -> LatLon and convergence

## Objective

- Validate UTM reverse projection convergence and scale metadata behaviour.
- Cover a central-meridian UTM to latitude/longitude edge case where convergence is zero and scale is the UTM central-meridian scale.
- Reject invalid convergence and scale metadata before it can be stored on UTM value types or converted `LatLonUtm` results.

## Modified Files

- `src/utm/utm.cpp`
  - Tightened optional scale metadata validation to require a positive finite scale factor.
- `src/latlon/latlon_utm.cpp`
  - Added finite convergence validation to `LatLonUtm::setConvergence()`.
  - Added positive finite scale validation to `LatLonUtm::setScale()`.
- `test/utm_unittest.cpp`
  - Added central-meridian reverse-projection coverage.
  - Added invalid convergence and scale metadata regression coverage.
- `docs/migration/phase_status.md`
  - Advanced Phase 8 from Phase 8B completed to Phase 8C completed.

## Added Or Modified Tests

- `utm_unittest.rejects_invalid_value_type_inputs`
- `utm_unittest.utm_to_latlon_preserves_central_meridian_convergence_and_scale`
- `utm_unittest.rejects_invalid_latlon_projection_metadata`

## Commands Executed

- `git -c safe.directory=D:/workspace/github/geodesy status --short`
- `cmake --build build --target utm_unittest --config Debug --parallel`
- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`
- `cmake --build build --target utm_unittest --config Debug --clean-first --parallel`
- `ctest --test-dir build -C Debug -R utm_unittest --output-on-failure`
- `cmake --build build --parallel`
- `ctest --test-dir build --output-on-failure`
- `ctest --test-dir build -C Debug --output-on-failure`
- `git -c safe.directory=D:/workspace/github/geodesy restore -- build`

## Test Results

- Initial `git status` without command-scoped `safe.directory` failed due to this checkout's dubious-ownership guard; command-scoped `safe.directory` status succeeded.
- Configure/build attempts inside the sandbox needed local Windows SDK metadata access under `C:\Users\jbzs_\AppData\Local\Microsoft SDKs`; reruns with that access succeeded.
- Initial focused UTM test run after adding Phase 8C tests failed as expected, exposing missing projection-metadata validation. The run also showed stale object output from the existing build tree; a clean focused rebuild was required.
- Focused UTM clean build: `cmake --build build --target utm_unittest --config Debug --clean-first --parallel` passed.
- Focused UTM test: `ctest --test-dir build -C Debug -R utm_unittest --output-on-failure` passed, 1/1 tests.
- Debug configure: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug` passed.
- Debug build: `cmake --build build --parallel` passed.
- Required raw CTest: `ctest --test-dir build --output-on-failure` failed because the Visual Studio multi-configuration build requires `-C Debug`; all 11 tests were reported `Not Run`.
- Actual Debug CTest: `ctest --test-dir build -C Debug --output-on-failure` passed, 11/11 tests.

## Unfinished Items

- The required raw CTest command remains incompatible with the current Visual Studio multi-configuration build layout without `-C Debug`.
- MGRS conversion remains open for Phase 9.

## Next Phase Handoff

Continue with Phase 9: MGRS. Start from `docs/migration/phase_status.md`, this Phase 8C report, and the migration plan. Keep scope to `src/utm/mgrs.h`, `src/utm/mgrs.cpp`, `src/utm/utm_mgrs.h`, `src/utm/utm_mgrs.cpp`, `src/latlon/latlon_utm_mgrs.h`, `src/latlon/latlon_utm_mgrs.cpp`, and required MGRS tests unless compilation requires a documented exception. Add focused tests first for MGRS parsing, formatting precision, illegal grid letters, illegal bands, UTM bridge conversion, and round trips, then implement narrowly and run focused MGRS tests plus full Debug validation.
