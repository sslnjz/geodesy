# AGENTS.md

## 1. Scope

This file applies to the entire repository unless a more specific `AGENTS.md` exists in a subdirectory.

The project is a C++17 implementation and refactoring of the JavaScript geodesy library originally maintained by Chris Veness. The target is not a loose rewrite: the C++ implementation must preserve the behaviour, numerical conventions, public concepts, and tested examples of the JavaScript library as closely as practical while using idiomatic C++17.

Primary module families to preserve:

- `dms`: degrees/minutes/seconds parsing, formatting, wrapping, compass points.
- `vector3d`: Cartesian vector operations used by n-vector and geocentric calculations.
- `latlon-spherical`: spherical earth distance, bearings, destination, midpoint, intersections, rhumb-line operations, polygon operations.
- `latlon-ellipsoidal`: geodetic latitude/longitude/height and geocentric Cartesian conversions.
- `latlon-ellipsoidal-datum`: historical datum definitions and Helmert transforms.
- `latlon-ellipsoidal-referenceframe`: dynamic reference frames and epoch-aware transforms.
- `latlon-ellipsoidal-vincenty`: Vincenty direct and inverse geodesic calculations.
- `latlon-nvector-spherical` and `latlon-nvector-ellipsoidal`: n-vector based geodesy.
- `utm` and `mgrs`: UTM and MGRS conversions and formatting.
- `osgridref`: British National Grid conversions.

## 2. Agent role

Act as a senior C++17 maintainer responsible for completing and hardening the port. Work from the current local C++ codebase, even if it is incomplete or inconsistent. Prefer targeted refactoring over wholesale replacement, but replace code that is numerically wrong, unsafe, or incompatible with the expected API.

Every change must be traceable to one of these goals:

1. Preserve or restore behaviour from the JavaScript reference implementation.
2. Improve C++17 correctness, type safety, build reliability, or test coverage.
3. Remove inconsistent, duplicated, or fragile code from the initial C++ version.
4. Clarify geodesy algorithms without changing their mathematical intent.

Do not introduce unrelated features, speculative APIs, or convenience wrappers that are not backed by the original library or by explicit project requirements.

## 3. Required workflow for each task

Before editing:

1. Inspect the current local implementation and identify the affected module, public API, and tests.
2. Compare the implementation with the corresponding JavaScript module and existing examples when available.
3. Identify whether the task is a bug fix, API alignment, refactor, numerical correction, or test addition.
4. Check `Coding_rules.md` and follow it for all code changes.
5. For migration work, read `docs/migration/2026-05-08-geodesy-js-to-cpp-migration.md` and `docs/migration/phase_status.md` before selecting the next phase.

While editing:

1. Keep changes minimal and module-local where possible.
2. Preserve existing public names unless there is a strong correctness reason to change them.
3. Add or update tests for every behavioural change.
4. Add comments to core geodesy calculations, constants, transformations, and edge-case handling.
5. Keep comments natural and technical. Do not mention tools, prompts, generated code, or automated assistance.
6. Use the required file header at the beginning of every new or modified project source file.

After editing:

1. Build the project.
2. Run all available tests.
3. Run focused tests for the changed module.
4. Check boundary cases involving poles, anti-meridian wrapping, zero distance, coincident points, near-antipodal points, invalid latitude/longitude, invalid UTM/MGRS zones, datum changes, and failed Vincenty convergence where relevant.
5. Summarise what changed, why it changed, what was tested, and any remaining limitation.

## 4. Build and test commands

Use the repository's existing build system first. If the repository uses CMake, these commands are the default baseline:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

For release checks:

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure
```

If the repository has custom scripts, package-manager targets, or CI commands, prefer those commands and record them in the task summary.

## 5. Source file header policy

Every project-owned source file must start with the following block exactly, before `#pragma once`, include guards, includes, namespace declarations, or any other content:

```cpp
/**********************************************************************************
*  MIT License                                                                    *
*                                                                                 *
*  Copyright (c) 2021 Binbin Song <ssln.jzs@gmail.com>                            *
*                                                                                 *
*  Geodesy tools for conversions between (historical) datums                      *
*  (c) Chris Veness 2005-2019                                                     *
*  www.movable-type.co.uk/scripts/latlong-convert-coords.html                     *
*  www.movable-type.co.uk/scripts/geodesy-library.html#latlon-ellipsoidal-datum   *
*                                                                                 *
*  Permission is hereby granted, free of charge, to any person obtaining a copy   *
*  of this software and associated documentation files (the "Software"), to deal  *
*  in the Software without restriction, including without limitation the rights   *
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      *
*  copies of the Software, and to permit persons to whom the Software is          *
*  furnished to do so, subject to the following conditions:                       *
*                                                                                 *
*  The above copyright notice and this permission notice shall be included in all *
*  copies or substantial portions of the Software.                                *
*                                                                                 *
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       *
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    *
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         *
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  *
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  *
*  SOFTWARE.                                                                      *
***********************************************************************************/
```

Apply this rule to `.h`, `.hpp`, `.hh`, `.cpp`, `.cc`, `.cxx`, and project-owned test source files. Do not add this project header to vendored third-party code.

## 6. Comment policy

Core code must be commented, especially where the implementation uses geodesy formulae, unit conversions, datum transformations, iteration, or non-obvious numerical tolerances.

Good comments explain:

- The formula or geodesy concept being implemented.
- Units of inputs and outputs, especially degrees vs radians and metres vs kilometres.
- Why a tolerance, wrap, clamp, or exceptional branch is necessary.
- The purpose of transformation parameters, ellipsoid constants, epochs, or reference frames.
- The reason for preserving a result compatible with the JavaScript reference implementation.

Avoid comments that:

- Mention AI, Codex, ChatGPT, generated output, prompts, or automated authorship.
- Explain trivial syntax.
- Apologise, speculate, or describe uncertainty.
- Use phrases such as "generated by", "AI converted", "the assistant", "as requested", or "TODO from prompt".
- Hide incorrect code behind vague comments.

Use concise, natural technical English for source comments unless an existing file consistently uses another language.

## 7. Numerical behaviour requirements

Geodesy calculations are sensitive to units, wrapping, and floating-point convergence. The following rules are mandatory:

1. Store angular public inputs and outputs in degrees unless the function name or documentation explicitly says radians.
2. Convert to radians only at the calculation boundary.
3. Normalise longitude consistently to the documented interval, normally `[-180, +180)` or the interval used by the original JavaScript function.
4. Normalise bearings to `[0, 360)` unless the original function specifies another range.
5. Preserve signed zero only when it affects documented formatting or test output.
6. Avoid integer truncation in zone, band, and grid calculations unless explicitly required by the algorithm.
7. Use explicit tolerances for floating-point comparisons; never compare non-integral geodesy results with raw `==` except for sentinel values.
8. Handle non-convergence in iterative algorithms explicitly. Do not return a plausible-looking value when the algorithm failed.
9. Retain ellipsoid and datum constants with sufficient precision to match the reference implementation.
10. Do not replace established geodesy formulae with approximations unless tests and documentation clearly justify the change.

## 8. API compatibility expectations

When designing or repairing C++ APIs:

1. Prefer value types for coordinates, vectors, datums, ellipsoids, and transform parameters.
2. Keep object names close to the source concepts: `Dms`, `Vector3d`, `LatLon`, `Cartesian`, `Datum`, `Ellipsoid`, `Transform`, `Utm`, `Mgrs`, `OsGridRef`.
3. Use `std::optional` for genuinely optional return values, not for ordinary error handling.
4. Use exceptions for invalid user inputs when the operation cannot produce a meaningful value.
5. Use `enum class` for closed sets such as format types, hemisphere values, or reference-frame identifiers.
6. Avoid macros except include guards or build-system integration where unavoidable.
7. Avoid hidden global mutable state. If a compatibility setting is required, document it and make it controlled.
8. Maintain clear namespace boundaries. Do not pollute the global namespace.
9. Do not expose implementation-only helpers in public headers unless required for tests or API compatibility.

## 9. Testing expectations

Tests must cover both ordinary examples and numerical edge cases. When a module is changed, add or update tests for that module before considering the task complete.

At minimum, preserve examples equivalent to the JavaScript test set for:

- DMS parsing and formatting.
- Spherical distance, bearing, midpoint, destination, intersection, rhumb-line, and polygon operations.
- Ellipsoidal Cartesian conversion and datum conversion.
- Reference-frame and epoch conversion.
- Vincenty direct and inverse calculations, including convergence failure cases.
- N-vector conversion and path operations.
- UTM conversion, including Norway/Svalbard zone exceptions.
- MGRS parsing, formatting, and precision handling.
- OS grid reference conversions.
- 3D vector operations.

Use tolerances appropriate to the result being checked. Tests should compare physical quantities in meaningful units, for example metres for distance and degrees or arc-seconds for angular values.

## 10. Refactoring constraints

Allowed:

- Split large files when it improves maintainability.
- Add small internal helpers for angle wrapping, radian conversion, tolerance checks, or string parsing.
- Replace duplicated formulas with well-named shared helpers.
- Move constants into `constexpr` definitions.
- Improve exception messages and input validation.
- Introduce tests that lock down corrected behaviour.

Not allowed without explicit task scope:

- Changing the library's conceptual model.
- Replacing C++17 with a newer language standard.
- Adding heavy third-party dependencies.
- Introducing runtime network access.
- Changing public formatting output without corresponding tests and documented rationale.
- Removing modules because they are incomplete or difficult.

## 11. Task summary format

For every completed task, report:

```text
Changed:
- <files and high-level changes>

Reason:
- <behaviour fixed, API aligned, or refactor goal>

Validation:
- <build command>
- <test command>
- <important focused tests>

Notes:
- <known limitation, follow-up, or "None">
```

Keep summaries factual. Do not mention tools or generated authorship in comments, source files, or commit messages.

## 12. Migration execution protocol

Migration work must follow the phase plan in `docs/migration/2026-05-08-geodesy-js-to-cpp-migration.md`.

1. Execute only one Phase per task unless the user explicitly widens the scope.
2. Start from the next open Phase recorded in `docs/migration/phase_status.md`.
3. At the end of every Phase, write a report under `docs/migration/reports/` using `docs/migration/templates/phase_report_template.md`.
4. Do not skip the required build, focused tests, or full test run. If a command cannot be run, record the exact reason in the Phase report and task summary.
5. Do not leave comments, documentation text, commit messages, or source text containing AI, Codex, ChatGPT, generated, prompt, or similar authorship traces.
