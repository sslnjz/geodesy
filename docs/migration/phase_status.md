# Geodesy Migration Phase Status

This file tracks execution status for the phased JavaScript-to-C++17 geodesy migration plan.

Source plan: `docs/migration/2026-05-08-geodesy-js-to-cpp-migration.md`

| Phase | Name | Status | Report | Notes |
|---|---|---|---|---|
| Phase 0 | Build And API Baseline | Completed with known baseline failures | `docs/migration/reports/phase-00-build-api-baseline-report.md` | Debug configure/build succeeded. Focused baseline reference test passed. Full Debug CTest runs with `-C Debug` but has existing spherical and Vincenty test failures. |
| Phase 1 | DMS And Angle Utilities | Completed with known baseline failures | `docs/migration/reports/phase-01-dms-report.md` | DMS focused test passed. Debug configure/build succeeded. Required raw CTest still needs `-C Debug`; full Debug CTest has existing spherical and Vincenty failures from the baseline. |
| Phase 2 | Vector3d | Completed with known baseline failures | `docs/migration/reports/phase-02-vector3d-report.md` | Focused Vector3d test passed. Debug configure/build succeeded. Required raw CTest still needs `-C Debug`; full Debug CTest has existing spherical and Vincenty failures from the baseline. |
| Phase 3 | Base LatLon | Completed with known baseline failures | `docs/migration/reports/phase-03-baselatlon-report.md` | Focused LatLon test passed. Debug configure/build succeeded. Required raw CTest still needs `-C Debug`; full Debug CTest has existing spherical and Vincenty failures from the baseline. |
| Phase 4 | Spherical LatLon | Completed with known baseline failures | `docs/migration/reports/phase-04c-spherical-polygon-report.md` | Phase 4A-4C completed spherical distance, bearing, path, rhumb-line, cross/along-track, intersection, and polygon-area hardening. Raw CTest still needs `-C Debug`; full Debug CTest still has the existing Vincenty encoding failure. |
| Phase 4B | Spherical rhumb and track operations | Completed with known baseline failures | `docs/migration/reports/phase-04b-spherical-rhumb-report.md` | Focused spherical test passed. Debug configure/build succeeded after allowing Windows SDK metadata access. Raw CTest still needs `-C Debug`; full Debug CTest still has the existing Vincenty encoding failure. |
| Phase 4C | Spherical polygon area | Completed with known baseline failures | `docs/migration/reports/phase-04c-spherical-polygon-report.md` | Focused spherical test passed. Debug configure/build succeeded after allowing Windows SDK metadata access. Raw CTest still needs `-C Debug`; full Debug CTest still has the existing Vincenty encoding failure. |
| Phase 5 | Ellipsoidal Coordinates And Cartesian | Not started | - | - |
| Phase 6 | Datum And Helmert 7-Parameter Transforms | Not started | - | - |
| Phase 7 | Vincenty Geodesics | Not started | - | - |
| Phase 8 | UTM | Not started | - | - |
| Phase 9 | MGRS | Not started | - | - |
| Phase 10 | OS Grid Reference | Not started | - | - |
| Phase 11 | Reference Frames And Epoch Transforms | Not started | - | - |
| Phase 12 | Spherical N-Vector | Not started | - | - |
| Phase 13 | Ellipsoidal N-Vector And NED | Not started | - | - |
| Phase 14 | API Compatibility Layer | Not started | - | - |
| Phase 15 | Full Library Convergence | Not started | - | - |
