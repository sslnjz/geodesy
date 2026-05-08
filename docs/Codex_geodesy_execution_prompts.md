# Codex 执行 Geodesy JS -> C++17 迁移计划的多会话 Prompt 集合

> 适用项目：将 Chris Veness geodesy JavaScript 库能力迁移到本地 C++17 代码库。  
> 执行方式：以 `AGENTS.md` + `Coding_rules.md` + 迁移计划文档为固定上下文，按 Phase 单独开 Codex 会话执行。  
> 建议将迁移计划放入仓库：`docs/migration/2026-05-08-geodesy-js-to-cpp-migration.md`。

---

## 0. 总体执行原则

1. 不让 Codex 在一个会话里完成全部迁移。
2. 每个会话只处理一个 Phase，复杂 Phase 再拆成 A/B/C 子会话。
3. 每个会话必须先读：
   - `AGENTS.md`
   - `Coding_rules.md`
   - `docs/migration/2026-05-08-geodesy-js-to-cpp-migration.md`
   - `docs/migration/phase_status.md`，如果已经存在
   - 上一个 Phase 的报告，路径形如 `docs/migration/reports/phase-XX-*.md`
4. 每个会话必须在结束前写入本次交接摘要。
5. 任何未通过的测试必须记录，不允许只说“已完成”。
6. 不允许提前实现后续 Phase，除非当前 Phase 的编译必须做最小改动。
7. 遇到旧代码混乱时，以测试驱动的小范围重构为准，不做一次性大改。
8. 每个源码文件开头必须保留项目要求的 MIT license 注释块。
9. 核心代码需要自然、工程化注释；不得出现 AI、Codex、ChatGPT、generated、prompt 等痕迹。
10. 修改必须可验证：至少运行聚焦测试；条件允许时运行全量 `ctest`。


## 0.1 Superpowers 启动指令

你的迁移计划本身要求使用 `superpowers:subagent-driven-development` 或 `superpowers:executing-plans`。因此，每个正式执行会话的 Prompt 最前面建议增加以下一行：

```text
REQUIRED SUB-SKILL: Use superpowers:executing-plans for this phase. If the phase is explicitly split across independent review or exploration tasks, use superpowers:subagent-driven-development and consolidate all findings before editing code.
```

对于本项目，默认使用 `superpowers:executing-plans`；只有代码审查、跨模块探索、并行失败排查时才使用 `superpowers:subagent-driven-development`。

---

## 1. 仓库准备 Prompt

用于第一个 Codex 会话。目标是让 Codex 把计划和执行机制写入仓库，而不是开始迁移代码。

```text
你是本仓库的 C++17 迁移执行工程师。现在只做执行机制准备，不实现 geodesy 功能。

请先阅读：
- AGENTS.md
- Coding_rules.md
- docs/migration/2026-05-08-geodesy-js-to-cpp-migration.md，如果该文件不存在，请从当前任务附件或本会话提供内容创建该文件

任务：
1. 检查 AGENTS.md 是否明确要求迁移任务读取 Coding_rules.md 和 docs/migration/2026-05-08-geodesy-js-to-cpp-migration.md。
2. 如缺失，请最小修改 AGENTS.md，加入迁移执行协议：
   - 每次只执行一个 Phase。
   - 每个 Phase 结束必须写报告。
   - 不允许跳过测试。
   - 不允许留下 AI/Codex/ChatGPT/generated/prompt 等注释痕迹。
3. 创建 docs/migration/phase_status.md，用于记录每个 Phase 的状态。
4. 创建 docs/migration/reports/ 目录。
5. 创建 docs/migration/templates/phase_report_template.md，包含：
   - Phase 名称
   - 本次目标
   - 修改文件
   - 新增/修改测试
   - 执行命令
   - 测试结果
   - 未完成事项
   - 下一 Phase 交接说明
6. 不修改业务代码，不修改测试代码，除非为了创建迁移文档目录。

完成前请运行：
- git status --short

输出：
- 实际创建/修改的文件
- 后续应从 Phase 0 开始执行
```

---

## 2. Phase 执行通用 Prompt 模板

每个 Phase 单独开启一个 Codex 会话时使用。把 `{PHASE_NUMBER}`、`{PHASE_NAME}`、`{PHASE_SCOPE}` 替换成当前 Phase。

```text
你是本仓库的 C++17 迁移执行工程师。请严格按迁移计划执行当前 Phase，不要跳到后续 Phase。

必须先阅读以下文件并遵守：
- AGENTS.md
- Coding_rules.md
- docs/migration/2026-05-08-geodesy-js-to-cpp-migration.md
- docs/migration/phase_status.md，如果存在
- docs/migration/reports/ 中最近一个 Phase 报告，如果存在

当前任务：执行 Phase {PHASE_NUMBER}: {PHASE_NAME}

当前 Phase 范围：
{PHASE_SCOPE}

执行规则：
1. 先检查当前 git 状态，不覆盖用户未提交修改。
2. 只修改当前 Phase 计划列出的文件；如果必须修改其他文件以保证编译，请在报告中说明原因。
3. 先补充或修正测试，再实现或重构代码。
4. 保持 C++17，不引入非计划依赖。
5. 核心算法保留自然、工程化注释，但注释不得出现 AI、Codex、ChatGPT、generated、prompt 等痕迹。
6. 每个新增或修改的源码文件开头必须有指定 MIT license 注释块。
7. 公共角度输入输出使用 degrees，内部三角函数使用 radians，距离和高度使用 metres。
8. 使用 std::invalid_argument、std::domain_error、std::runtime_error 等标准异常，不抛裸指针，不使用非确定行为。
9. 不做与当前 Phase 无关的大规模格式化。

验证要求：
1. 运行当前 Phase 的聚焦测试。
2. 如构建允许，运行：
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
   cmake --build build --parallel
   ctest --test-dir build --output-on-failure
3. 如果某条命令失败，必须定位原因并尽力修复；如仍失败，报告中写清失败命令、失败原因、下一步建议。

收尾要求：
1. 更新 docs/migration/phase_status.md。
2. 写入 docs/migration/reports/phase-{PHASE_NUMBER}-{PHASE_NAME}-report.md。
3. 最后输出：
   - 完成内容
   - 修改文件
   - 测试命令和结果
   - 未完成事项
   - 下一 Phase 的建议启动 Prompt 摘要
```

---

## 3. Phase 0 Prompt：Build And API Baseline

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 0: Build And API Baseline。

必须先阅读：
- AGENTS.md
- Coding_rules.md
- docs/migration/2026-05-08-geodesy-js-to-cpp-migration.md

当前 Phase 目标：
建立可复现构建、测试与 JavaScript reference examples 的测试组织，不开始迁移算法。

允许检查/修改：
- CMakeLists.txt
- src/CMakeLists.txt
- include/geodesy/*.h
- test/CMakeLists.txt
- test/js_reference_examples_*_unittest.cpp
- docs/migration/phase_status.md
- docs/migration/reports/phase-00-build-api-baseline-report.md

任务：
1. 检查 Debug configure/build/ctest 是否可运行。
2. 建立 JavaScript reference examples 的测试组织，后续 Phase 可持续添加 golden cases。
3. 不做源代码大重构。
4. 记录当前失败测试和构建问题。
5. 更新 Phase 状态和报告。

必须运行：
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

---

## 4. Phase 1 Prompt：DMS And Angle Utilities

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 1: DMS And Angle Utilities。

先阅读：
- AGENTS.md
- Coding_rules.md
- docs/migration/2026-05-08-geodesy-js-to-cpp-migration.md
- docs/migration/phase_status.md
- docs/migration/reports/phase-00-build-api-baseline-report.md，如果存在

目标：
迁移/修复 dms.js 对应的 Dms 与 angle utility 行为。

允许修改：
- src/dms/dms.h
- src/dms/dms.cpp
- src/utils/algorithm.h
- test/dms_unittest.cpp
- test/js_reference_examples_*_unittest.cpp，如用于 DMS golden cases
- docs/migration/phase_status.md
- docs/migration/reports/phase-01-dms-angle-utilities-report.md

必须覆盖：
- parse
- toDms
- toLat
- toLon
- toBearing
- compassPoint
- wrap90
- wrap180
- wrap360

测试必须包含：
- decimal degrees
- DMS with hemisphere suffixes
- negative values
- unicode degree/prime/double-prime
- signed zero 相关输出
- wrapping
- compass precision

验证：
先运行聚焦 DMS 测试，再运行全量 Debug build + ctest。完成后更新状态和报告。
```

---

## 5. Phase 2 Prompt：Vector3d

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 2: Vector3d。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 1 报告。

允许修改：
- src/vector/vector3d.h
- src/vector/vector3d.cpp
- test/vector3d_unittest.cpp
- docs/migration/phase_status.md
- docs/migration/reports/phase-02-vector3d-report.md

任务：
1. 修正 public getters 和 operator declarations，确保 header 是可移植 C++17。
2. 实现或验证 plus、minus、scalar multiplication、scalar division、unary minus、dot、cross、unit、angle、rotation。
3. 明确处理 zero-length normalization。
4. 添加 unary minus、cross-product sign、zero vector、angle sign、rotation 测试。
5. 运行聚焦 vector 测试和全量 ctest。

不得提前迁移 Cartesian、n-vector 或 LatLon 算法。
```

---

## 6. Phase 3 Prompt：Base LatLon

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 3: Base LatLon。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 2 报告。

允许修改：
- src/latlon/latlon.h
- src/latlon/latlon.cpp
- test/latlon_unittest.cpp
- docs/migration/phase_status.md
- docs/migration/reports/phase-03-base-latlon-report.md

任务：
1. 将 base coordinate 明确为 value type，latitude/longitude 存储为 degrees。
2. 使用 explicit overloads 或 named factories 支持 JS 风格 parse/constructor 行为。
3. 实现或验证 toString、toGeoJSON、equality、longitude normalization。
4. 添加 numeric input、string input、invalid coordinates、anti-meridian longitudes、formatting 测试。
5. 保持 API 接近 JS 概念，但不要用宽泛 std::variant 构造器模拟 JS 动态输入。
6. 运行聚焦测试和全量 ctest。
```

---

## 7. Phase 4 拆分 Prompt：Spherical LatLon

Phase 4 较大，建议拆成 3 个会话。

### Phase 4A：基础球面路径

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 4A: Spherical LatLon basic geodesic methods。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 3 报告。

总 Phase 是 Phase 4，但本会话只做 4A。

允许修改：
- src/latlon/latlon_spherical.h
- src/latlon/latlon_spherical.cpp
- test/latlon_spherical_unittest.cpp
- test/js_reference_examples_*_unittest.cpp，如用于 spherical golden cases
- docs/migration/phase_status.md
- docs/migration/reports/phase-04a-spherical-basic-report.md

任务：
实现或验证：
- distanceTo
- initialBearingTo
- finalBearingTo
- midpointTo
- intermediatePointTo
- destinationPoint

测试：
- JavaScript reference examples
- poles
- anti-meridian crossing
- coincident points
- zero distance

不得处理 rhumb、cross-track、polygon area，除非已有代码需要最小修复以通过编译。
```

### Phase 4B：rhumb 与 track 方法

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 4B: Spherical LatLon rhumb and track methods。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 4A 报告。

允许修改：
- src/latlon/latlon_spherical.h
- src/latlon/latlon_spherical.cpp
- test/latlon_spherical_unittest.cpp
- docs/migration/phase_status.md
- docs/migration/reports/phase-04b-spherical-rhumb-track-report.md

任务：
实现或验证：
- rhumbDistanceTo
- rhumbBearingTo
- rhumbDestinationPoint
- rhumbMidpointTo
- crossTrackDistanceTo
- alongTrackDistanceTo
- intersection，如当前 API 中存在或计划中要求

测试：
- JS reference examples
- anti-meridian rhumb line
- coincident points
- path endpoint cases
- numerical tolerances

运行聚焦 spherical 测试和全量 ctest。
```

### Phase 4C：polygon area 与 Phase 4 收敛

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 4C: Spherical polygon area and Phase 4 convergence。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 4B 报告。

允许修改：
- src/latlon/latlon_spherical.h
- src/latlon/latlon_spherical.cpp
- test/latlon_spherical_unittest.cpp
- docs/migration/phase_status.md
- docs/migration/reports/phase-04c-spherical-polygon-convergence-report.md

任务：
1. 实现或验证 polygon area，包括 pole 和 anti-meridian handling。
2. 对 Phase 4A/4B 的 spherical 行为做一次收敛检查。
3. 不引入椭球、UTM、datum 功能。
4. 运行 spherical 聚焦测试和全量 ctest。
5. 在 phase_status.md 中把 Phase 4 标记为完成或部分完成，并说明原因。
```

---

## 8. Phase 5 拆分 Prompt：Ellipsoidal Coordinates And Cartesian

### Phase 5A：Ellipsoid 与 forward ECEF

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 5A: Ellipsoidal constants and geodetic to ECEF conversion。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 4C 报告。

允许修改：
- src/latlon/ellipsoids.h
- src/latlon/latlon_ellipsoidal.h
- src/latlon/latlon_ellipsoidal.cpp
- src/vector/cartesian.h
- src/vector/cartesian.cpp
- test/latlon_ellipsoidal_unittest.cpp
- docs/migration/phase_status.md
- docs/migration/reports/phase-05a-ellipsoidal-forward-ecef-report.md

任务：
1. 定义 Ellipsoid 常量，保留足够 reference precision。
2. 实现或验证 geodetic latitude/longitude/height to ECEF Cartesian。
3. 添加 WGS84、non-zero height、poles、equator 测试。
4. 不实现 datum transform。
```

### Phase 5B：reverse ECEF 与 round-trip

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 5B: ECEF Cartesian to geodetic and Phase 5 convergence。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 5A 报告。

允许修改：
- src/latlon/ellipsoids.h
- src/latlon/latlon_ellipsoidal.h
- src/latlon/latlon_ellipsoidal.cpp
- src/vector/cartesian.h
- src/vector/cartesian.cpp
- test/latlon_ellipsoidal_unittest.cpp
- docs/migration/phase_status.md
- docs/migration/reports/phase-05b-ellipsoidal-reverse-ecef-report.md

任务：
1. 实现或验证 Cartesian to geodetic conversion。
2. 添加 round-trip conversion 测试。
3. 检查 degrees/radians 边界与 height semantics。
4. 运行聚焦测试和全量 ctest。
```

---

## 9. Phase 6 Prompt：Datum And Helmert 7-Parameter Transforms

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 6: Datum And Helmert 7-Parameter Transforms。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 5B 报告。

允许修改：
- src/latlon/ellipsoids.h
- src/latlon/latlon_ellipsoidal_datum.h
- src/latlon/latlon_ellipsoidal_datum.cpp
- src/vector/cartesian_datum.h
- src/vector/cartesian_datum.cpp
- test/latlon_ellipsoidal_datum_unittest.cpp
- docs/migration/phase_status.md
- docs/migration/reports/phase-06-datum-helmert-report.md

任务：
1. 定义 Datum 和 Transform 为 immutable value-oriented data structures。
2. 实现或验证 WGS84、OSGB36、NAD83 等 JS reference datums。
3. 实现或验证 convertDatum，使用 Cartesian conversion 和 Helmert 7-parameter transform。
4. rotations arc-seconds 与 scale ppm 的单位转换必须只在计算边界发生一次。
5. 添加 WGS84 <-> OSGB36、identity datum、invalid datum path 测试。
6. 运行聚焦测试和全量 ctest。
```

---

## 10. Phase 7 拆分 Prompt：Vincenty Geodesics

### Phase 7A：Inverse

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 7A: Vincenty inverse solution。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 6 报告。

允许修改：
- src/latlon/latlon_ellipsoidal_vincenty.h
- src/latlon/latlon_ellipsoidal_vincenty.cpp
- test/latlon_ellipsoidal_vincenty_unittest.cpp
- docs/migration/phase_status.md
- docs/migration/reports/phase-07a-vincenty-inverse-report.md

任务：
1. 定义或修正 typed inverse result structs。
2. 实现或验证 Vincenty inverse solution。
3. 添加 distance、initial bearing、final bearing、coincident points、near-antipodal non-convergence 测试。
4. 非收敛不得返回 plausible-looking values；必须显式失败。
```

### Phase 7B：Direct 与收敛

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 7B: Vincenty direct solution and convergence。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 7A 报告。

允许修改：
- src/latlon/latlon_ellipsoidal_vincenty.h
- src/latlon/latlon_ellipsoidal_vincenty.cpp
- test/latlon_ellipsoidal_vincenty_unittest.cpp
- docs/migration/phase_status.md
- docs/migration/reports/phase-07b-vincenty-direct-report.md

任务：
1. 定义或修正 typed direct result structs。
2. 实现或验证 Vincenty direct solution。
3. 添加 destination、final bearing、round-trip、edge cases 测试。
4. 对 Phase 7A/7B 做一次聚焦收敛验证。
5. 运行聚焦测试和全量 ctest。
```

---

## 11. Phase 8 拆分 Prompt：UTM

### Phase 8A：UTM value type、parse、format

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 8A: UTM value type, parsing and formatting。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 7B 报告。

允许修改：
- src/utm/utm.h
- src/utm/utm.cpp
- test/CMakeLists.txt
- test/utm_unittest.cpp，如不存在则创建
- docs/migration/phase_status.md
- docs/migration/reports/phase-08a-utm-parse-format-report.md

任务：
1. 定义 Utm value type，包括 zone、Hemisphere、easting、northing、datum、convergence、scale。
2. 实现或验证 UTM parsing 和 formatting。
3. 添加 invalid zones、invalid easting、invalid northing 测试。
4. 不实现 lat/lon projection，除非已有 API 编译需要最小占位修复。
```

### Phase 8B：LatLon -> UTM

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 8B: latitude/longitude to UTM conversion。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 8A 报告。

允许修改：
- src/utm/utm.h
- src/utm/utm.cpp
- src/latlon/latlon_utm.h
- src/latlon/latlon_utm.cpp
- test/utm_unittest.cpp
- docs/migration/phase_status.md
- docs/migration/reports/phase-08b-latlon-to-utm-report.md

任务：
1. 实现或验证 latitude/longitude to UTM conversion。
2. 添加 ordinary zones、Norway exception、Svalbard exceptions、boundary latitude/longitude 测试。
3. 检查 Kruger series indexing。
4. 运行聚焦 UTM 测试。
```

### Phase 8C：UTM -> LatLon 与收敛

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 8C: UTM to latitude/longitude conversion and convergence。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 8B 报告。

允许修改：
- src/utm/utm.h
- src/utm/utm.cpp
- src/latlon/latlon_utm.h
- src/latlon/latlon_utm.cpp
- test/utm_unittest.cpp
- docs/migration/phase_status.md
- docs/migration/reports/phase-08c-utm-to-latlon-report.md

任务：
1. 实现或验证 UTM to latitude/longitude conversion。
2. 添加 round-trip conversion 测试。
3. 对 Phase 8A/8B/8C 做收敛验证。
4. 运行聚焦 UTM 测试和全量 ctest。
```

---

## 12. Phase 9 Prompt：MGRS

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 9: MGRS。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 8C 报告。

允许修改：
- src/utm/mgrs.h
- src/utm/mgrs.cpp
- src/utm/utm_mgrs.h
- src/utm/utm_mgrs.cpp
- src/latlon/latlon_utm_mgrs.h
- src/latlon/latlon_utm_mgrs.cpp
- test/CMakeLists.txt
- test/mgrs_unittest.cpp，如不存在则创建
- docs/migration/phase_status.md
- docs/migration/reports/phase-09-mgrs-report.md

任务：
1. 定义 MGRS grid zone、band、100 km square、easting、northing、precision fields。
2. 实现或验证 MGRS parsing 和 formatting。
3. 实现或验证 MGRS <-> UTM conversion。
4. 添加 precision levels、illegal grid letters、illegal bands、round-trip conversion 测试。
5. 注意 precision digits 必须按规则截断，不得错误四舍五入。
6. 运行聚焦 MGRS 测试和全量 ctest。
```

---

## 13. Phase 10 Prompt：OS Grid Reference

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 10: OS Grid Reference。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 9 报告。

允许修改：
- src/utm/osgridref.h
- src/utm/osgridref.cpp
- src/latlon/latlon_osgridref.h
- src/latlon/latlon_osgridref.cpp
- test/CMakeLists.txt
- test/osgridref_unittest.cpp，如不存在则创建
- docs/migration/phase_status.md
- docs/migration/reports/phase-10-osgridref-report.md

任务：
1. 定义 Airy 1830、true origin、false origin、scale 等 national grid constants。
2. 实现或验证 OS grid reference parsing 和 formatting。
3. 实现或验证 WGS84/OSGB36 conversion bridge through datum conversion。
4. 添加 multiple formatting precisions 和 malformed grid references 测试。
5. 运行聚焦 OS Grid 测试和全量 ctest。
```

---

## 14. Phase 11 拆分 Prompt：Reference Frames And Epoch Transforms

### Phase 11A：数据结构和 epoch 语义

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 11A: Reference frame data structures and decimal epoch semantics。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 10 报告。

允许修改：
- src/latlon/latlon_ellipsoidal_referenceframe_txparams.h
- src/latlon/latlon_ellipsoidal_referenceframe.h
- src/vector/cartesian_referenceFrame.h
- test/latlon_ellipsoidal_referenceframe_unittest.cpp
- docs/migration/phase_status.md
- docs/migration/reports/phase-11a-referenceframe-data-epoch-report.md

任务：
1. 将 string-only epoch arithmetic 替换为 decimal-year semantics。
2. 定义 14-parameter Helmert transform data，使用 named fields 或清晰的 std::array<double, 7> pairs。
3. 修复 optional misuse 和 pointer-throw 行为。
4. 添加 epoch decimal handling 和 metadata preservation 基础测试。
5. 不实现完整 conversion chains。
```

### Phase 11B：转换链和收敛

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 11B: Reference frame conversion chains and convergence。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 11A 报告。

允许修改：
- src/latlon/latlon_ellipsoidal_referenceframe_txparams.h
- src/latlon/latlon_ellipsoidal_referenceframe.h
- src/latlon/latlon_ellipsoidal_referenceframe.cpp
- src/vector/cartesian_referenceFrame.h
- src/vector/cartesian_referenceFrame.cpp
- test/latlon_ellipsoidal_referenceframe_unittest.cpp
- docs/migration/phase_status.md
- docs/migration/reports/phase-11b-referenceframe-conversion-report.md

任务：
1. 实现或验证 reference-frame conversion chains。
2. 添加 direct transform、chained transform、epoch rate adjustment、unsupported transforms、metadata preservation 测试。
3. unsupported path 必须 predictable fail。
4. 运行聚焦 reference-frame 测试和全量 ctest。
```

---

## 15. Phase 12 Prompt：Spherical N-Vector

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 12: Spherical N-Vector。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 11B 报告。

允许修改：
- src/vector/nvector_spherical.h
- src/vector/nvector_spherical.cpp
- src/latlon/latlon_nvector_spherical.h
- src/latlon/latlon_nvector_spherical.cpp
- test/latlon_nvector_spherical_unittest.cpp
- docs/migration/phase_status.md
- docs/migration/reports/phase-12-spherical-nvector-report.md

任务：
1. 实现或验证 lat/lon <-> n-vector conversion。
2. 实现或验证 great-circle、path intersection、nearest point、triangulation、trilateration、polygon enclosure、area、mean point operations。
3. 添加 poles、anti-meridian paths、segment endpoints、sign conventions 测试。
4. 注意 great-circle direction 和 cross-product sign。
5. 运行聚焦 n-vector spherical 测试和全量 ctest。
```

---

## 16. Phase 13 Prompt：Ellipsoidal N-Vector And NED

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 13: Ellipsoidal N-Vector And NED。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 12 报告。

允许修改：
- src/vector/nvector_ellipsoidal.h
- src/vector/nvector_ellipsoidal.cpp
- src/vector/nvector_cartesian.h
- src/vector/nvector_cartesian.cpp
- src/vector/ned.h
- src/vector/ned.cpp
- src/latlon/latlon_nvector_ellipsoidal.h
- src/latlon/latlon_nvector_ellipsoidal.cpp
- test/latlon_nvector_ellipsoidal_unittest.cpp
- docs/migration/phase_status.md
- docs/migration/reports/phase-13-ellipsoidal-nvector-ned-report.md

任务：
1. 实现或验证 ellipsoidal lat/lon <-> n-vector conversion。
2. 实现或验证 n-vector Cartesian conversion。
3. 实现或验证 NED delta 和 destination from delta。
4. 添加 height、datum、round-trip conversion、local NED movements 测试。
5. 明确 NED 是 local north/east/down metres，不是 ECEF x/y/z。
6. 运行聚焦测试和全量 ctest。
```

---

## 17. Phase 14 Prompt：API Compatibility Layer

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 14: API Compatibility Layer。

先阅读 AGENTS.md、Coding_rules.md、迁移计划、phase_status.md、Phase 13 报告。

允许修改：
- include/geodesy/*.h
- corresponding src/**/*.h
- public header compile tests 或 install-tree tests
- docs/migration/phase_status.md
- docs/migration/reports/phase-14-api-compatibility-report.md

任务：
1. 保留已经成为 C++ public API 的 JavaScript-style method names。
2. 对 parse 和 conversion entrypoints 添加 explicit named factories，避免 ambiguous constructors。
3. 统一 exception taxonomy：std::invalid_argument、std::domain_error、std::runtime_error。
4. 确保每个 public header 可以独立编译。
5. 不引入宽泛 std::variant constructor 去模拟所有 JS input shapes。
6. 运行 public header compile tests、全量 build 和 ctest。
```

---

## 18. Phase 15 Prompt：Full Library Convergence

```text
你是本仓库的 C++17 迁移执行工程师。请执行 Phase 15: Full Library Convergence。

先阅读：
- AGENTS.md
- Coding_rules.md
- docs/migration/2026-05-08-geodesy-js-to-cpp-migration.md
- docs/migration/phase_status.md
- docs/migration/reports/ 下所有 Phase 报告摘要

当前 Phase 目标：
不做新功能扩展，只做行为对齐、稳定性收敛、测试补齐和文档化剩余限制。

允许修改：
- 仅修改 failing tests 或 documented compatibility gaps 所需文件
- docs/migration/phase_status.md
- docs/migration/reports/phase-15-full-library-convergence-report.md
- 项目文档中用于记录 remaining limitations 的文件

任务：
1. 运行 Debug configure/build/ctest。
2. 运行 Release configure/build/ctest。
3. 对照 JavaScript reference examples 检查所有模块输出。
4. 移除或隔离 touched modules 中 fragile global mutable state。
5. 记录 remaining limitations。
6. 不混入新功能。

必须运行：
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure

最终输出：
- 全库完成度
- 未完成/限制项
- Debug 和 Release 测试结果
- 建议提交说明
```

---

## 19. 失败修复 Prompt

当某个 Phase 失败或 Codex 改动过大时使用。

```text
当前仓库在上一个 Codex 会话后存在失败。请不要继续后续 Phase，只做失败修复。

先阅读：
- AGENTS.md
- Coding_rules.md
- docs/migration/2026-05-08-geodesy-js-to-cpp-migration.md
- docs/migration/phase_status.md
- 最近的 docs/migration/reports/phase-*.md

任务：
1. 运行 git status --short，列出当前修改。
2. 运行上一报告中失败的命令，复现失败。
3. 定位失败到最小文件和最小函数。
4. 只修复导致失败的问题，不添加新功能，不进入下一个 Phase。
5. 如果上次改动范围明显超出 Phase，尽量回退或隔离无关改动。
6. 运行聚焦测试和必要的全量 ctest。
7. 更新上一 Phase 报告中的 failure/fix section，或创建 phase-XX-failure-repair-report.md。

输出：
- 根因
- 修复内容
- 验证命令
- 剩余风险
```

---

## 20. 代码审查 Prompt

每 2-3 个 Phase 后建议跑一次，或者在重要模块完成后跑。

```text
请对当前分支相对 main 的改动做代码审查。不要修改代码，除非发现明显编译错误且能做极小修复；否则只输出审查报告。

先阅读：
- AGENTS.md
- Coding_rules.md
- docs/migration/2026-05-08-geodesy-js-to-cpp-migration.md
- docs/migration/phase_status.md

审查重点：
1. 是否违反 Coding_rules.md。
2. 是否有 AI/Codex/ChatGPT/generated/prompt 等痕迹。
3. 每个源码文件是否有指定 MIT license 注释块。
4. public API 是否稳定，header 是否可独立编译。
5. degrees/radians、metres、height semantics 是否混淆。
6. datum transform、Vincenty、UTM 等高风险算法是否有测试覆盖。
7. 是否存在过度重构、跨 Phase 修改、隐藏全局状态。
8. 是否存在异常类型不规范、裸指针 throw、可疑 optional misuse。

输出：
- Critical issues
- High risk issues
- Medium risk issues
- Suggested fixes
- 是否建议继续下一个 Phase
```

---

## 21. 上下文压缩/交接 Prompt

当 Codex 会话上下文太长时，让它先写交接摘要，然后在新会话继续。

```text
当前会话上下文已经较长。请停止新增代码，只写一个可供下一 Codex 会话继续执行的交接摘要。

请创建或更新：
- docs/migration/reports/phase-XX-handoff.md

交接摘要必须包含：
1. 当前执行的 Phase 和子阶段。
2. 已完成的具体任务。
3. 已修改文件列表。
4. 新增/修改测试列表。
5. 已运行命令及结果。
6. 当前失败或未验证的命令。
7. 关键实现决策。
8. 不要重复犯的坑。
9. 下一会话应该读取的文件。
10. 下一会话应该使用的启动 Prompt。

不要继续写代码，不要开始新的 Phase。
```

---

## 22. 新会话恢复 Prompt

从交接摘要继续某个未完成 Phase。

```text
你是本仓库的 C++17 迁移执行工程师。请从上一个会话的交接摘要继续，不要重新开始，不要进入后续 Phase。

必须先阅读：
- AGENTS.md
- Coding_rules.md
- docs/migration/2026-05-08-geodesy-js-to-cpp-migration.md
- docs/migration/phase_status.md
- docs/migration/reports/phase-XX-handoff.md
- 当前 Phase 相关的最近报告

任务：
1. 复述当前 Phase 剩余任务，但不要输出冗长计划。
2. 检查 git status，确认已有改动。
3. 从 handoff 中的“下一步”继续执行。
4. 只修改当前 Phase 范围内文件。
5. 运行聚焦测试和必要的全量测试。
6. 更新 phase_status.md 和当前 Phase 报告。
```

---

## 23. 推荐执行批次

为了减少上下文污染，推荐按以下批次执行：

| 批次 | 会话 | 内容 |
|---|---|---|
| 0 | prep | 仓库准备，文档与状态文件 |
| 1 | phase-00 | Build/API baseline |
| 2 | phase-01 | DMS |
| 3 | phase-02 | Vector3d |
| 4 | phase-03 | Base LatLon |
| 5 | phase-04a | Spherical basic |
| 6 | phase-04b | Spherical rhumb/track |
| 7 | phase-04c | Spherical polygon/convergence |
| 8 | phase-05a | Ellipsoidal forward ECEF |
| 9 | phase-05b | ECEF reverse/round-trip |
| 10 | phase-06 | Datum/Helmert |
| 11 | phase-07a | Vincenty inverse |
| 12 | phase-07b | Vincenty direct/convergence |
| 13 | phase-08a | UTM parse/format |
| 14 | phase-08b | LatLon -> UTM |
| 15 | phase-08c | UTM -> LatLon/convergence |
| 16 | phase-09 | MGRS |
| 17 | phase-10 | OS Grid Reference |
| 18 | phase-11a | ReferenceFrame data/epoch |
| 19 | phase-11b | ReferenceFrame conversion |
| 20 | phase-12 | Spherical n-vector |
| 21 | phase-13 | Ellipsoidal n-vector/NED |
| 22 | phase-14 | API compatibility |
| 23 | phase-15 | Full convergence |

每完成 3-4 个批次，建议运行一次“代码审查 Prompt”。

---

## 24. 给 Codex 的最小上下文包

每次新会话只提供以下内容，不要复制整份迁移计划：

```text
上下文包：
1. 本仓库有 AGENTS.md 和 Coding_rules.md，必须先读。
2. 迁移计划在 docs/migration/2026-05-08-geodesy-js-to-cpp-migration.md。
3. 状态文件在 docs/migration/phase_status.md。
4. 报告目录在 docs/migration/reports/。
5. 当前只执行 Phase X / Phase XA，不执行后续阶段。
6. 结束前必须更新状态和报告。
```

如果 Codex 无法访问这些文件，再把当前 Phase 对应的 prompt 粘贴进去。

---

## 25. 建议提交策略

1. 每个 Phase 或子 Phase 一个 commit。
2. 提交信息不要出现 AI/Codex/ChatGPT/generated/prompt。
3. 推荐格式：
   - `migration: establish build and test baseline`
   - `migration: implement dms parsing and formatting`
   - `migration: stabilize vector3d operations`
   - `migration: add spherical distance and bearing support`
   - `migration: add vincenty inverse geodesics`
4. 如果一个 Phase 只修测试，提交信息可以是：
   - `test: add dms reference cases`
   - `test: cover utm zone edge cases`

