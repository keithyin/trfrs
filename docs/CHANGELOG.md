# Changelog

本项目遵循 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/) 风格，
版本号遵循 [Semantic Versioning](https://semver.org/lang/zh-CN/)。

## [Unreleased]

## [5.3.0] - 2026-11-20

### Changed

- `--ssr`：重叠消解从「任何重叠都整条丢弃」改为「**实质重叠**（交集 ≥ 较小 span 的
  一半）淘汰低优先级 tract + 仅**边界相接**时把后一条按整周期右移裁剪」。裁剪后不足
  3 个完整拷贝的 tract 整条丢弃。这避免了 TRF 常见的 1–2 bp 边界漂移误删紧邻的真实
  SSR 位点；同 cyclic word 的嵌套语义与「完整拷贝数多者优先、起点靠前者优先」的优先级
  保持不变。（`src/ssr.rs`）
- `scripts/ssr_check.py`：level 3 的 `interference` 只按实质重叠判定；参考 stretch 分类
  新增「边界裁剪」预期类别（含裁剪后整条丢弃），并在命令行与 HTML 报告中展示。
- 文档：同步更新 `README.md`、`docs/resolve-candidates-tutorial.html`（四步流程、
  50% 分界、整周期裁剪公式与走查、测试表、FAQ）与
  `docs/ssr-containment-resolution-bug.md`（历史说明）。

## [5.2.0] - 2026-11-20

### Added

- 新增 `--ssr` 模式与 `Params::ssr_preset()`：把 STR 报告改写为干净的
  `motif × 完整拷贝数` tract，支持 motif 最小周期归一化、相位重锚、同一区域多 tract
  提取与候选消解。（`src/ssr.rs`）
- 新增 `scripts/ssr_check.py`：`--ssr` 输出的独立参考实现、三层一致性校验与单文件
  HTML 报告。

### Docs

- 新增 `docs/resolve-candidates-tutorial.html` 候选消解教程。
- 新增 `docs/ssr-containment-resolution-bug.md` 跨 motif 包含消解的回归记录。
