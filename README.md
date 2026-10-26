# trfrs — Tandem Repeats Finder (Rust port)

[![CI](https://github.com/keithyin/trfrs/actions/workflows/ci.yml/badge.svg)](https://github.com/keithyin/trfrs/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/keithyin/trfrs)](https://github.com/keithyin/trfrs/releases)
[![License](https://img.shields.io/badge/license-AGPL--3.0--or--later-blue)](#许可证)
[![Rust](https://img.shields.io/badge/rust-2021-dea584?logo=rust)](https://www.rust-lang.org/)
[![Platforms](https://img.shields.io/badge/platforms-linux%20%7C%20macos%20%7C%20windows-blue)](#环境要求)

`trfrs` 是 Tandem Repeats Finder (TRF) 参考 C 实现 4.10.0 的 Rust 移植版本。
命令行接口、控制流程与输出文件均与 C 版本保持一致，输出逐字节相同。

- 依赖只有 `clap`（命令行解析）与 `libc`（Unix 下恢复 `SIGPIPE` 默认处理），产物为单个可执行文件
- 支持 FASTA 单序列/多序列文件、标准输入管道
- 保留上游全部选项与退出码（含 `getopt_long_only` 风格的长选项与短选项合并）

原始 C 程序作者：Gary Benson（1999 年首发，2020 年最后更新），详见 <https://tandem.bu.edu/trf/trf.unix.help.html>。

## 环境要求

- Rust 工具链（支持 edition 2021，建议使用当前 stable）
- 无外部系统依赖

```bash
rustc --version
cargo --version
```

## 构建

```bash
# 调试构建
cargo build

# 发布构建（推荐；Cargo.toml 已启用 fat LTO，扫描性能接近 C 实现）
cargo build --release
```

产物为 `target/release/trfrs`。

## 运行

参数格式与 C 版完全一致，前 8 个参数是必需的：

```bash
trfrs File Match Mismatch Delta PM PI Minscore MaxPeriod [options]
```

| 参数 | 含义 |
| --- | --- |
| `File` | 输入序列文件（FASTA 格式）；`-` 表示从标准输入读取（需配合 `-ngs`） |
| `Match` | 匹配得分（正整数） |
| `Mismatch` | 错配罚分（正整数） |
| `Delta` | 插入/缺失罚分（正整数） |
| `PM` | 匹配概率（整数，通常 80） |
| `PI` | 插入/缺失概率（整数，通常 10） |
| `Minscore` | 报告重复所需的最小比对分数（通常 50） |
| `MaxPeriod` | 最大重复单元长度，1–2000（通常 2000） |

常用人类 DNA 参数组合：`2 5 7 80 10 50 2000`。

### 示例

```bash
# 使用仓库自带测试数据
./target/release/trfrs tests/data/small_test.fasta 2 5 7 80 10 50 2000

# 同时生成 .dat 数据文件
./target/release/trfrs tests/data/multi.fasta 2 5 7 80 10 50 2000 -d

# 紧凑 .dat 格式直接输出到屏幕（适合流程管道）
./target/release/trfrs tests/data/small_test.fasta 2 5 7 80 10 50 2000 -ngs

# 从标准输入读取
cat tests/data/small_test.fasta | ./target/release/trfrs - 2 5 7 80 10 50 2000 -ngs

# 版本信息 / 用法
./target/release/trfrs -v
./target/release/trfrs            # 参数不足时打印用法
```

也可以用 `cargo run` 运行，注意用 `--` 分隔 cargo 与程序参数：

```bash
cargo run --release -- tests/data/small_test.fasta 2 5 7 80 10 50 2000 -ngs
```

### 选项

| 选项 | 说明 |
| --- | --- |
| `-m` | 额外输出掩码序列文件（`.mask`） |
| `-f` | 输出重复两侧的侧翼序列（默认 500 bp） |
| `-d` | 输出 `.dat` 数据文件 |
| `-h` | 关闭 HTML 输出（同时隐含 `-d`） |
| `-r` | 关闭冗余重复剔除 |
| `-l <n>` | 预期的最大 TR 长度，单位百万（如 `-l 3` 表示 300 万；HG38 需 `-l 6`） |
| `-ngs` | 多序列文件使用紧凑 `.dat` 输出到 stdout，成功返回 0 |

长选项同样可用且接受唯一前缀缩写：`--dat`、`--mask`、`--flank`、`--html-off`、`--redund-off`、`--ngs`、`--maxlength`、`--version`、`--help`。

短选项可以合并，选项值可以紧邻或使用 `=`；短选项字母不区分大小写（`-D` 与 `-d` 相同），
单短横线的长选项写法（`-ngs`、`-maxlength 3`）也照旧可用：

```bash
trfrs input.fasta 2 5 7 80 10 50 2000 -df      # 等价于 -d -f
trfrs input.fasta 2 5 7 80 10 50 2000 -l3      # 等价于 -l 3 或 -l=3
trfrs input.fasta 2 5 7 80 10 50 2000 -ngs     # 等价于 --ngs
```

命令行解析由 `clap` 完成；选项写在 8 个位置参数之前或之后均可，
未知选项按 `clap` 的格式报错并以状态码 2 退出。

## 输出文件

输出文件写在当前工作目录，前缀取输入文件名（不含目录部分），
`<参数串>` 为 `Match.Mismatch.Delta.PM.PI.Minscore.MaxPeriod`。

单序列输入：

| 文件 | 内容 |
| --- | --- |
| `<prefix>.<参数串>.1.html` | 每个重复的 HTML 注释页 |
| `<prefix>.<参数串>.1.txt.html` | 比对（alignment）明细 |
| `<prefix>.<参数串>.dat` | `-d` 时的数据文件 |
| `<prefix>.<参数串>.mask` | `-m` 时的掩码序列 |

多序列输入：每个序列的 HTML 文件带 `.s<序号>`，并额外生成 `<prefix>.<参数串>.summary.html` 汇总表。

`-ngs` 模式不写文件，只向 stdout 打印 `.dat` 内容；此时不会再输出其他提示信息。
管道下游提前关闭时进程按 `SIGPIPE` 终止（shell 中退出状态 141），与 C 版本一致，不会产生 panic。

序列描述按 FASTA 头中的原始字节写入输出（即使不是合法 UTF-8），与 C 版本逐字节相同。

## 测试

`tests/ngs_output.rs` 使用 `tests/data/` 中的 FASTA 样本与 C 参考程序生成的期望输出做逐字节比对：

```bash
cargo test
```

包含 8 个用例，覆盖 `-ngs` 输出、`.dat` 与 summary 输出、短选项合并、缺少 `-l` 参数、
非 UTF-8 序列描述的字节保留、管道提前关闭时的 `SIGPIPE` 行为等场景。

## 代码结构

| 参考 C 源文件 | Rust 模块 |
| --- | --- |
| `trf.c` | `src/cli.rs`（命令行解析）、`src/main.rs`（入口与参数校验） |
| `trfrun.h` | `src/control.rs`（控制流程） |
| `tr30dat.h` / `tr30dat.c` | `src/state.rs`、`src/algo.rs`、`src/scan.rs` |
| `trfclean.h` | `src/clean.rs`（冗余剔除、HTML 输出） |
| FASTA 解析 | `src/fasta.rs` |
| 生成表 | `src/tables.rs` |
| `TRFPARAMSET` | `src/params.rs` |

C 版本使用文件级全局变量保存状态；Rust 版本将其收拢到 `state::Trf`，每个 FASTA 序列创建一个实例。

## 许可证

AGPL-3.0-or-later，与原 TRF 程序一致。详见各源文件头部声明。
