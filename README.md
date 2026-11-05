# trfrs — Tandem Repeats Finder (Rust port)

[![CI](https://github.com/keithyin/trfrs/actions/workflows/ci.yml/badge.svg)](https://github.com/keithyin/trfrs/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/keithyin/trfrs)](https://github.com/keithyin/trfrs/releases)
[![crates.io](https://img.shields.io/crates/v/trfrs)](https://crates.io/crates/trfrs)
[![License](https://img.shields.io/badge/license-AGPL--3.0--or--later-blue)](#许可证)
[![Rust](https://img.shields.io/badge/rust-2021-dea584?logo=rust)](https://www.rust-lang.org/)
[![Platforms](https://img.shields.io/badge/platforms-linux%20%7C%20macos%20%7C%20windows-blue)](#环境要求)

`trfrs` 是 Tandem Repeats Finder (TRF) 参考 C 实现 4.10.0 的 Rust 移植版本。
搜索算法与报告内容与 C 版本一致；命令行接口与输出文件布局按 Rust/Unix 惯例重新设计。

- 依赖只有 `clap`（命令行解析）、`thiserror`（错误类型派生）与 `libc`（Unix 下恢复 `SIGPIPE` 默认处理），产物为单个可执行文件
- 支持 FASTA / FASTQ 单序列/多序列文件、标准输入管道（FASTQ 为本移植扩展）
- 常规 GNU 风格选项（`--ngs`、`--html`）、类型化参数解析与范围校验、标准 `--help`/`--version`；退出码：`0` 成功、`1` 运行错误、`2` 命令行错误

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

`FILE` 是唯一必需的位置参数，搜索参数都是带默认值的 key-value 选项：

```bash
trfrs FILE [--match N] [--mismatch N] [--indel N] [--pm N] [--pi N] [--minscore N] [--minperiod N] [--maxperiod N] [OPTIONS]
```

### 搜索参数

| 参数 | 含义 | 默认值 |
| --- | --- | --- |
| `FILE` | 输入序列文件（FASTA 或 FASTQ，按每条记录的首字符自动识别）；`-` 表示从标准输入读取（需配合 `--ngs`） | 必需 |
| `--match` | 匹配得分（正整数） | `2` |
| `--mismatch` | 错配罚分（正整数） | `5` |
| `--indel` | 插入/缺失罚分（正整数） | `7` |
| `--pm` | 匹配概率（仅支持 75 与 80） | `80` |
| `--pi` | 插入/缺失概率（0–100） | `10` |
| `--minscore` | 报告重复所需的最小比对分数 | `50` |
| `--minperiod` | 最小重复单元长度，1 到 maxperiod | `1` |
| `--maxperiod` | 最大重复单元长度，1–2000 | `2000` |
| `--str` | STR 模式：等价于 `--minperiod 1 --maxperiod 15`，显式区间参数优先 | — |

默认值即常用的人类 DNA 参数组合 `2 5 7 80 10 50 2000`，因此只给文件名即可运行。

如果检测 repeats 少的，比如 unit2 * 2. 那么要调整 --minscore 为 8 。

### 输出选项

| 选项 | 说明 |
| --- | --- |
| `-o, --output <DIR>` | 输出文件写入的目录（不存在时自动创建，默认当前目录） |
| `-d, --dat` | 输出 `.dat` 数据文件 |
| `--csv` | 输出 `.csv` 表格（本移植新增，与 `.dat` 字段一一对应） |
| `-m, --mask` | 输出掩码序列文件（`.mask`） |
| `-f, --flank` | 报告中包含重复两侧的侧翼序列（默认 500 bp） |
| `--html` | 输出 HTML 报告（重复表格与比对明细；默认关闭） |
| `--no-redundancy` | 关闭冗余重复剔除 |
| `-n, --ngs` | 紧凑单行记录直接输出到 stdout，不写文件（适合流程管道） |
| `-v, --verbose` | 打印进度消息（默认关闭） |

### 其他选项

| 选项 | 说明 |
| --- | --- |
| `-l, --max-length <N>` | 预期的最大 TR 长度，单位百万（如 `-l 3` 表示 300 万；HG38 需 `-l 6`） |
| `-h, --help` | 打印帮助 |
| `-V, --version` | 打印版本 |

参数值在解析时即做类型与范围检查（如 `--maxperiod 0`、`--pm 79` 都是命令行错误，
退出码 2 并在 stderr 给出诊断）；`--minperiod` 不得超过 `--maxperiod`。
短选项可以合并（`-df` 等价于 `-d -f`），选项值支持 `--match 2` 与 `--match=2` 两种写法。

### 示例

```bash
# 使用仓库自带测试数据（全部参数取默认值；输出 HTML 报告）
./target/release/trfrs tests/data/small_test.fasta --html

# 显式给出全部搜索参数
./target/release/trfrs tests/data/small_test.fasta --html \
    --match 2 --mismatch 5 --indel 7 --pm 80 --pi 10 --minscore 50 --maxperiod 2000

# 同时生成 .dat 数据文件与 .csv 表格（便于 spreadsheet / pandas 直接读取）
./target/release/trfrs tests/data/multi.fasta -d --csv

# 输出写入指定目录
./target/release/trfrs tests/data/multi.fasta -d --csv -o results

# 紧凑单行格式直接输出到屏幕（适合流程管道）
./target/release/trfrs tests/data/small_test.fasta --ngs

# 从标准输入读取
cat tests/data/small_test.fasta | ./target/release/trfrs - --ngs

# 显示进度消息（默认安静）
./target/release/trfrs tests/data/small_test.fasta --ngs --verbose

# 版本信息 / 用法
./target/release/trfrs --version
./target/release/trfrs --help
```

也可以用 `cargo run` 运行，注意用 `--` 分隔 cargo 与程序参数：

```bash
cargo run --release -- tests/data/small_test.fasta --ngs
```

## 输入格式（FASTA / FASTQ）

`File` 可以是 FASTA（`>` 开头）或 FASTQ（`@` 开头）文件。格式按每条记录的首字符自动识别，
因此同一文件内混用两种记录也可以；`-` 标准输入同样适用。

```text
@read1 description
ACCCCTCAGGGTACCCCTCAGGGTACCCCTCAGGGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

FASTQ 支持多行 sequence 与多行 quality，quality 字符串按长度跳过、不参与打分，
序列按与 FASTA 完全相同的规则归一化（只保留 ASCII 字母并转大写）。输出中的序列描述
取自 `@` 之后的头部文本，与 FASTA 的 `>` 头部处理一致。

## STR 模式（period 1–15）

只关心短串联重复（STR）时，用 `--str` 把报告的重复单元长度限制在 1–15 bp（包含
period=1 的单碱基重复）：

```bash
trfrs sample.fasta --str --ngs
# 等价于
trfrs sample.fasta --minperiod 1 --maxperiod 15 --ngs
```

`--str` 只提供 period 区间的默认值，不改变打分参数；显式写出的
`--minperiod` / `--maxperiod` 优先生效：

```bash
trfrs sample.fasta --str --maxperiod 2000 --ngs   # 仍按 1–2000 报告
```

区间过滤发生在冗余剔除之前（与 C 版 `RemoveBySize()` 的时机相同），因此
30 bp 这类 period 为 15 倍数的重复会先被移除，15 bp 的 STR 不会因为
“冗余”而丢失。

使用该模式时需要注意：

- 默认 `--minscore 50` 要求约 25 个匹配碱基，`(AC)5` 这类短 STR 不会被报告。
  检测短 STR 通常要把阈值调低，例如 `--minscore 20`（更短的重复可用 8）：
  `trfrs sample.fasta --str --minscore 20 --ngs`。
- 默认开启的冗余剔除会去掉与更大周期重复重叠的较短 STR；需要完整目录时加 `--no-redundancy`。
- 该模式只做精确过滤，**不会**降低扫描量：内部 `max_distance` 仍有 500 的下限，
  因此 `--str` 的耗时与 `--maxperiod 500` 相当；打分矩阵仍约为 604 B/碱基，
  百万碱基级序列需要数百 MB，大基因组请配合 `-l` 估算或分块处理。

## 作为库使用

同一 crate 同时产出 `lib` 与 `bin`，库与命令行共用同一套算法，结果与 `.dat` 记录逐字段一致。
库调用完全在内存中完成：不写文件、不打印、不启动进程。

```toml
[dependencies]
trfrs = "4.10"
```

```rust
use trfrs::{Params, TrfFinder};

fn main() -> Result<(), trfrs::Error> {
    // 默认参数即 2 5 7 80 10 50 2000
    let finder = TrfFinder::new(&Params::default())?;

    let sequence = "TTTTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGG";
    for repeat in finder.find(sequence)? {
        // 区域为 0-based 半开区间，可直接切片
        println!(
            "[{}, {}) period={} copies={:.1} pattern={}",
            repeat.start, repeat.end, repeat.period, repeat.copies, repeat.pattern
        );
    }

    // 同一个 finder 可反复用于不同序列
    assert!(finder.find("ACGT")?.is_empty());
    Ok(())
}
```

### 参数

`Params` 的字段与命令行一一对应（`match_weight`、`mismatch`、`indel`、
`match_probability_percent`、`indel_probability_percent`、`minscore`、`min_period`、
`maxperiod`），另加 `remove_redundancy`（默认 `true`，等价于命令行未加 `-r`）。
只检测 STR 时用 `Params::str_preset()`，它等价于
`Params { min_period: 1, maxperiod: 15, ..Params::default() }`，与命令行 `--str` 语义一致：

```rust
use trfrs::{Params, TrfFinder};

let finder = TrfFinder::new(&Params::str_preset())?;
// 需要更短/更低拷贝数的 STR 时同时降低 minscore，与命令行一样：
let sensitive = TrfFinder::new(&Params {
    minscore: 20,
    ..Params::str_preset()
})?;
# Ok::<(), trfrs::Error>(())
```

`TrfFinder::new` 做与命令行相同的校验，非法值返回 `Error::InvalidParam`，不会退出进程：
`min_period` 必须为正且不大于 `maxperiod`，`maxperiod` 必须在 1–2000 之间；
`match_probability_percent` 仅支持 75 与 80（上游只提供这两种统计表）。

### 结果字段

| 字段 | 含义 |
| --- | --- |
| `start` / `end` | 0-based 半开区间，区域为 `&seq[start..end]`；也可用 `Repeat::sequence(&seq)` |
| `period` | 重复单元长度 |
| `copies` | 拷贝数（区域长度 / period） |
| `consensus_size` | consensus 长度，含 indel 时与 `period` 不同 |
| `score` | 比对分数 |
| `percent_matches` / `percent_indels` | 匹配 / indel 百分比（0–100 整数，与 `.dat` 列相同） |
| `percent_a` … `percent_t` | 重复区域碱基组成百分比（0–100 整数） |
| `entropy` | consensus pattern 的熵（0–2） |
| `pattern` | consensus pattern |

返回顺序与 `.dat` 一致（发现顺序）；需要按坐标排序时由调用方自行 `sort_by_key`。

### 输入与边界

- 接受 `&str`、`String`、`&[u8]`、`Vec<u8>` 等 `AsRef<[u8]>`；只保留 ASCII 字母并转大写，
  空白与换行被忽略，规则与 FASTA 读取一致。
- 非 `A/C/G/T`（如 `N`）按错配参与打分，与命令行行为一致。
- 归一化后为空序列返回 `Error::EmptySequence`。
- 打分权重（`match_weight`/`mismatch`/`indel`）与 `minscore` 必须落在 32 位有符号整数范围内，
  超出返回 `Error::InvalidParam`（算法内部以 `i32` 参与打分，不做静默截断）。
- 参数层面的错误始终以 `Err` 返回；内部算法不变式失败或大块内存分配失败为不可恢复错误，
  与 C 版本一致会终止进程，调用方应按 abort 对待。
- `min_period` 只做输出阶段的区间过滤，不影响打分，也不写入 `.dat` 头部与参数串；
  period 区间过滤在冗余剔除之前完成，因此 period 为其倍数的更长重复不会挤掉目标 STR。
- 内存：与 C 实现相同，打分矩阵约 604 字节/碱基，且与 `maxperiod` 无关；百万碱基级
  序列需要数百 MB；超大序列请自行分块，跨块的重复会被截断。`str_preset()` 只缩小报告
  范围，不会降低扫描量与内存占用。
- `TrfFinder` 只保存参数，自动 `Send + Sync`，可在多线程中共享 `&TrfFinder` 并发调用 `find`，
  每次调用独立分配扫描状态。

## 输出文件

输出文件默认写在当前工作目录，`--output/-o <DIR>` 可指定输出目录（不存在时自动创建）。
前缀取输入文件名（不含目录部分），文件名不包含搜索参数；
参数不同的运行会互相覆盖，需要保留多组结果时请分别使用不同的 `-o` 目录。

输出内容按需生成：`--html` 输出 HTML 报告，`-d`/`--csv`/`-m` 分别输出数据文件、
表格与掩码序列；一个都不给时只做搜索、不写任何文件（程序会给出提示）。

单序列输入：

| 文件 | 内容 |
| --- | --- |
| `<prefix>.1.html` | 每个重复的 HTML 注释页（`--html`） |
| `<prefix>.txt.1.html` | 比对（alignment）明细（`--html`） |
| `<prefix>.dat` | `-d` 时的数据文件 |
| `<prefix>.csv` | `--csv` 时的表格文件 |
| `<prefix>.mask` | `-m` 时的掩码序列 |

多序列输入：每个序列的文件带 `.s<序号>`（如 `<prefix>.s3.dat`），并额外生成
`<prefix>.summary.html` 汇总表（`--html`）；HTML 内的表格分页按 `<prefix>.<页号>.html` 编号。

`--csv` 的表格与 `.dat` 字段一一对应，首行为表头，每条重复一行；
多序列输入只写一个文件，用 `sequence_index` 列区分序列（编号与 `.s<序号>` 一致），
没有重复的序列不产生任何行：

| 列 | 含义 |
| --- | --- |
| `sequence_index` | 序列在输入文件中的序号（单序列恒为 1） |
| `sequence_name` | FASTA 头中的序列描述（含逗号或引号时按 RFC 4180 加引号） |
| `start` / `end` | 重复起止坐标（1-based，闭区间） |
| `period` | 重复单元长度 |
| `copy_number` | 拷贝数（1 位小数） |
| `pattern_size` | 共识模板长度 |
| `percent_matches` / `percent_indels` | 匹配 / indel 百分比 |
| `score` | 比对分数 |
| `percent_a` … `percent_t` | 重复区域内 A/C/G/T 的百分比 |
| `entropy` | 熵（2 位小数） |
| `repeat_unit` | 共识重复单元 |
| `repeat_sequence` | 重复区域的实际序列 |

`--ngs` 模式不写文件，只向 stdout 打印紧凑记录；此时不会再输出其他提示信息。
管道下游提前关闭时进程按 `SIGPIPE` 终止（shell 中退出状态 141），与 C 版本一致，不会产生 panic。

序列描述按 FASTA / FASTQ 头中的原始字节写入输出（即使不是合法 UTF-8）。
`.dat` 与 HTML 报告的内容（`Parameters:` 头、表格、比对明细）与 C 版本保持一致。

## 测试

`tests/ngs_output.rs` 使用 `tests/data/` 中的 FASTA 样本与 C 参考程序生成的期望输出做逐字节比对：

```bash
cargo test
```

覆盖 `--ngs` 输出、`.dat` 与 summary 输出、`--csv` 与 `.dat` 逐字段一致、短选项合并（`-df`）、
`-o` 输出目录、`--html`、非法参数与缺少参数值的报错、stdin 未配合 `--ngs` 的报错、
非 UTF-8 序列描述的字节保留、管道提前关闭时的 `SIGPIPE` 行为等场景。
`tests/data/small_test.fastq` 与同名 FASTA 夹具序列相同（含多行 sequence/quality），
验证 FASTQ 输入与 FASTA 输出逐字节一致。
`tests/data/str_short.fasta` 是专用于 STR 区间的夹具（period 1、2、3、5、15、16），
验证 `--str` 与 `--minperiod 1 --maxperiod 15`、`--maxperiod 15` 输出逐字节一致，
16 bp 单元被过滤、显式区间参数覆盖预设，以及非法 `--minperiod` 的报错与退出码。

`tests/lib_api.rs` 校验库接口：把 `TrfFinder::find` 的结果与 `tests/data/*.out`
（C 参考程序的 `-ngs` 输出）逐字段比对，并覆盖坐标换算、输入归一化、参数校验、
空序列、finder 复用与冗余剔除开关。`Params::str_preset()` 的结果同时与
命令行 `--str` 的记录逐字段比对，`min_period` 的过滤与校验也在其中。

## 代码结构

| 参考 C 源文件 | Rust 模块 |
| --- | --- |
| `trf.c` | `src/cli.rs`（命令行解析）、`src/driver.rs`（驱动与参数校验） |
| `trfrun.h` | `src/control.rs`（控制流程） |
| `tr30dat.h` / `tr30dat.c` | `src/state.rs`、`src/algo.rs`、`src/scan.rs` |
| `trfclean.h` | `src/clean.rs`（冗余剔除、HTML 输出） |
| FASTA / FASTQ 解析 | `src/fasta.rs` |
| 生成表 | `src/tables.rs` |
| `TRFPARAMSET` | `src/params.rs` |
| 库入口 | `src/lib.rs`（模块组织、公共导出）、`src/api.rs`（`TrfFinder`/`Params`/`Repeat`）、`src/main.rs`（二进制入口） |

C 版本使用文件级全局变量保存状态；Rust 版本将其收拢到 `state::Trf`，每个 FASTA 序列创建一个实例。

## 许可证

AGPL-3.0-or-later，与原 TRF 程序一致。详见各源文件头部声明。
