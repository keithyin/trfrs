# trfrs — Tandem Repeats Finder (Rust port)

[![CI](https://github.com/keithyin/trfrs/actions/workflows/ci.yml/badge.svg)](https://github.com/keithyin/trfrs/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/keithyin/trfrs)](https://github.com/keithyin/trfrs/releases)
[![crates.io](https://img.shields.io/crates/v/trfrs)](https://crates.io/crates/trfrs)
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

`File` 是唯一必需的位置参数，7 个搜索参数都是带默认值的 key-value 选项：

```bash
trfrs File [--match n] [--mismatch n] [--indel n] [--pm n] [--pi n] [--minscore n] [--maxperiod n] [options]
```

| 参数 | 含义 | 默认值 |
| --- | --- | --- |
| `File` | 输入序列文件（FASTA 格式）；`-` 表示从标准输入读取（需配合 `-ngs`） | 必需 |
| `--match` | 匹配得分（正整数） | `2` |
| `--mismatch` | 错配罚分（正整数） | `5` |
| `--indel` | 插入/缺失罚分（正整数） | `7` |
| `--pm` | 匹配概率（整数） | `80` |
| `--pi` | 插入/缺失概率（整数） | `10` |
| `--minscore` | 报告重复所需的最小比对分数 | `50` |
| `--maxperiod` | 最大重复单元长度，1–2000 | `2000` |

默认值即常用的人类 DNA 参数组合 `2 5 7 80 10 50 2000`，因此只给文件名即可运行。

如果检测 repeats 少的，比如 unit2 * 2. 那么要调整 --minscore 为 8 。


### 示例

```bash
# 使用仓库自带测试数据（全部参数取默认值）
./target/release/trfrs tests/data/small_test.fasta

# 显式给出全部搜索参数（key-value 形式，也支持 --match=2 与 -match 2）
./target/release/trfrs tests/data/small_test.fasta \
    --match 2 --mismatch 5 --indel 7 --pm 80 --pi 10 --minscore 50 --maxperiod 2000

# 同时生成 .dat 数据文件
./target/release/trfrs tests/data/multi.fasta -d

# 紧凑 .dat 格式直接输出到屏幕（适合流程管道）
./target/release/trfrs tests/data/small_test.fasta -ngs

# 从标准输入读取
cat tests/data/small_test.fasta | ./target/release/trfrs - -ngs

# 版本信息 / 用法
./target/release/trfrs -v
./target/release/trfrs -u         # 无参数时也会打印用法
```

也可以用 `cargo run` 运行，注意用 `--` 分隔 cargo 与程序参数：

```bash
cargo run --release -- tests/data/small_test.fasta -ngs
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

长选项同样可用且接受唯一前缀缩写：`--dat`、`--mask`、`--flank`、`--html-off`、`--redund-off`、`--ngs`、`--maxlength`、`--match`、`--mismatch`、`--indel`、`--pm`、`--pi`、`--minscore`、`--maxperiod`、`--version`、`--help`。

短选项可以合并，选项值可以紧邻或使用 `=`；短选项字母不区分大小写（`-D` 与 `-d` 相同），
单短横线的长选项写法（`-ngs`、`-maxlength 3`）也照旧可用：

```bash
trfrs input.fasta -df                           # 等价于 -d -f
trfrs input.fasta -l3                           # 等价于 -l 3 或 -l=3
trfrs input.fasta -ngs                          # 等价于 --ngs
trfrs input.fasta -match 2 -minscore 50         # 等价于 --match 2 --minscore 50
```

命令行解析由 `clap` 完成；选项与文件名位置任意，
未知选项按 `clap` 的格式报错并以状态码 2 退出。

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
`match_probability_percent`、`indel_probability_percent`、`minscore`、`maxperiod`），
另加 `remove_redundancy`（默认 `true`，等价于命令行未加 `-r`）。
`TrfFinder::new` 做与命令行相同的校验，非法值返回 `Error::InvalidParam`，不会退出进程；
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
- 内存：与 C 实现相同，打分矩阵约 604 字节/碱基，百万碱基级序列需要数百 MB；
  超大序列请自行分块，跨块的重复会被截断。
- `TrfFinder` 只保存参数，自动 `Send + Sync`，可在多线程中共享 `&TrfFinder` 并发调用 `find`，
  每次调用独立分配扫描状态。

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

覆盖 `-ngs` 输出、`.dat` 与 summary 输出、短选项合并、缺少 `-l` 参数、
非 UTF-8 序列描述的字节保留、管道提前关闭时的 `SIGPIPE` 行为等场景。

`tests/lib_api.rs` 校验库接口：把 `TrfFinder::find` 的结果与 `tests/data/*.out`
（C 参考程序的 `-ngs` 输出）逐字段比对，并覆盖坐标换算、输入归一化、参数校验、
空序列、finder 复用与冗余剔除开关。

## 代码结构

| 参考 C 源文件 | Rust 模块 |
| --- | --- |
| `trf.c` | `src/cli.rs`（命令行解析）、`src/driver.rs`（驱动与参数校验） |
| `trfrun.h` | `src/control.rs`（控制流程） |
| `tr30dat.h` / `tr30dat.c` | `src/state.rs`、`src/algo.rs`、`src/scan.rs` |
| `trfclean.h` | `src/clean.rs`（冗余剔除、HTML 输出） |
| FASTA 解析 | `src/fasta.rs` |
| 生成表 | `src/tables.rs` |
| `TRFPARAMSET` | `src/params.rs` |
| 库入口 | `src/lib.rs`（模块组织、公共导出）、`src/api.rs`（`TrfFinder`/`Params`/`Repeat`）、`src/main.rs`（二进制入口） |

C 版本使用文件级全局变量保存状态；Rust 版本将其收拢到 `state::Trf`，每个 FASTA 序列创建一个实例。

## 许可证

AGPL-3.0-or-later，与原 TRF 程序一致。详见各源文件头部声明。
