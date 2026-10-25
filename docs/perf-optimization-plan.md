# trfrs 性能分析与优化方案

状态：待实施
日期：2026-10-25

## 1. 背景

trfrs 是 Tandem Repeats Finder (TRF) 4.10.0 C 实现的 Rust 移植，输出要求逐字节一致。
本文档记录 Rust 版与 C 版的实测性能差异、差距来源分析，以及对应的优化方案。

## 2. 实测性能对比

### 测试环境

- 256 核 Linux
- Rust 版：`cargo build --release`（`Cargo.toml` 已启用 fat LTO + `codegen-units = 1`）
- C 版：项目默认 `CFLAGS = -g -O2`（无 `-O3`、无 LTO）

### 测试 1：10 Mb 合成序列

数据：随机 ACGT 背景，每隔 500–2000 bp 以 30% 概率插入串联重复
（重复单元取自 `AT / CAG / TTAGGG / ACGTTACGTT / ATCATCATCATCGGT`，拷贝数 10–300）。

参数：`2 5 7 80 10 50 2000 -ngs`，各跑 3 轮：

| 版本 | 耗时（3 轮） | 相对 |
| --- | --- | --- |
| C (TRF 4.10.0) | 190.5 ± 0.1 s | 1.00x |
| Rust (trfrs) | 215.7 ± 0.3 s | 1.13x（慢约 13%） |

`-ngs` 输出逐字节一致（约 2.4 MB），正确性无损。

### 测试 2：2 Mb 纯随机序列（定位瓶颈）

| 版本 | 耗时 |
| --- | --- |
| C | 0.68 s |
| Rust | 0.82 s |

关键结论：tuple 扫描循环（`newtupbo`）只占总耗时约 3%，
约 97% 的时间花在候选触发的 wraparound DP（`narrowbandwrap`/`newwrap`）
及 `get_statistics` 上。因此总差距的 13% 主要来自 **DP 路径**，而非扫描循环。

## 3. 差距来源分析

DP 热路径是主要瓶颈，按影响排序：

1. **C 版 `narrowbandwrap` 用寄存器指针游走**
   （`TRF/src/tr30dat.c:182`，`register int *pup, *pdiag, *pcurr` 逐单元递减，
   循环体全部由宏展开；C 项目是单翻译单元，`-O2` 下全程内联）。
   Rust 版 `src/algo.rs:187-582` 用 `Vec<i32>` 索引访问 + bounds check，
   且经 `get_s`/`put_s` 间接（虽已 `#[inline]`）。
2. **Rust 热循环全部通过 `&mut self` 访问大结构体字段**（`src/state.rs` 的 `Trf`）：
   `self.s`、`self.diag`、`self.up` 之间可能被编译器判为潜在别名，
   热循环内重复 load；C 版是文件级全局数组，别名分析更自由。
3. **回溯循环里的 `Option<(i32, i32)>` + match**（`algo.rs:132-135, 1014` 等）
   替代 C 宏，最内层存在分支开销。
4. **`get_top_periods` 每个候选分配 2 个 Vec**（`scan.rs:66-67`），
   候选多时反复堆分配。
5. **`i32 → usize` 转换遍布内层循环**（`algo.rs`/`scan.rs`），
   部分 bounds check 可能未被 LLVM 消除。
6. 次要项：
   - FASTA 逐字节 `read(&mut [u8; 1])`（`fasta.rs:77-81`；C 版 `getc` 同样逐字节，差距不大）；
   - 多序列模式下 dat/mask/summary 用裸 `File` 无 `BufWriter`（`control.rs:422-462`）；
   - 输出内容用 `format!` 拼接大 `String`。

## 4. 优化方案

### 4.1 优先级 P0（预计将差距缩到 ~5% 以内，甚至反超）

#### 方案 A：`narrowbandwrap`/`newwrap` 内层循环改指针游走

- 位置：`src/algo.rs:187-582`（`narrowbandwrap`）、`src/algo.rs:585-800`（`newwrap`）。
- 做法：进入循环前把 `self.s`、`self.diag`、`self.up` 取为局部切片/原始指针，
  用 `get_unchecked`/裸指针完全复刻 C 的指针推进模式；
  消除内层循环的 bounds check 与 `self.` 字段重复 load。
- 风险：涉及 `unsafe`；需保持输出逐字节一致（用现有测试验证）。

#### 方案 B：`get_top_periods` 复用缓冲区

- 位置：`src/scan.rs:59-125`（调用点 `scan.rs:379`）。
- 做法：把每次调用分配的 `vec![0f64; len]` / `vec![0i32; len]`
  提升为 `Trf` 的可复用缓冲字段（`Vec::clear()` 后重填），消除每候选 2 次堆分配。
- 风险：低，纯机械改动。

### 4.2 优先级 P1

#### 方案 C：回溯循环用哨兵值替代 `Option`

- 位置：`src/algo.rs` 回溯相关（`backward_trace_cell`/`forward_trace_cell`
  及 `get_*_pair_alignment_with_copynumber`）。
- 做法：用 `-1` 之类哨兵值替代 `Option<(i32, i32)>`，去掉最内层 match 判空。

### 4.3 优先级 P2（顺手做的小项）

#### 方案 D：FASTA 读取改块读

- 位置：`src/fasta.rs:59-149`。
- 做法：`ByteReader` 内部改为块缓冲 + `memchr(b'>')` 找头，
  序列正文批量扫描后再逐字节分发，保持 `getc`/`ungetc` 语义不变。

#### 方案 E：多序列输出句柄包 `BufWriter`

- 位置：`src/control.rs:422-462`。
- 做法：summary/mask/dat 的 `File` 包一层 `BufWriter`，逐序列 `write_all` 后 flush。

#### 方案 F：构建参数（可选，注意公平性）

- 本机构建可加 `RUSTFLAGS="-C target-cpu=native"`；
  对比测试时 C 版需同加 `-march=native` 才公平。

## 5. 验证方法

1. 正确性：`cargo test`（`tests/ngs_output.rs` 逐字节比对 C 参考输出），
   另用 10 Mb 基准数据比对 `-ngs` 输出逐字节一致。
2. 性能：同一合成数据（10 Mb 含重复 + 2 Mb 纯随机），
   各跑 3 轮取中位数，对比 C 版 `-O2` 与 Rust 优化前后耗时。
3. 回归基线（2026-10-25 实测）：
   - 10 Mb 含重复：C 190.5 s，Rust 215.7 s（1.13x）
   - 2 Mb 纯随机：C 0.68 s，Rust 0.82 s
