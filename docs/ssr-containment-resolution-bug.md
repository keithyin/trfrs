# SSR 提取 Bug 2：跨 motif 的包含消解错误——发现与处理

> 适用版本：trfrs 5.2.0（`--ssr` 首次引入）。文中行号以该版本为准，后续改动可能使行号偏移，
> 以函数名为准。
>
> **后续变更**：第 ③ 步已从"任何重叠都整条丢弃"改为"实质重叠（交集 ≥ 较小 span 的一半）
> 淘汰高优先级落败者 + 仅边界相接时按整周期裁剪后一条"。本文保留的是当时修复的语义与回归
> 数据；当前行为见 `docs/resolve-candidates-tutorial.html` 与 `src/ssr.rs` 中
> `resolve_candidates` / `trim_boundary_overlaps` 的注释。

## 摘要

`--ssr` 过滤提取的候选消解阶段，旧实现把"被更大 tract 完全包含"一律当作嵌套产物丢弃，
**不区分两者是否属于同一个 cyclic word（旋转等价的 motif）**。这导致精确 `ATATATATC×3`
区域内部真实存在的 `AT×4` 精细重复被错误吞掉（数据集 seq 55 丢失 4 条 tract）。
修复后：只有**同一 cyclic word** 的嵌套（如 `(ATG)2 ⊂ (GAT)3`）才按相位嵌套产物丢弃；
不同 motif 之间的包含与重叠交给贪心重叠消解（完整拷贝数多者优先、起点靠前者优先）处理。

## 一、处理链路：这段逻辑在程序的什么位置

`--ssr` 模式下，一条输入序列的完整处理链路与各处理点的代码位置如下：

```
TRF 扫描产出 STR 记录
        │
        ▼
src/clean.rs :: Trf::trf_clean()                       (clean.rs:511)
  ① remove_by_period()      —— period 区间过滤          (clean.rs:514)
  ② remove_redundancy()     —— TRF 既有冗余剔除          (clean.rs:516-518)
  ③ params.ssr 分支：把 STR 记录整体改写为干净 SSR tract   (clean.rs:519-530)
        │  调用 crate::ssr::replace_with_ssr_tracts()
        ▼
src/ssr.rs（本 bug 所在模块，整条过滤逻辑都在这里）
  ④ replace_with_ssr_tracts()                           (ssr.rs:67)
     入口纯函数：对每条存活记录调用 ⑤，收集候选 tract
  ⑤ extract_record_tracts()                             (ssr.rs:111)
     归一化 motif（minimal_period，ssr.rs:88）→ 枚举旋转锚定
     → 逐拷贝延伸 → 产出候选 Tract（本 bug 的"受害者"在此生成）
  ⑥ resolve_candidates()                                (ssr.rs:190)
     候选消解，三步：                                    (ssr.rs:190-220)
       a. 相同 span 去重                                 (ssr.rs:191-192)
       b. 包含丢弃 ← ★ 本 bug 的错误处理位置（修复前
          不区分 motif；现加 same_motif 条件）            (ssr.rs:194-205)
       c. 贪心重叠消解                                    (ssr.rs:207-218)
     same_motif() 判定助手                               (ssr.rs:103)
  ⑦ build_record()                                      (ssr.rs:223)
     每个存留 tract 构造 RepeatRecord（matches=100、
     score=(copies-1)×period×match 等）
        │
        ▼
写回 self.records，后续 .dat/.csv/.mask/HTML 输出全部自动跟随
```

Python 校验脚本（发现该 bug 的工具）中对应的处理位置：
`scripts/ssr_check.py` —— `reference_stretches()` 生成全序列穷举参考；
main() 的 level 3 循环逐 stretch 分类（`canonical tract missing` 即本 bug 的信号）。

## 二、现象：一致性校验如何发现它

### 信号

对测试数据 `str-seqs.fa`（70 条序列）跑三层校验时，level 3 报出若干
`canonical tract missing`：某个参考 stretch 被**同 motif 的 STR 区域完全覆盖**、
也没有更大 tract 干扰，但 `--ssr` 输出里没有与参考 tract 精确一致的记录。其中
seq 55 一组四条尤为整齐（1-based 坐标）：

```
参考 stretch（穷举精确扫描）          --ssr 输出
AT [496..503]  ×4                    （缺失）
AT [505..512]  ×4                    （缺失）
AT [522..529]  ×4                    （缺失）
AT [531..538]  ×4                    （缺失）
```

这四条都是 level 1 可独立验证的精确 `motif × copies` 区域，不可能凭空错。

### 定位过程

把 `src/ssr.rs` 的算法用 Python 逐行镜像（`extract_record_tracts` + `resolve_candidates`），
对 seq 55 分别打印**候选生成**与**消解结果**两层：

- 候选层：四条 AT tract **都在**（AT 记录的旋转扫描正常找到它们）；
- 消解层：四条全部消失——消失位置恰好是"包含丢弃"步骤。

由此把问题锁定在 `resolve_candidates` 的第②步，而不是提取本身。
（同一轮校验还报了 seq 24 的 TTC/TCT 差异，镜像比对证明那是贪心消解的**设计行为**——
候选都在、赢者按规则产生——不属于 bug。区分"候选没生成"与"候选被规则丢掉"是
这次定位的关键动作。）

### 修复前的错误行为（seq 55，1-based）

该位点的序列结构是 `ATATATATC ×3`（period 9），TRF 在同一位置报了 period 2（motif AT）
和 period 9（motif ATATATATC）两条记录：

```
495  503            504  512            521  529            530  538
A T A T A T A T  C  A T A T A T A T  C  A T A T A T A T  C ...
└── AT×4 ──┘        └── AT×4 ──┘        └── AT×4 ──┘
└──────────── ATATATATC×3（span [486..513]）────────────┘
└──────────── TATATATCA×3（span [514..540]，另一旋转）──┘
```

旧包含丢弃按（first 升序、last 降序）扫描并维护"已达右端 reach"：
period-9 候选 `[487..513]` 先入（reach=513），四条 AT tract 的右端全部 ≤ 513，
**整批被当成嵌套产物丢弃**；随后 `[514..540]` 同样吞掉 `[522..529]`、`[531..538]`。
最终输出只剩 period-9 与更右侧的 AT tract，四条 AT×4 凭空消失。

## 三、根因分析

### 原规则与它依赖的论证

旧的包含丢弃规则：

```rust
// 修复前（示意）
if kept.first <= tract.first && tract.last <= kept.last {
    // 当作嵌套产物丢弃 —— 不看 motif
}
```

当时的论证是"容器的完整拷贝数一定 ≥ 被包含者"，因此丢弃被包含者不损失信息。
这个论证依赖 **Fine–Wilf 周期定理**：一个长度足够的窗口同时具有周期 q 与 p 时，
其最小周期整除 gcd(q, p)。**该定理只对"被包含窗口与容器是同一个 cyclic word"
的情形给出想要的结论**（最小周期相同 ⇒ 容器 span 更长 ⇒ 拷贝数更多）。

### 反例：不同 cyclic word 的包含

9-mer `ATATATATC` 的 9 个碱基中前 8 个本身就是 AT-周期的，所以它的精确区域内
**恰好嵌着**一段与它词族完全不同的 `AT×4` 精确重复。此时：

- 被包含者 period 2、4 拷贝；容器 period 9、仅 3 拷贝——"容器拷贝数更多"不成立；
- 被包含者不是容器的相位平移子窗口，而是**真实存在的另一类 SSR**；
- 报告它符合 `--ssr` 的语义（干净 motif×repeats），丢弃它则是漏报。

换言之，"包含 ⇒ 嵌套产物"只在同词情形成立；跨词包含必须交给后续的
贪心重叠消解按统一规则裁决，而不能预先武断丢弃。

## 四、当前处理（修复后的语义与代码）

### 修复点 1：包含丢弃加 `same_motif` 条件

`src/ssr.rs :: resolve_candidates()`（ssr.rs:190，包含逻辑在 194-205 行）：

```rust
candidates.sort_by_key(|tract| (tract.first, std::cmp::Reverse(tract.last)));
let mut outer: Vec<Tract> = Vec::with_capacity(candidates.len());
for tract in candidates {
    let nested = outer.iter().any(|kept| {
        kept.first <= tract.first
            && tract.last <= kept.last
            && same_motif(&kept.pattern, &tract.pattern)   // ← 修复点
    });
    if !nested {
        outer.push(tract);
    }
}
```

`same_motif()`（ssr.rs:103）判定两个（已归一化的）motif 是否互为旋转：

```rust
fn same_motif(a: &str, b: &str) -> bool {
    a.len() == b.len() && b.repeat(2).contains(a)
}
```

等长 + 互为旋转 ⟹ 同一最小周期、同一 cyclic word。tract 的 `pattern` 是
归一化单元在 tract 起点处的旋转，因此该判定对"AT vs TA"、"GAT vs ATG"这类
同一词的相位差成立，对 "AT vs ATATATATC"（长度不同）不成立。

### 修复点 2：跨 motif 的重叠交给贪心消解

`resolve_candidates` 第③步（ssr.rs:207-216）按（完整拷贝数 desc、起点 asc）
贪心接受互不重叠的 tract。上面最小复现的完整走查：

```
候选：AT [0..8)×4、AT [9..17)×4、AT [18..26)×4、
      TA [1..7)×3、TA [10..16)×3、TA [19..25)×3、
      ATATATATC [0..27)×3

① span 去重：无变化
② 包含丢弃（仅同词）：
     TA [1..7)  ⊂ AT [0..8)   同词 ✓ → 丢弃
     TA [10..16) ⊂ AT [9..17)  同词 ✓ → 丢弃
     TA [19..25) ⊂ AT [18..26) 同词 ✓ → 丢弃
     AT [0..8) ⊂ ATATATATC [0..27)? 不同词 → 保留 ★ 修复点
③ 贪心（4 拷贝先于 3 拷贝）：
     接受 AT [0..8)、AT [9..17)、AT [18..26)
     ATATATATC [0..27) 与它们重叠 → 丢弃
输出：三条 AT×4（与参考一致）★ 修复前只输出 9-mer 一条
```

### 为什么"同词嵌套才丢弃"仍然正确

- 同词包含：容器 span ⊇ 被包含 span，两者最小周期相同（归一化保证），
  故容器拷贝数 ≥ 被包含者；被包含者只是同一重复的相位平移子窗口，
  双报告纯冗余。典型：`(ATG)2 ⊂ (GAT)3`。
- 检查对象是"已保留列表"而非全体候选：若 t ⊆ k 且 k 被丢弃（k ⊆ m），
  则 t ⊆ m，m 更早处理且已保留（或又被更大者包含），由 ⊆ 的传递性，
  对保留列表检查与对全体候选检查等价。
- 不同词包含：无周期定理背书，两者是独立检测到的真实 SSR，
  保留候选、由贪心按用户规则（完整拷贝数多者优先、平局取最左）裁决。

## 五、回归验证

### 单元测试（`src/ssr.rs` `mod test`）

| 测试 | 位置 | 覆盖 |
| --- | --- | --- |
| `a_finer_period_tract_survives_a_different_motif_container` | ssr.rs:479 | 本 bug 的最小复现：序列 `ATATATATCATATATATCATATATATCA` + 记录 `AT`/`ATATATATC`，断言三条 AT×4 全部保留 |
| `a_nested_tract_loses_against_its_container` | ssr.rs:468 | 同词嵌套仍丢弃：`(ATG)2 ⊂ (GAT)3`，只报 GAT×4 |

### 数据集三层校验（`scripts/ssr_check.py`）

修复前后 seq 55 的 `--ssr` 输出对比（1-based，区域 [445..553]）：

| 修复前 | 修复后 |
| --- | --- |
| 445..454 AT×5、456..463 AT×4、465..474 AT×5、478..483 AT×3、485..494 AT×5、**514..540 TATATATCA×3**、542..553 AT×6（**缺 4 条 AT×4**） | 445..454 AT×5、456..463 AT×4、465..474 AT×5、478..483 AT×3、485..494 AT×5、**496..503 AT×4、505..512 AT×4、514..519 TA×3、522..529 AT×4、531..538 AT×4**、542..553 AT×6（11 条） |

全数据集三层校验最终 0 违例：level 1（每条输出是精确 motif×copies 区域）、
level 2（每条 tract 落在同周期参考 span 内）、level 3（可达的参考 stretch 全部报告；
不可达者归入贪心消解 / 区域 motif 不相关 / TRF 检测阶段未发现等预期类别）。
seq 55 相关 stretch 现分类为 `已报告` 或 `被更大 tract 消解`。

## 六、相关

- Bug 1（run 跳进：`s = e` 跳过 `[e-period+1, e)` 区间的锚点）由同一轮校验发现，
  修复位于 `extract_record_tracts()`（ssr.rs:111，跳进语句在 ssr.rs:175），
  回归测试 `a_short_anchor_run_does_not_hide_a_longer_one`（ssr.rs:390）。
- 校验与报告工具：`scripts/ssr_check.py`（`--html` 可生成可视化报告），
  用法见 README"SSR 结果的一致性校验"一节。
