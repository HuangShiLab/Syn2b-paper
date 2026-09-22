# HPC 任务清单

更新：2026-09-17。整理所有必须在工作站/HPC 上执行、本地无法完成的任务，
按优先级排序。每项含：前置条件、命令、产出、QC 门槛、计算量、以及它关闭
稿件中的哪一条 open item。

---

## 执行状态（2026-09-21 更新）

| 任务 | Job ID | 状态 |
|---|---|---|
| 基因组下载 | — | ✅ 完成（304,905 在盘；12,564 个 ENA+NCBI 双缺，覆盖 96.0%） |
| T2/T3/T4 | 4071690 等 | ✅ 完成，结果已回灌两篇稿件 |
| t1a reps skani 三角矩阵 | 4076217 | ✅ 完成（46 分钟） |
| ANI 列通道（逐簇 skani） | 4086086 已停；4096924 | 🔄 4096924 幂等续跑（20,909/22,535 簇） |
| 预消化 v1 | 4092774 | ✅ 完成但发现 ENA 头部缺陷（见下） |
| ENA 无效 TGT 清除 | 4093181 | ✅ 完成（~22.6 万个 ID 均为字面量 "ENA" 的 TGT 已删） |
| 重消化（净化头部后） | 4096453（48 shards） | ✅ 完成（2026-09-22 07:54–07:58） |
| 普查主阵列（mem=180G, mega=12） | 4096923 取消；4096932 | 🔄 修复缺失 FASTA 处理后重启；0/24,063 完成 |

**2026-09-21 事故与修复**：ENA browser API 的 FASTA 头以 `>ENA|...` 开头，
导致所有 ENA 源基因组的 TGT ID 均为字面量 "ENA" — syn2b synteny 在含 ≥2 个
此类基因组的目录上立即报 duplicate genome id 退出。修复：digest 前剥离
`ENA|` 前缀（digest_all.py / census_worker.py），受污染 TGT 全部清除重消化；
同时 shard 进程增加单任务异常隔离（坏任务不再杀死整个 shard），普查阵列
内存提升至 180G、mega 并发门 12。

---

## 执行状态（2026-09-17，已通过 SSH 提交）

| 任务 | Job ID | 状态 |
|---|---|---|
| T2 43,334 对 syn2bani v0.1.1 重算（8×24h） | 4071690 | RUNNING（~46% @ 16 min） |
| T3+T4 队列 skani 0.3.2 + syn2bani triangle | 4071698 | RUNNING |
| T1a reps skani 65,703 三角矩阵（预筛层） | 4071749 | PENDING（排队） |
| 对照组（renamed-copy=0；双工具口径记录） | — | 完成 controls.md |

**T1 普查前置发现**：`genomes_all` 只含 65,703 个代表基因组；种内普查还需
251,839 个成员基因组，其中 239,370 个可从 NCBI 下载（12,469 个已被撤回，
QC 将如实记录；普查覆盖 96.1%）。下载清单已生成：
`data/gtdb-r207/download_missing/urls_missing_genomes.tsv`（URL+accession）。
**I/O 节点下载脚本**：`download_missing/io_download_missing.sh`（8 并发、
断点续传、跳过已有、失败重试，进度每 5 分钟记录；登录节点 15 分钟限制与
I/O 节点不可 SSH，故需经 OnDemand 在 I/O 节点启动，或明确授权分块登录节点
下载）。二进制：syn2b 64717ab / syn2bani v0.1.1 已构建，
skani 0.3.2 在 `tools/bin/`。

---

## 环境速查（来自现有脚本的实际路径）

| 项 | 值 |
|---|---|
| 存储 | `/lustre1/g/aos_shihuang/` |
| 四仓 | `$L/{Syn2b, Syn2bANI, Syn2b-paper, Syn2bANI-paper}`（`L=/lustre1/g/aos_shihuang`） |
| GTDB R207 基因组 | `/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all`（`{accession}.fna` 命名） |
| SynTracker 队列组装 | `/lustre1/g/aos_shihuang/data/syntracker_validation/assemblies/` |
| 分区 | `amd`（历史脚本）；SLURM array 上限 ~95 个任务 |
| 约束 | **少长作业优于多短作业**；运行中须审计磁盘（配额） |

**先构建二进制（每个 clone 一次）：**
```bash
cd $L/Syn2b      && git pull && cargo build --release   # -> target/release/syn2b
cd $L/Syn2bANI   && git pull && cargo build --release   # -> target/release/syn2bani (v0.1.1)
# skani >= 0.2.1（解决 skani 版本 open item）：
cargo install skani --version 0.3.2   # 或模块加载系统里的 >= 0.2.1
skani --version                       # 确认
```

---

## T1 · GTDB R207 种内结构普查（旗舰，优先级最高）

**目标**：全部 711,020,841 个种内配对的 junction/inverted fraction/ANI。
**关闭**：REVIEW_5 V1（应用规模）；新 Results §7 + Zenodo 资源。
**计算量**：结构 ~3,560 core-hours + ANI ~40–80 + digest ~4；存储 ~110 GB TGT + ~50 GB 输出（审计器管住）。

```bash
# 0) 同步并准备（一次性，登录节点即可）
cd $L/Syn2b-paper && git pull
WORK=$L/syn2b_census; mkdir -p $WORK/logs
cp results/census/tasks.jsonl $WORK/                 # 24,063 个任务（本地已生成）
python3 scripts/gtdb_census/prepare_workdir.py \
    --taxonomy data/gtdb_metadata/accession_taxonomy_r207.tsv.gz \
    --genome-dir $L/data/gtdb-r207/genomes_all \
    --tasks $WORK/tasks.jsonl --workdir $WORK
#    -> manifest.json / cluster_accessions/ / missing_genomes.txt
#    检查 missing_genomes.txt 数量，记录到 $WORK/PREP_NOTES.md

# 1) 对照组（必做，几分钟；写进 $WORK/controls.md）
SYN=$L/Syn2b/target/release/syn2b
$SYN digest -i <EDL933.fna> -o /tmp/e.tgt -e BcgI,AlfI,AloI,FalI
$SYN digest -i <EDL933 改名副本.fna> -o /tmp/e2.tgt -e BcgI,AlfI,AloI,FalI
mkdir /tmp/ctl && cd /tmp/ctl && cp /tmp/e.tgt /tmp/e2.tgt . && $SYN synteny --input . --output m
#    期望：对照自身 = 0 junctions
#    EDL933 vs Sakai = 2 junctions（两套实现交叉验证的锚点）

# 2) 提交普查（6 个长任务；断点续跑，可重复 sbatch）
sed -i 's|WORKDIR=.*|WORKDIR='"$WORK"'|; s|SYN2B=.*|SYN2B='"$SYN"'|' \
    scripts/gtdb_census/census_long_job.slurm
#    按实际配额改 SOFT_GB / HARD_GB / MIN_FREE_GB（默认 250/320/50）
sbatch scripts/gtdb_census/census_long_job.slurm

# 3) ANI 列（可与 2 并行，幂等）
python3 scripts/gtdb_census/skani_ani_pass.py \
    --workdir $WORK --skani ~/.cargo/bin/skani --workers 8

# 4) dnadiff 抽查（T5 的 200 对；可与 2 并行）

# 5) 合并 + QC
python3 scripts/gtdb_census/merge_census.py --workdir $WORK --out $L/Syn2b-paper/results/census
#    QC 门槛：unique pairs 与 ΣC(n,2) 一致（缺基因组除外，逐簇列出）；
#    重复行 = 0；对照 = 0/2。把 census_qc.md 回传仓库。
```

**运行中监控**：`tail -f $WORK/worker_*.log | grep AUDIT`（PAUSE/RESUME 事件）；
`du -sh $WORK`；配额触发时审计器自动回收已完成簇的 TGT。

---

## T2 · 43,334 对 Syn2bANI 结构重算（撤回项，优先级 2）

**关闭**：Syn2b-paper §5 撤回注记（r=0.414 偏相关）；Syn2bANI-paper 同项。
**计算量**：参照 s13 历史 ~760 core-hours；用 v0.1.1 二进制。

现有 `run_s2b_slice.sh` 是幂等的（跳过已有输出行），直接复用，但按
「少长作业」改成 8×24h：

```bash
cd $L/Syn2bANI-paper && git pull
WORK=$L/Syn2bANI-paper/results/gtdb50k
export GTDB50K_WORK=$WORK GTDB50K_GENOMES=$L/data/gtdb-r207/genomes_all
export SYN2BANI=$L/Syn2bANI/target/release/syn2bani   # v0.1.1
# 旧输出清空后重算（旧列含 bug 统计量，勿混用）：
mv $WORK/s2b_out $WORK/s2b_out_prefix_v0.1.1.bak
# 以 NSLICES=8 起 8 个长任务（脚本按 NSLICES 切片）：
sbatch --array=0-7%8 --cpus-per-task=4 --mem=8G --time=1-00:00:00 \
    --wrap 'bash scripts/gtdb50k/run_s2b_slice.sh $SLURM_ARRAY_TASK_ID'
# 完成后合并（复用既有 merge 脚本或逐行 cat s2b_out/*.tsv）
# 然后在 Syn2b-paper 重算偏相关（撤回项的正式数字）：
cd $L/Syn2b-paper && python3 scripts/sv_reanalysis.py results/gtdb50k
```

**QC 门槛**：自比较/近重复对的 breakpoint_count = 0；高 ANIm 子集的
junction 数与 Syn2b 通道一致性抽查（同 pair 两工具 ±2 内）。

---

## T3 · skani ≥0.2.1 重跑（优先级 3）

**关闭**：skani 版本 open item（两篇稿）。

```bash
# a) SynTracker 四队列（Syn2b-paper Figure 4 的 ANI 轴）
SK=~/.cargo/bin/skani   # >= 0.2.1
ASM=$L/data/syntracker_validation/assemblies
for c in Escherichia_coli_hypermutator Helicobacter_pylori \
         Neisseria_gonorrhoeae Streptomyces_rimosus; do
    $SK triangle $ASM/${c}/*.fna > skani_${c}.tsv
done
#    -> 覆盖 data/syntracker_validation/skani/*.tsv（保留旧版为 .v010.bak），
#    本地重新生成 Figure 4 并检查数值漂移（>=94.6% 区间预期 <0.05%）

# b) Syn2bANI-paper 的仿真阶梯与小规模比较：按其 README 的 skani 命令
#    用新版重跑并替换结果文件（结果文件名不变，旧版备份）。
```

---

## T4 · SynTracker 四队列的 syn2bani 通道（优先级 3）

**关闭**：Syn2bANI-paper open item（目前只重算了 Syn2b 通道）。

```bash
S2B=$L/Syn2bANI/target/release/syn2bani   # v0.1.1
ASM=$L/data/syntracker_validation/assemblies
for c in <四个队列目录名>; do
    $S2B triangle --verbose $ASM/${c}/*.fna > syn2bani_${c}.tsv
done
#    -> data/syntracker_validation/syn2bani/（旧文件是 bug 前产物，已停用）。
#    QC：每个队列自比较 breakpoint_count = 0；中位数与 Syn2b 通道一致
#    （0/3/7/10 量级）。
```

---

## T5 · 普查 dnadiff 抽查（200 对，优先级 2，随 T1 并行）

```bash
cd $L/Syn2b-paper
# 从普查输出抽样 200 对（分层：>=97 ANI 100 对、95-97 50 对、<95 50 对），
# 复用 Syn2bANI-paper/scripts/gtdb50k/run_dnadiff_slice.sh 的逐对 dnadiff 模式，
# 对比 junction 与 dnadiff 事件数的秩相关（预期 rho >= 0.6，同 §6 的 AUC 量级）。
```

**产出**：`results/census/dnadiff_spotcheck.tsv` + 一段 §7 的验证语句。

---

## T6（可选，T1 落地后）· 完整基因组大规模验证（V3）

RefSeq complete 种内配对全量 + 10% 抽样 dnadiff 真值 → §5 完整基因组面板
升级为主图。复用 T1 管线，输入换成 complete-level 清单（本地
`fetch_vignette_genomes.py` 的过滤器可扩展为全量清单）。

---

## 本地（等外置盘，非 HPC）

**cagPAI 三分类分层检验**：挂载 `/Volumes/MoneyCat` 后（Seagate Basic 盘
已查无该 metadata），用队列 metadata 对 island_internal / island_boundary /
island_spanned 三类做疾病分期分层检验。Syn2bANI-paper 的
`case_studies/h_pylori_cagpai/results/cagpai_states_extended_filtered.tsv`
已在仓库，只缺 metadata。

---

## 计算量汇总

| 任务 | core-hours | 形态 |
|---|---:|---|
| T1 结构普查 | ~3,560 | 6×16核×96h（少长作业） |
| T1 ANI 列 | 40–80 | 单节点数小时 |
| T2 43,334 重算 | ~760 | 8×4核×24h |
| T3 skani 重跑 | <20 | 单节点 |
| T4 syn2bani 队列 | <10 | 单节点 |
| T5 dnadiff 抽查 | ~5 | 单节点 |
| **合计** | **~4,400–4,500** | |

## 完成判定（对稿件）

| 任务 | 稿件条目 |
|---|---|
| T1 | 新 Results §7「GTDB 结构普查」+ Zenodo 资源 + 摘要重写 |
| T2 | §5 撤回注记替换为正式偏相关数字 |
| T3 | Methods skani 版本注记删除，换实测 |
| T4 | Syn2bANI-paper Fig 6 双通道 |
| T5 | §7 dnadiff 锚定段落 |
| T6 | §5 完整基因组面板升主图 |
