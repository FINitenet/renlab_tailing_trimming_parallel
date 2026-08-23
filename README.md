<p align="right"><strong>中文</strong> | <a href="README_EN.md">English</a></p>

# miR3End

**A parallel workflow for profiling miRNA 3′-end trimming and tailing**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)](VERSION)

## 概述

miRNA 3′ 末端的非模板加尾（non-templated tailing）与核苷酸剪切（trimming）是
影响 miRNA 稳定性、降解及功能调控的重要转录后修饰。本流程面向植物单端 small
RNA-seq 数据，旨在系统识别并定量 miRNA 3′ 末端的加尾和剪切事件，同时描述不同
修饰长度、末端碱基组成及其在样本间的分布特征。

流程首先对原始 FASTQ reads 进行接头去除、质量控制和长度筛选，并通过与
tRNA/snoRNA 参考序列的精确比对去除潜在污染。过滤后的 reads 随后比对至 miRNA hairpin
参考序列。对于无法直接比对的 reads，流程从其 3′ 端依次剪除 1–10 nt，并在每次
剪切后重新比对，以恢复可能含有非模板末端修饰的 miRNA reads。根据 read 的原始
长度、剪切长度、比对起点以及注释 miRNA 的成熟序列边界，构建每条 miRNA 的
trimming–tailing profile，并进一步统计 5GMC、末端碱基类型和长度分布。

除 miRNA 末端修饰分析外，本流程还提供常规及 UMI small RNA 文库的基因组比对、
mapping statistics 汇总、全长 reads 与 5GMC 长度分布，以及 tail-base 统计表和
可视化结果。流程还可利用 ShortStack 对 18–28 nt reads 进行零错配基因组定位和
多重比对分配，再通过 featureCounts 与 GFF3 biotype 注释计算文库 RNA type 组成。
样本级并行调度与单样本线程控制相互独立，可用于多样本 small RNA-seq 项目的批量分析。

本项目由 `renlab_tailing_trimming_20240111` 重构而来。在保持核心统计定义和主要
分析参数不变的基础上，将原先分散于多个目录的 Python 和 R 分析逻辑整合到当前
仓库，并以单遍 SAM 扫描替代重复扫描策略。Bowtie 索引和物种参考文件由用户按
`resources/README.md` 准备，仓库不再依赖原服务器目录结构。

### 分析流程

```mermaid
flowchart TD
    A[单端 small RNA-seq FASTQ] --> B[样本识别与并行调度]

    B --> C1[接头与低质量序列过滤<br/>Q ≥ 20，长度 12–30 nt]
    C1 --> C2[tRNA/snoRNA 精确比对<br/>去除潜在污染 reads]
    C2 --> C3[miRNA hairpin 精确比对]
    C3 -->|直接比对成功| C5[合并并标准化比对结果]
    C3 -->|未比对 reads| C4[从 3′ 端依次剪除 1–10 nt<br/>每轮重新比对 miRNA hairpin]
    C4 --> C5
    C5 --> D[单遍扫描 SAM]
    D --> D1[11 × 11 trimming–tailing profile]
    D --> D2[5GMC reads 与 sequence-logo 输入]
    D --> D3[全长 reads 与 5GMC 长度分布]
    D1 --> D4[每个 miRNA 一页的<br/>多样本气泡矩阵图]

    B --> E1{文库类型}
    E1 -->|常规文库| E2[接头过滤]
    E1 -->|UMI 文库| E3[UMI 提取、合并与计数]
    E1 -->|UMI 文库| E2
    E3 --> E7[UMI reads 基因组比对<br/>及长度分布]
    E2 --> E4[筛选 18–28 nt reads]
    E4 --> E5[基因组精确比对]
    E5 --> E6[Mapping statistics]
    E4 --> E8[ShortStack 零错配基因组定位<br/>多重比对 reads 分配]
    E8 --> E9[featureCounts biotype 注释]
    E9 --> E10[RNA type 优先级分类<br/>组成与长度分布图]

    D1 --> F[Tail-base 统计与可视化]
    D2 --> F
    E6 --> F
```

### 主要输出

- 每条 miRNA 的 11 × 11 trimming–tailing profile
- 多样本并排的 trimming–tailing 气泡矩阵多页 PDF
- miRNA 总 reads、1–10 nt 加尾及剪切事件的汇总统计
- 5GMC reads 及 sequence-logo 输入文件
- 全长 reads 和 5GMC 的长度分布工作簿及样本级双面板 PDF
- 不同末端碱基及修饰长度的计数、比例和 RPM 结果
- 常规或 UMI small RNA 的基因组比对结果及 mapping statistics
- ShortStack 定位后的 RNA type 计数、比例、长度分布及多样本组成图

## 相比原版本的提升

| 项目 | 原版本 | 新版本 |
| --- | --- | --- |
| 项目结构 | 运行时引用 `TRMRNAseqTools/0_workflow_for_srna_seq.py` 等外部脚本 | 已将相关逻辑和 R 脚本整合到本仓库；外部参考资源及其参数见 `resources/README.md` |
| 多样本运行 | 主要按样本串行处理 | 使用 `--jobs` 控制样本级并行，可同时处理多个样本 |
| 资源控制 | 工具线程数分散在不同脚本中 | 使用 `--threads-per-sample` 统一控制单样本线程数，总线程数约为 `jobs × threads-per-sample` |
| profile 统计 | Perl 针对每条 miRNA 重复扫描 SAM，共进行 538 次全文件扫描 | Python 单遍扫描 SAM，同时生成 profile、summary 和 5GMC 结果 |
| 中断续跑 | 需要手动判断已完成步骤 | `--resume` 根据最终输出跳过已完成的样本和步骤 |
| 运行前检查 | 缺少统一预览入口 | `--dry-run` 可预览样本、参数、索引及待执行命令 |
| 外部依赖 | 依赖外部 Python/Perl 脚本，并使用部分非必要工具 | 移除 Perl、`seqkit` 和 `dnapilib` 依赖；保留流程必需的软件及 Bowtie 索引 |
| 日志与报错 | 多层脚本调用，定位失败步骤较困难 | 统一入口、分步骤日志和样本级错误信息，便于排查及重跑 |

### 性能优化

耗时最明显的优化位于 `profile` 步骤。旧实现会针对 538 个 miRNA 分别重读一次
SAM；新实现只遍历一次 SAM，并在同一过程中完成全部统计。测试中，一个约 1.2 GB
的 SAM 文件可在约 37 秒内完成该步骤。新旧版本生成的 profile、summary 和 5GMC
三个关键结果已进行逐字节比较，结果一致。

FASTQ 和 SAM 均采用流式读取，避免把整个大文件一次性载入内存。R 汇总脚本也对
稀疏碱基统计进行了调整，降低无效数据展开造成的额外开销。

### 并行方式

并行以“样本”为单位：`--jobs` 表示最多同时处理的样本数，
`--threads-per-sample` 表示每个样本调用 Bowtie 等工具时可使用的线程数。例如：

```bash
python3 main.py -i 1_rawdata -o results --jobs 2 --threads-per-sample 8
```

该配置最多同时运行 2 个样本，总 CPU 需求约为 16 线程。这样既能提升多样本吞吐量，
也能避免为每个样本无限制启动进程。

### 独立性和可复现性

原流程引用的外部小 RNA workflow、统计逻辑和 R 脚本已经纳入当前仓库，不再需要
保持原服务器上的脚本目录结构。所有可执行程序默认从当前 Conda 环境或 `PATH`
查找；所有参考索引和注释均通过命令行参数或 `MIR3END_*` 环境变量提供。仓库中不再
包含 `/bios-store1`、`/home` 或 `/usr/local` 的机器专属默认路径。外部参考资源的
要求记录在 `resources/README.md`，分析时应同时保存其来源、版本和校验值。

`test/` 提供不依赖外部基因组或比对程序的合成最小数据和可直接运行的 profile
smoke test。正式分析前仍建议先用少量代表性样本验证本机的软件版本、索引路径和
资源配置。

## 安装与运行环境

推荐使用版本化 Conda 环境；`environment.yml` 同时包含 Python、R 和命令行依赖：

```bash
conda env create -f environment.yml
conda activate mir3end-0.1.0
python main.py --version
python test/run_smoke_test.py
```

仅安装 Python 依赖时可使用 `requirements.txt`，但完整流程仍需要 `bowtie`、
`trim_galore`、`ShortStack`、`featureCounts` 和 `Rscript`。这些工具均从 `PATH`
查找，也可用相应参数覆盖。

## 发布与引用

当前发布候选版本为 `0.1.0`。软件采用 [MIT License](LICENSE)，引用元数据见
[`CITATION.cff`](CITATION.cff)。仓库同时提供 `.zenodo.json`；合并发布分支后，在
Zenodo 中连接该 GitHub 仓库并创建 `v0.1.0` Release，即可获得版本 DOI。DOI 生成后
应回填 `CITATION.cff`、README 和论文摘要中的软件 URL。

## 快速开始

```bash
python3 main.py \
  -i /path/to/1_rawdata \
  -o /path/to/project \
  --mir-hairpin /path/to/hairpin_index_prefix \
  --trsno /path/to/trsno_index_prefix \
  --genome-index /path/to/genome_index_prefix \
  --genome-fasta /path/to/genome.fa \
  --rnatype-annotation /path/to/annotation.gff3 \
  --jobs 2 \
  --threads-per-sample 8
```

总线程需求约为 `jobs × threads-per-sample`。建议第一次使用 `--jobs 2`，确认
内存和 CPU 负载后再增加。

输入文件默认匹配 `*.fastq.gz`。样本名由文件名去除 `--suffix` 后得到，并将连字符
`-` 统一转换为下划线 `_`。如果文件名为 `sample_R1.fastq.gz`，可使用：

```bash
python3 main.py -i 1_rawdata -o project_output --suffix _R1.fastq.gz
```

## 输出目录

新项目默认将可重建的中间文件放入 `work/`，将用于检查、统计和发表的结果放入
`results/`。两条分析分支分别存放，避免编号目录相互混杂：

```text
project_output/
├── work/
│   ├── tailing_trimming/
│   │   ├── trimmed_fasta/
│   │   └── remapping/<sample>/
│   └── srna/
│       ├── trimmed/
│       ├── umi/
│       └── genome_mapping/
└── results/
    ├── tailing_trimming/
    │   ├── profiles/
    │   ├── 5gmc/
    │   ├── sequence_logo/
    │   ├── length/
    │   │   ├── <sample>_len_dist.xlsx
    │   │   └── plots/<sample>_length_distribution.pdf
    │   ├── bubble/tailing_trimming_bubble.pdf
    │   └── tailbase/
    │       ├── tables/<sample>/
    │       ├── plots/<sample>/
    │       └── workspace/
    └── srna/
        ├── mapping_summary.tsv
        └── rnatype/
            ├── rnatype_summary.tsv
            ├── rnatype_summary.xlsx
            ├── tables/<sample>_rnatype_by_length.tsv
            └── plots/
                ├── rnatype_composition.pdf
                ├── rnatype_composition.png
                └── <sample>_rnatype_by_length.pdf
```

`--layout auto` 是默认设置：新项目采用上述结构；如果检测到旧版的
`3_remapping/` 或 `4_163.results/`，则继续使用原目录，避免重复计算。也可以显式
使用 `--layout organized` 或 `--layout legacy`。

## 分步运行、功能与结果解读

流程可拆分为 8 个步骤。下表中的“直接图形”仅表示该步骤本身是否产图；不直接产图的
步骤通常生成后续统计和绘图所需的标准化中间结果。

| 步骤 | 主要功能 | 前置结果 | 直接图形 |
| --- | --- | --- | --- |
| `preprocess` | miRNA 分支质控、污染过滤和逐级 3′ remapping | 原始 FASTQ | FastQC HTML |
| `profile` | 构建 11 × 11 修剪–加尾矩阵、summary、5GMC 和 sequence-logo 输入 | `preprocess` | 不直接产图 |
| `bubble` | 按 miRNA 绘制多样本修剪–加尾气泡矩阵 | `profile` | 多页 PDF |
| `srna` | 常规/UMI 文库处理、18–28 nt 筛选和基因组比对 | 原始 FASTQ | FastQC HTML |
| `mapping_summary` | 汇总 Bowtie 日志中的总 reads、mapped reads 和比对率 | `srna` | 不直接产图 |
| `rnatype` | ShortStack 定位、featureCounts 注释和唯一 RNA type 分类 | `srna` | RNA type 组成图和长度分布图 |
| `length` | 统计 miRNA 全长 reads 与 5GMC 的长度分布 | `preprocess` | Excel 内嵌图和样本级 PDF |
| `tailbase` | 汇总 trimming、总体 tailing 和非模板 tailing 的长度及碱基组成 | `profile`、`mapping_summary` | 每个样本 2 个 PDF |

```mermaid
flowchart LR
    P[preprocess] --> R[profile]
    R --> B[bubble]
    P --> L[length]
    R --> T[tailbase]
    S[srna] --> M[mapping_summary]
    S --> Y[rnatype]
    M --> T
```

以下图片均由固定的匿名示意数据生成，仅用于说明图形结构，不代表真实样本或预期的
生物学分布。可运行 `python3 docs/generate_demo_figures.py` 重新生成 README 示例图。

### 1. `preprocess`：miRNA 末端修饰预处理与 remapping

该步骤使用 Trim Galore 去接头、去除低质量及含 `N` 的 reads，并保留 12–30 nt
序列。reads 首先精确比对 tRNA/snoRNA 参考以去除潜在污染，再精确比对 miRNA
hairpin；未比对序列从 3′ 端逐次剪除 1–10 nt 后重新比对。最终合并为带有剪切信息的
SAM，供后续 profile 与长度分析使用。

```bash
python3 main.py -i 1_rawdata -o project_output --steps preprocess -j 2
```

主要结果：

- `work/tailing_trimming/trimmed_fasta/<sample>_trimmed_fastqc.html`：去接头后 FastQC 报告
- `work/tailing_trimming/trimmed_fasta/<sample>_trimmed.fasta`：格式化的 reads
- `work/tailing_trimming/remapping/<sample>/merged-alignment-sorted.sam`：合并后的核心比对结果
- `work/tailing_trimming/remapping/<sample>/remapping.log`：tRNA/snoRNA 和 miRNA 比对日志

FastQC 报告用于检查碱基质量、接头残留、序列长度和高频序列。README 中的简化示例如下；
正式分析应打开每个样本对应的交互式 HTML 报告。

![QC output example](docs/images/qc_example.png)

### 2. `profile`：修剪–加尾矩阵与 5GMC 提取

该步骤单遍扫描每个样本的合并 SAM，并依据 mature miRNA 的起止位置，将 reads 归入
0–10 nt trimming 与 0–10 nt tailing 的 11 × 11 矩阵。同时提取 5′ genomic matching
component（5GMC）reads，并生成 sequence-logo 输入和按 miRNA 汇总的计数表。

```bash
python3 main.py -o project_output --steps profile --resume
```

主要结果：

- `results/tailing_trimming/profiles/<sample>.txt`：样本完整 profile
- `results/tailing_trimming/profiles/<sample>.summary.txt`：总量、trimming 和 1–10 nt tailing 汇总
- `results/tailing_trimming/profiles/<sample>/<miRNA>.txt`：单个 miRNA 的 11 × 11 矩阵
- `results/tailing_trimming/5gmc/<sample>.5GMC`：5GMC reads
- `results/tailing_trimming/sequence_logo/<sample>/`：每个 miRNA 的 sequence-logo 输入

本步骤不直接生成图片。单 miRNA 矩阵由 `bubble` 可视化，5GMC 则由 `length` 和
`tailbase` 用于后续统计。

### 3. `bubble`：多样本 trimming–tailing 气泡图

该步骤将同一 miRNA 在多个样本中的 11 × 11 profile 并排绘制。横轴为 trimming
长度，纵轴为 tailing 长度，圆面积表示该组合在当前 miRNA profile 中的相对 read
比例。PDF 每页对应一个 miRNA，适合逐个比较样本间的末端修饰模式。

```bash
python3 main.py \
  -o project_output \
  --steps bubble \
  --filelist filelist \
  --bubble-cols 2
```

主要结果：`results/tailing_trimming/bubble/tailing_trimming_bubble.pdf`。`--bubble-cols`
控制每行样本数，`--bubble-output` 可覆盖输出路径。

![Trimming-tailing bubble plot example](docs/images/bubble_example.png)

### 4. `srna`：small RNA 文库处理与基因组比对

普通文库完成接头过滤后，保留 18–28 nt reads 并使用 Bowtie 对基因组进行零错配比对。
UMI 文库还会识别 linker、提取 12 nt UMI、合并相同 insert–UMI 组合，并另行输出 UMI
reads 的比对和长度统计。`--umi-flag 1` 表示 UMI 文库，默认值 `2` 表示普通文库。

```bash
python3 main.py -i 1_rawdata -o project_output --steps srna --umi-flag 2 -j 2
```

主要结果：

- `work/srna/trimmed/<sample>_trimmed_fastqc.html`：small RNA 分支 FastQC 报告
- `work/srna/genome_mapping/<sample>_trimmed_len.fq.gz`：18–28 nt reads
- `work/srna/genome_mapping/<sample>_aligned.sam`：基因组比对结果
- `work/srna/genome_mapping/<sample>.mapresults.txt`：Bowtie 比对日志
- UMI 模式下的 `work/srna/umi/<sample>_umi.fa.gz`、比对结果和 `<sample>_umi_dist.csv`

本步骤的直接图形仍为 FastQC HTML，其图形结构可参考步骤 1；UMI 长度输出为 TSV/CSV
统计表，不自动生成 PDF。

### 5. `mapping_summary`：基因组比对统计汇总

该步骤解析所有样本的 Bowtie 日志，统一汇总输入 reads 数、至少一次成功比对的 reads
数、比对率和 Bowtie reported alignment 数。它用于项目级质控，也是 `tailbase` 计算
标准化指标时的输入之一。

```bash
python3 main.py -o project_output --steps mapping_summary
```

主要结果为 `results/srna/mapping_summary.tsv`，包含 `Sample`、`Total_reads`、
`Mapped_reads`、`Mapped_rate(%)` 和 `Reported` 五列。本步骤不直接绘图，以便用户按实验
设计自行制作跨样本 QC 图或纳入统计报告。

### 6. `rnatype`：RNA type 分类与组成可视化

该步骤将 18–28 nt reads 经 ShortStack 零错配定位和多重比对分配后，使用
featureCounts 根据 GFF3 `biotype` 注释，并按预设优先级为每条已定位 read 确定唯一
RNA type。图中颜色与图例类别一一对应，堆叠从上到下遵循图例顺序；纵坐标统一为
`Percentage (%)`。

```bash
python3 main.py \
  -i 1_rawdata \
  -o project_output \
  --steps srna,rnatype \
  --jobs 2 \
  --threads-per-sample 8
```

主要结果：

- `results/srna/rnatype/rnatype_summary.tsv` 和 `.xlsx`：跨样本 RNA type 汇总
- `results/srna/rnatype/tables/<sample>_rnatype_by_length.tsv`：样本按 read 长度分类结果
- `results/srna/rnatype/plots/rnatype_composition.pdf` 和 `.png`：多样本总体组成
- `results/srna/rnatype/plots/<sample>_rnatype_by_length.pdf`：样本内按长度的组成

![RNA type output example](docs/images/rnatype_example.png)

### 7. `length`：全长 reads 与 5GMC 长度分布

该步骤从 miRNA remapping SAM 统计原始 read 长度及对应的 5GMC 长度，并以 mapped
miRNA reads 总数为分母计算比例。每个样本同时产生 Excel 工作簿和双面板 PDF，用于
判断完整 reads 与去除末端修饰后的基因组匹配部分是否具有不同长度峰。

```bash
python3 main.py -o project_output --steps length --resume
```

主要结果：

- `results/tailing_trimming/length/<sample>_len_dist.xlsx`
- `results/tailing_trimming/length/plots/<sample>_length_distribution.pdf`

工作簿包含 `total reads mapped to miRNA`、`distribution of whole reads` 和
`distribution of 5GMC` 三个工作表，首页嵌入两幅图。如果旧工作簿已存在但缺图，使用
`--steps length --resume` 可读取现有统计并补图，无需重新扫描 SAM。

![Whole-read and 5GMC length output example](docs/images/length_example.png)

### 8. `tailbase`：末端碱基和修饰长度汇总

该步骤整合 profile summary、5GMC、mapping summary 和 miRNA 注释，对 trimming、总体
tailing 以及非模板 tailing 分别统计长度、A/T/C/G 碱基组成、百分比和 RPM。总体图包含
所有推断的加尾事件；非模板图进一步排除与参考序列相同的模板性延伸，更适合描述真正的
non-templated tailing。

```bash
python3 main.py -o project_output --steps tailbase --resume
```

主要结果：

- `results/tailing_trimming/tailbase/tables/<sample>/<sample>_tail_summary.xlsx`
- `results/tailing_trimming/tailbase/tables/<sample>/<sample>_tail_nt_summary.xlsx`
- `results/tailing_trimming/tailbase/tables/<sample>/<sample>_trim_summary.xlsx`
- `results/tailing_trimming/tailbase/plots/<sample>/<sample>_tailing_trimming_length.pdf`
- `results/tailing_trimming/tailbase/plots/<sample>/<sample>_tailing_trimming_nt_length.pdf`

![Tail-base output example](docs/images/tailbase_example.png)

从已有 `work/tailing_trimming/remapping/<sample>/` 继续运行时，可以省略原始 FASTQ，
或用 `--filelist` 指定样本。旧目录通过 `--layout legacy` 继续支持。`--resume` 根据每一步
的最终输出跳过已完成任务；`--dry-run` 用于检查样本、参数、索引和即将执行的命令。

### RNA type 分类规则与可配置参数

`rnatype` 步骤复用 `srna` 步骤产生的 18–28 nt FASTQ。ShortStack 使用零错配、
`mmap=u`、`bowtie_m=1000` 和 `ranmax=50`，与旧版 `0_workflow.sh` 保持一致；随后
featureCounts 使用 GFF3 中 `gene,ncRNA_gene` 的 `biotype` 属性进行注释。重叠注释按
miRNA primary transcript、phasi/tasiRNA、hc-siRNA、transposable element、lncRNA、
protein-coding、snoRNA、snRNA、rRNA、tRNA 的顺序确定唯一 RNA type。未比对到基因组的
reads 不进入分母；已比对但没有匹配注释的 reads 归入 `unassigned`。

```bash
python3 main.py \
  -i 1_rawdata \
  -o project_output \
  --steps srna,rnatype \
  --jobs 2 \
  --threads-per-sample 8
```

如果 `srna` 的长度过滤结果已经存在，可只运行 `--steps rnatype --resume`。默认参考为
TAIR10，可通过 `--genome-fasta` 和 `--rnatype-annotation` 替换。featureCounts 默认沿用
旧流程的非链特异设置 `--rnatype-stranded 0`；若文库方向已通过实验设计或独立检查确认，
可设置为 `1`（正向）或 `2`（反向）。自定义 RNA type 重叠优先级可通过
`--rnatype-priority` 指定，每行一个类别、由上到下优先级递减。

## 目录内容

```text
main.py                 统一入口和步骤调度
pipeline/               Python 实现；样本级受控并行
pipeline/bubble.py      多样本 trimming–tailing 气泡矩阵图
pipeline/rnatype.py     ShortStack/featureCounts RNA type 分类与绘图
assets/tail_base_summary.R  本地 R 汇总脚本
docs/                   README 匿名示意图及其可重复生成脚本
resources/              R/Python 所需的本地参考表（Git 忽略）
requirements.txt        Python 依赖
```

默认的拟南芥 Bowtie 索引可以通过 `--mir-hairpin`、`--trsno` 和
`--genome-index` 覆盖。运行 `python3 main.py --help` 查看全部参数。

## 补充材料：定量标准化与 trimming–tailing 对角线判读

### S1. 不同统计量及其分母

本流程中的 raw count、profile 内百分比和 RPM 回答的问题不同，不能相互替代。设
`C(s,m,t,a)` 为样本 `s` 中 miRNA（或相同成熟序列家族）`m` 在 trimming 长度 `t`、
tailing 长度 `a` 格子的 read count，则各方法定义如下。

| 方法 | 计算方式 | 主要用途 | 局限与推荐用法 |
| --- | --- | --- | --- |
| raw count | `C(s,m,t,a)` | 检查支持某个修饰状态的实际 read 数 | 受文库深度影响，不直接用于跨样本定量；应始终随百分比或 RPM 一起检查 |
| GitHub 原始 bubble 面积 | `C(s,m,t,a) / sum(t,a) C(s,m,t,a) × 1500` | 比较同一 miRNA 在样本内的 trimming–tailing 构成 | 分母是该 miRNA 的整个 11 × 11 profile；这是 profile 内比例而不是 RPM。低丰度 miRNA 的少量 reads 也可能形成大气泡，不能据此判断绝对丰度回复 |
| 全部 reads 构成百分比 | `100 × C / 该成熟序列家族全部分配 reads` | 展示 canonical 与各种修饰状态的完整组成 | 推荐作为 bubble 颜色；适合回答“该修饰占这个 miRNA 的多少” |
| modified-only 百分比 | `100 × C / 该家族全部 modified reads`，并排除 `(0,0)` | 放大观察已发生修饰的 reads 如何分布 | 只描述修饰谱，不能与 canonical 比较；当 modified raw count 很低时不稳定，必须同时报告 raw count |
| all-miRNA RPM | `10^6 × C / 样本全部已分配 miRNA reads` | 跨样本比较 miRNA 或修饰状态丰度 | 当前推荐的主定量分母；不依赖 ShortStack RNA-type 注释，但仍需检查 HEN1 缺失是否改变总体 miRNA 池 |
| PTGS-proxy RPM | `10^6 × C / phasi_tasiRNA reads` | 检查结果是否依赖 all-miRNA 分母 | 作为敏感性分析；`phasi_tasiRNA` 只是可获得的 PTGS-siRNA proxy，不代表全部 PTGS-siRNA |
| rRNA-derived RPM | `10^6 × C / rRNA reads` | 第二套外部参考分母敏感性分析 | 容易受 rRNA 降解、污染和建库差异影响，不建议单独作为主结论 |
| combined stable-pool RPM | `10^6 × C / (all-miRNA + phasi/tasiRNA + rRNA-derived reads)` | 给 bubble 面积提供跨样本可比的丰度尺度 | 推荐用于候选 bubble 的面积；颜色仍使用家族内百分比，以同时表达“组成”和“丰度” |

在 `nrpd1` 导致 hc-siRNA 整体塌陷的实验中，不应使用总 small-RNA reads 或 hc-siRNA
总量作为主要分母，否则分母本身的强烈生物学变化会人为放大其他 RNA 类别。建议采用：

1. 以 **all-miRNA RPM** 作为跨样本丰度的主结果；
2. 同时报告 PTGS-proxy RPM、rRNA-derived RPM 和 combined stable-pool RPM 的分母敏感性；
3. bubble 颜色使用家族内百分比，面积使用 combined stable-pool RPM，并在图旁保留 raw count；
4. 若不同分母下效应方向及候选排序一致，结论才视为对标准化选择稳健。

当前描述性比较不使用 edgeR/TMM，也不计算 P-value 或 FDR。先对 raw count 按上述参考池
计算 RPM，再计算各生物学组的 mean RPM；组间倍数变化定义为
`log2(mean RPM_A / mean RPM_B)`。不添加 pseudocount 时，零均值产生的
`Inf`、`-Inf` 或 `NaN` 应原样保留，并同时报告每组 raw count 和重复间离散程度。
`|log2FC| >= 1` 仅表示描述性的至少两倍变化，不等于统计显著。

### S2. 对角线点的含义

在 bubble 图中，横轴为 trimming 长度，纵轴为 tailing 长度。因此对角线
`trimming = tailing = k` 表示成熟 miRNA 的 3′ 端缺失了 `k` 个参考碱基，同时 read
末端又多出 `k` 个碱基。在 5′ 起点相同的前提下，这类 read 的总长度与注释成熟 miRNA
相同，可称为“长度补偿型 trim-and-tail”或“3′ 端替换”。例如成熟序列末端为 `...ACGA`，
read 末端为 `...ACUU` 时，可被记为 `(trim=2, tail=2)`。

对角线不能自动解释为真实生物学修饰。当前 GitHub 流程先要求 read 与 hairpin 零错配
比对；未比对 read 再从 3′ 端逐个剪除 1–10 nt 并重新零错配比对。因此，一条与成熟
miRNA 等长、但末端连续 `k` 个碱基不同的 read，会在算法上自然落到 `(k,k)`。对角线
信号可能是以下来源的混合：

- 真实的先 trimming、后 non-templated tailing；在 `hen1` 背景下，以 U 为主并受
  `heso2` 基因型影响时具有较强生物学合理性；
- 3′ 端测序错误或末端低质量，尤其容易形成 `(1,1)`；
- mature miRNA 注释末端偏差、alternative DCL processing 或同源家族成员的末端差异；
- adapter 残留，或 hairpin/genome 上的模板性延伸被误判为非模板 tail；
- Bowtie `-a` 多重比对使同一 read 被计入多个相关 hairpin，从而重复强化相同模式；
- 低 raw count 格子因 profile 内百分比显示方式而在视觉上被放大。

### S3. 对角线判定与报告标准

建议对每个候选的对角线格子提取原始 read 序列和实际 tail sequence，并按下表判读。

| 检查项目 | 更支持真实 trim-and-tail | 更支持技术或算法伪影 |
| --- | --- | --- |
| tail 组成 | 以 U 为主，或存在明确的酶偏好 | 碱基混杂、富集 adapter motif |
| 基因型方向 | `hen1` 富集，并在 `heso2` 背景中按预期改变 | WT 和所有突变体近似一致 |
| 组织与表型 | 花序中出现与育性回复一致的变化，叶片差异较弱 | 组织和表型方向不一致 |
| 生物学重复 | 至少多个重复方向一致且有足够 raw count | 单个重复驱动，或仅有少量 reads |
| 3′ 碱基质量 | 末端质量正常 | 末端质量明显下降，特别是 `(1,1)` |
| precursor/genome 检查 | tail 不存在于对应 hairpin 或基因组延伸序列 | 所谓 tail 实际可由模板解释 |
| 替代比对策略 | 允许少量错配或使用显式 3′ 端解析后仍存在 | 改变精确匹配策略后对角线明显塌陷 |
| 家族与多重比对 | 按相同成熟序列合并家族、每条 read 总权重为 1 后仍存在 | 只在重复注释成员中同时出现 |

最低限度应同时报告：对角线 raw count、占全部家族 reads 的百分比、占 modified reads
的百分比、RPM、tail 碱基组成以及三个生物学重复。可增加两个汇总指标：

```text
Diagonal fraction = sum(k >= 1) C(trim=k, tail=k) / 全部 modified reads
U-diagonal fraction = 对角线上 U-tail reads / 全部对角线 reads
```

对角线点只有在 raw count 充分、重复稳定、tail 以 U 为主、排除模板性延伸，并且基因型
及组织变化与表型方向一致时，才适合作为候选分子回复证据。单独的 `(1,1)` 点或仅由
profile 内大气泡支持的点，应优先视为待验证信号，而不是直接下生物学结论。
