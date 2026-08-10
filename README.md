# Parallel miRNA tailing/trimming workflow

这是从 `renlab_tailing_trimming_20240111` 重新整理出的独立版本。Python、Perl、R
脚本和运行期数据文件都保存在本目录，不再调用 `TRMRNAseqTools` 或其他脚本目录。

## 运行环境

- Python 3.8+
- Python 包见 `requirements.txt`
- 命令行工具：`trim_galore`、`bowtie`、`Rscript`
- R 包：`tidyverse`、`openxlsx`、`data.table`、`reshape2`、`lubridate`
- Bowtie 索引仍属于大型参考数据，通过参数指定，不复制进代码目录

## 快速开始

```bash
python3 main.py \
  -i /path/to/1_rawdata \
  -o /path/to/project \
  --jobs 2 \
  --threads-per-sample 8
```

总线程需求约为 `jobs × threads-per-sample`。建议第一次使用 `--jobs 2`，确认
内存和 CPU 负载后再增加。

输入文件默认匹配 `*.fastq.gz`。样本名由文件名去除 `--suffix` 后得到，并将连字符
`-` 统一转换为下划线 `_`。如果文件名为 `sample_R1.fastq.gz`，可使用：

```bash
python3 main.py -i 1_rawdata -o results --suffix _R1.fastq.gz
```

## 分步运行

步骤顺序为：

1. `preprocess`：接头过滤以及 miRNA tailing/trimming remapping
2. `profile`：单遍扫描 SAM，生成 163 profile 和 5GMC 汇总；替代旧版对每条
   miRNA 重读一次 SAM 的 Perl 实现
3. `srna`：常规或 UMI 小 RNA genome mapping
4. `mapping_summary`：汇总 Bowtie 日志
5. `length`：生成每个样本的长度分布工作簿
6. `tailbase`：生成 tail base R 统计和图表

例如只运行前两个步骤：

```bash
python3 main.py -i 1_rawdata -o results --steps preprocess,profile -j 2
```

从已有 `3_remapping/<sample>/` 继续运行时，可以省略原始 FASTQ，或用
`--filelist` 指定样本。`--resume` 会根据每一步的最终输出跳过已完成样本；
`--dry-run` 用于检查样本、参数、索引和即将执行的命令。

UMI 文库使用 `--umi-flag 1`，普通小 RNA 文库使用默认值 `2`。

## 目录内容

```text
main.py                 统一入口和步骤调度
pipeline/               Python 实现；样本级受控并行
assets/tail_base_summary.R  本地 R 汇总脚本
resources/              R/Python 所需的小型参考表
requirements.txt        Python 依赖
```

默认的拟南芥 Bowtie 索引可以通过 `--mir-hairpin`、`--trsno` 和
`--genome-index` 覆盖。运行 `python3 main.py --help` 查看全部参数。
