# 水稻转录组耐碱性分析程序

## 简介

本程序用于自动分析水稻在碱性胁迫下的转录组数据，通过标准化的流程从原始测序数据到差异表达分析和功能富集分析，帮助研究人员快速识别与耐碱性相关的基因。

## 功能特点

- **数据质量控制**：使用FastQC评估原始数据质量
- **参考基因组比对**：使用HISAT2进行高效比对
- **基因表达定量**：使用featureCounts进行准确计数
- **差异表达分析**：识别碱性胁迫下的差异表达基因
- **功能注释和富集分析**：分析差异基因的生物学功能
- **结果可视化**：生成MA图、火山图和热图等可视化结果

## 安装依赖

程序需要以下软件和库：

### 命令行工具
- FastQC: 用于质量控制
- HISAT2: 用于参考基因组比对
- Samtools: 用于BAM文件处理
- featureCounts: 用于基因表达定量

### Python库
- Python 3.6+
- pandas: 用于数据处理
- numpy: 用于数值计算
- matplotlib: 用于绘图
- seaborn: 用于高级可视化
- scipy: 用于统计分析
- scikit-learn: 用于聚类分析
- pyyaml: 用于读取配置文件

## 安装方法

### 安装命令行工具

**Linux/macOS:**
```bash
# 使用conda安装
conda install -c bioconda fastqc hisat2 samtools subread

# 或使用apt安装 (Ubuntu/Debian)
sudo apt-get install fastqc hisat2 samtools
```

**Windows:**
- 下载并安装 [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/install)
- 在WSL中按照Linux方法安装

### 安装Python库

```bash
pip install pandas numpy matplotlib seaborn scipy scikit-learn pyyaml
```

## 使用方法

1. **准备配置文件**
   - 复制并修改 `config.yaml` 文件，设置样本信息和参考基因组路径

2. **运行分析**
   ```bash
   python rice_alkaline_resistance_analysis.py --config config.yaml
   ```

## 配置文件说明

配置文件 `config.yaml` 包含以下部分：

- **output_dir**: 输出结果的目录
- **sample_info**: 样本信息，包括fastq文件路径和分组
- **reference**: 参考基因组信息，包括基因组索引和GTF文件
- **analysis**: 分析参数，如线程数和差异表达阈值

## 输出结果

分析完成后，输出目录将包含以下子目录：

- **qc**: 质量控制结果
- **alignment**: 比对结果（BAM文件）
- **quantification**: 基因表达定量结果
- **differential_expression**: 差异表达分析结果
- **functional_analysis**: 功能富集分析结果
- **visualization**: 可视化结果（MA图、火山图、热图）

## 示例配置

```yaml
# 水稻转录组耐碱性分析配置文件

# 输出目录
output_dir: output

# 样本信息
sample_info:
  control_1:
    fastq1: data/control_1_R1.fastq.gz
    fastq2: data/control_1_R2.fastq.gz
    group: control
  control_2:
    fastq1: data/control_2_R1.fastq.gz
    fastq2: data/control_2_R2.fastq.gz
    group: control
  treated_1:
    fastq1: data/treated_1_R1.fastq.gz
    fastq2: data/treated_1_R2.fastq.gz
    group: treated
  treated_2:
    fastq1: data/treated_2_R1.fastq.gz
    fastq2: data/treated_2_R2.fastq.gz
    group: treated

# 参考基因组信息
reference:
  genome_index: reference/hisat2_index/rice
  gtf_file: reference/rice_gtf.gtf

# 分析参数
analysis:
  threads: 4
  qvalue_threshold: 0.05
  fold_change_threshold: 2
```

## 注意事项

1. 确保所有依赖软件都已正确安装并添加到系统PATH
2. 参考基因组需要提前构建HISAT2索引
3. 样本文件路径必须正确设置
4. 分析大型数据集时，建议增加线程数以提高速度
5. 对于功能富集分析，建议使用专业的工具如GOATOOLS或clusterProfiler

## 故障排除

- **FastQC失败**: 检查fastq文件路径是否正确
- **HISAT2比对失败**: 检查基因组索引是否正确构建
- **featureCounts定量失败**: 检查GTF文件格式是否正确
- **内存不足**: 对于大型数据集，建议增加系统内存或减少并行线程数

## 扩展功能

本程序可以通过以下方式扩展：

1. 添加更多的差异表达分析方法（如DESeq2、edgeR）
2. 集成专业的功能富集分析工具
3. 添加更多的可视化选项
4. 支持更多的测序数据类型（如single-cell RNA-seq）

## 引用

如果使用本程序发表研究，请引用：

```
Rice Alkaline Resistance Transcriptome Analysis Pipeline. (2026). GitHub: https://github.com/heihuqi333/rice-alkaline-resistance-analysis
```

## 联系方式

如有问题或建议，请联系：
- 邮箱: 1751200760@qq.com
- GitHub: https://github.com/heihuqi333-beep

---

**版本**: 1.0.0
**发布日期**: 2026-02-27
**作者**: Qi Hu
