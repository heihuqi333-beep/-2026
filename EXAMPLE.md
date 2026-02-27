# 使用示例

本文件提供了一个完整的水稻转录组耐碱性分析示例，帮助您快速上手使用本程序。

## 准备工作

### 1. 安装依赖

首先，确保您已安装所有必要的依赖：

```bash
# 安装命令行工具
conda install -c bioconda fastqc hisat2 samtools subread

# 安装Python库
pip install pandas numpy matplotlib seaborn scipy scikit-learn pyyaml
```

### 2. 准备数据

- **原始数据**：将您的测序数据（fastq.gz文件）放入 `data/` 目录
- **参考基因组**：
  - 下载水稻参考基因组序列（如IRGSP-1.0）
  - 构建HISAT2索引：
    ```bash
    hisat2-build rice_genome.fa reference/hisat2_index/rice
    ```
  - 下载对应的GTF注释文件

### 3. 配置文件设置

修改 `config.yaml` 文件，设置您的样本信息和参考基因组路径：

```yaml
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

## 运行分析

在命令行中执行以下命令：

```bash
python rice_alkaline_resistance_analysis.py --config config.yaml
```

## 预期输出

分析完成后，您将在 `output/` 目录中看到以下结果：

### 1. 质量控制结果
- `output/qc/` 目录包含每个样本的FastQC报告

### 2. 比对结果
- `output/alignment/` 目录包含每个样本的BAM文件

### 3. 基因表达定量
- `output/quantification/clean_counts.txt` - 基因计数矩阵

### 4. 差异表达分析
- `output/differential_expression/de_results.csv` - 完整的差异分析结果
- `output/differential_expression/significant_genes.csv` - 显著差异表达基因

### 5. 功能富集分析
- `output/functional_analysis/go_enrichment.csv` - GO富集分析结果

### 6. 可视化结果
- `output/visualization/ma_plot.png` - MA图
- `output/visualization/volcano_plot.png` - 火山图
- `output/visualization/heatmap.png` - 差异基因热图

## 结果解读

### 差异表达基因
- 显著差异表达基因定义为：FDR < 0.05 且 |log2(fold change)| > 1
- 正的log2(fold change)表示在处理组中上调
- 负的log2(fold change)表示在处理组中下调

### 功能富集分析
- 富集分析结果按FDR值排序
- 低FDR值表示该GO term在差异基因中显著富集

### 可视化结果
- **MA图**：展示基因表达水平与差异倍数的关系
- **火山图**：展示差异倍数与显著性的关系
- **热图**：展示差异基因在不同样本中的表达模式

## 自定义分析

### 添加更多样本
只需在 `config.yaml` 文件的 `sample_info` 部分添加更多样本信息。

### 调整分析参数
可以修改 `analysis` 部分的参数：
- `threads`：增加线程数以提高分析速度
- `qvalue_threshold`：调整差异表达的显著性阈值
- `fold_change_threshold`：调整差异表达的倍数阈值

### 扩展功能
- **添加新的分析步骤**：在 `RiceAlkalineResistanceAnalysis` 类中添加新的方法
- **修改可视化选项**：调整 `visualization` 方法中的绘图参数
- **集成专业工具**：替换或扩展现有的分析方法

## 故障排除

### 常见错误及解决方法

1. **FastQC失败**
   - 错误信息：`FastQC failed: sample_name`
   - 解决方法：检查fastq文件路径是否正确，确保文件存在且格式正确

2. **HISAT2比对失败**
   - 错误信息：`Alignment failed: sample_name`
   - 解决方法：检查基因组索引是否正确构建，确保索引文件存在

3. **featureCounts定量失败**
   - 错误信息：`Quantification failed`
   - 解决方法：检查GTF文件格式是否正确，确保文件存在

4. **内存不足**
   - 错误信息：`MemoryError`
   - 解决方法：减少并行线程数，或增加系统内存

5. **Python依赖缺失**
   - 错误信息：`ModuleNotFoundError: No module named 'xxx'`
   - 解决方法：使用pip安装缺失的库

## 性能优化

- **大型数据集**：增加 `threads` 参数以提高并行处理速度
- **内存限制**：对于内存有限的系统，减少同时处理的样本数
- **磁盘空间**：确保有足够的磁盘空间存储临时文件和输出结果

## 示例命令

### 1. 构建基因组索引
```bash
hisat2-build Oryza_sativa.IRGSP-1.0.dna.toplevel.fa reference/hisat2_index/rice
```

### 2. 运行分析
```bash
python rice_alkaline_resistance_analysis.py --config config.yaml
```

### 3. 查看结果
```bash
ls -la output/
```

---

通过本示例，您应该能够成功运行水稻转录组耐碱性分析程序，并获得有意义的结果。如果您遇到任何问题，请参考README.md中的故障排除部分或联系作者。