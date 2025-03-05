# Kraken Tools

A comprehensive Python package for analyzing taxonomic profiles from Kraken2 and Bracken, with tools for data processing, statistical analysis, and visualization.

## Table of Contents
- [Installation](#installation)
    - [Conda Environment Setup](#conda-environment-setup)
    - [Standard Installation](#standard-installation)
- [Overview](#overview)
- [Command-Line Usage](#command-line-usage)
- [Python API Usage](#python-api-usage)
- [Processing Kraken2/Bracken Files](#processing-kraken2bracken-files)
- [Downstream Analysis](#downstream-analysis)
- [Differential Abundance Analysis](#differential-abundance-analysis)
- [Visualization](#visualization)
- [Statistical Testing](#statistical-testing)
- [Advanced Configuration](#advanced-configuration)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)

---

## Installation

### Conda Environment Setup

It is **highly recommended** to use a Conda environment to manage the dependencies for `kraken_tools`.

**Steps:**

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/haslamdb/kraken_tools.git
   cd kraken_tools
   ```

2. **Create a Conda Environment:**
   ```bash
   conda create -n kraken_tools python=3.12 -y
   conda activate kraken_tools
   ```

3. **Install Dependencies:**
   ```bash
   conda install -c conda-forge -c bioconda pandas numpy scipy scikit-bio scikit-learn scikit-posthocs statsmodels matplotlib seaborn matplotlib-venn tqdm psutil
   conda install -c bioconda kraken2 bracken
   ```

4. **Install the Package:**
   ```bash
   pip install -e .
   ```

### Standard Installation

If you prefer not to use Conda, you can install directly using pip:

```bash
# Clone the repository
git clone https://github.com/yourusername/kraken_tools.git
cd kraken_tools

# Install in development mode
pip install -e .
```

Or install directly using pip (once published to PyPI):
```bash
pip install kraken-tools
```

**Required Dependencies:**
- pandas
- numpy
- scipy
- scikit-bio
- scikit-learn
- scikit-posthocs
- statsmodels
- matplotlib
- seaborn
- matplotlib-venn
- tqdm
- psutil

**External Requirements:**
- Kraken2 and Bracken must be installed separately and available in your PATH

---

## Overview

**Kraken Tools** provides an end-to-end solution for taxonomic analysis:

1. **Preprocessing**: Run quality control and host depletion with KneadData
2. **Taxonomic Classification**: Process raw sequences through Kraken2
3. **Abundance Estimation**: Estimate abundances with Bracken
4. **Processing**: Merge, normalize, and filter abundance data
5. **Downstream Analysis**: Perform statistical tests, visualization, and differential abundance analysis

---

## Command-Line Usage

The package provides a command-line tool `kraken-tools` for running the full analysis pipeline.

### End-to-End Pipeline (Raw Sequences to Analysis)

```bash
kraken-tools --run-preprocessing --input-fastq reads_1.fastq reads_2.fastq --paired \
    --kneaddata-dbs /path/to/kneaddata_db \
    --kraken-db /path/to/kraken_db \
    --bracken-db /path/to/kraken_db/database150mers.kmer_distrib \
    --sample-key /path/to/metadata.csv \
    --output-dir /path/to/output \
    --group-col "Group" \
    --threads 8
```

### Analysis of Existing Kraken2/Bracken Files

```bash
# Full pipeline with defaults
kraken-tools --sample-key /path/to/metadata.csv \
    --kreport-dir /path/to/kreport_dir \
    --bracken-dir /path/to/bracken_dir \
    --output-dir /path/to/output \
    --group-col "Group"
```

### Common Options

```bash
# Skip Kraken or Bracken processing
kraken-tools --sample-key metadata.csv --kreport-dir kreports/ --bracken-dir bracken/ \
    --output-dir results/ --skip-kraken

# Skip downstream analysis
kraken-tools --sample-key metadata.csv --kreport-dir kreports/ --bracken-dir bracken/ \
    --output-dir results/ --skip-downstream

# Run differential abundance analysis
kraken-tools --sample-key metadata.csv --kreport-dir kreports/ --bracken-dir bracken/ \
    --output-dir results/ --run-diff-abundance
```

### Parallel Processing for Large Datasets

```bash
kraken-tools --run-preprocessing --input-fastq reads_*.fastq --paired \
    --kneaddata-dbs /path/to/kneaddata_db \
    --kraken-db /path/to/kraken_db \
    --bracken-db /path/to/kraken_db/database150mers.kmer_distrib \
    --sample-key /path/to/metadata.csv \
    --output-dir /path/to/output \
    --group-col "Group" \
    --use-parallel \
    --threads-per-sample 4 \
    --max-parallel 8
```

---

## Python API Usage

You can also use the Python API for more flexibility and integration with your own scripts.

### Run End-to-End Preprocessing and Analysis

```python
from kraken_tools import run_preprocessing_and_analysis

results = run_preprocessing_and_analysis(
    input_fastq=["reads_1.fastq", "reads_2.fastq"],
    sample_key="metadata.csv",
    output_dir="results",
    paired=True,
    threads=8,
    kneaddata_dbs="/path/to/kneaddata_db",
    kraken_db="/path/to/kraken_db",
    bracken_db="/path/to/bracken_db",
    group_col="Group"
)
```

### Process and Analyze Existing Kraken2/Bracken Files

```python
from kraken_tools import run_full_pipeline

abundance_file, success = run_full_pipeline(
    sample_key="/path/to/sample_key.csv",
    kreport_dir="/path/to/kreport_dir",
    bracken_dir="/path/to/bracken_dir",
    output_dir="/path/to/output",
    group_col="Group",
    run_diff_abundance=True,
    log_file="kraken_analysis.log"
)

if success:
    print("Analysis completed successfully!")
    print(f"Abundance file: {abundance_file}")
```

### Processing Files Only

```python
from kraken_tools import process_kraken_files_only

abundance_file = process_kraken_files_only(
    sample_key="/path/to/sample_key.csv",
    kreport_dir="/path/to/kreport_dir",
    bracken_dir="/path/to/bracken_dir",
    output_dir="/path/to/output",
    log_file="kraken_processing.log"
)
```

### Running Differential Abundance Analysis

```python
from kraken_tools import run_taxonomic_differential_abundance

results = run_taxonomic_differential_abundance(
    abundance_file="/path/to/abundance_file.tsv",
    sample_key="/path/to/sample_key.csv",
    output_dir="/path/to/results",
    group_col="Treatment",
    methods=["aldex2", "ancom-bc"],
    min_abundance=0.01,
    min_prevalence=0.1
)

# Access the results
if 'aldex2' in results:
    significant = results['aldex2'][results['aldex2']['q_value'] < 0.05]
    print(f"Found {len(significant)} significant features with ALDEx2")
```

---

## Processing Kraken2/Bracken Files

Kraken Tools automates several key processing steps for Kraken2/Bracken output files:

- **Merging**: Combines files from multiple samples into a single table
- **Normalization**: Converts abundance values to relative abundance (0-1) or other formats
- **Filtering**: Removes low-abundance taxa based on thresholds

The output directory structure will look like:

```
output_dir/
├── BrackenProcessed/
│   ├── ProcessedFiles_bracken_s_merged.tsv
│   └── ProcessedFiles_bracken_s_relabundance.tsv
├── KrakenProcessed/
│   ├── ProcessedFiles_kraken_s_merged.tsv
│   └── ProcessedFiles_kraken_s_relabundance.tsv
└── ...
```

---

## Downstream Analysis

The downstream analysis module performs several types of analyses on the processed data:

- **Taxonomic Visualizations**: Heatmaps, bar plots, etc.
- **Principal Component Analysis (PCA)**: Visualizes similarities between samples
- **Diversity Metrics**: Alpha diversity metrics (Richness, Shannon, Simpson, etc.)
- **Statistical Tests**: Performs Kruskal-Wallis tests followed by Dunn's post-hoc tests
- **Differential Abundance Analysis**: Identifies significantly different taxa between groups

Results are organized in the `DownstreamAnalysis` directory:

```
output_dir/
├── DownstreamAnalysis/
│   ├── taxonomy_heatmap.svg
│   ├── taxonomy_pca.svg
│   ├── taxonomy_bars_by_group.svg
│   ├── diversity_metrics.tsv
│   ├── diversity_boxplots.svg
│   ├── StatisticalTests/
│   │   ├── kruskal_wallis_results.csv
│   │   ├── significant_comparisons.csv
│   │   └── dunn_posthoc_tests/
│   │       ├── dunn_taxon1.csv
│   │       ├── dunn_taxon2.csv
│   │       └── ...
│   └── ...
└── ...
```

---

## Differential Abundance Analysis

Kraken Tools includes implementations of three popular differential abundance methods:

- **ALDEx2**: Uses a Dirichlet-multinomial model to account for compositional data; works only with two-group comparisons
- **ANCOM**: Analysis of Composition of Microbiomes, based on log-ratio testing; works with multiple groups
- **ANCOM-BC**: ANCOM with bias correction, addressing limitations of the original ANCOM; works with multiple groups

### Running Differential Abundance Analysis

```bash
# Basic usage
kraken-tools --sample-key samples.csv --kreport-dir kreports/ --bracken-dir bracken/ \
    --output-dir results/ --group-col Group --run-diff-abundance

# Specify methods
kraken-tools --sample-key samples.csv --kreport-dir kreports/ --bracken-dir bracken/ \
    --output-dir results/ --group-col Group --run-diff-abundance \
    --diff-methods aldex2,ancom
```

### Output Files

Differential abundance analysis produces the following files in `output_dir/DifferentialAbundance/`:

- `aldex2_results.csv`: Results from ALDEx2 analysis
- `ancom_results.csv`: Results from ANCOM
- `ancom_bc_results.csv`: Results from ANCOM-BC
- `method_comparison.txt`: Comparison of significant features across methods
- `venn_diagram.png`: Venn diagram showing overlap of significant features
- `aldex2_volcano.png`: Volcano plot of ALDEx2 results
- `ancom_top_features.png`: Bar plot of top features by ANCOM W-ratio

### Interpreting Results

- **ALDEx2**: Features with *q-value* < 0.05 are considered significant
- **ANCOM**: Features with *W-ratio* > 0.7 are considered significant
- **ANCOM-BC**: Features with *q-value* < 0.05 are considered significant

---

## Visualization

Kraken Tools generates several visualization types:

- **Heatmaps**: Show taxonomic abundances across samples with clustering
- **PCA Plots**: Show the relationships between samples based on their taxonomic profiles
- **Bar Plots**: Display taxonomic composition by group
- **Diversity Plots**: Visualize alpha diversity metrics across samples
- **Volcano Plots**: Show effect size vs. significance for differentially abundant taxa
- **Venn Diagrams**: Show overlap of significant features identified by different methods

---

## Statistical Testing

### Kruskal-Wallis Test

A non-parametric alternative to ANOVA for comparing taxa abundance across multiple groups. Results include:

- Test statistic
- p-value
- Adjusted p-value (q-value) using Benjamini-Hochberg FDR correction
- Binary indicator of significance (Reject_H0)

### Dunn's Post-hoc Test

For taxa with significant Kruskal-Wallis results, a Dunn's post-hoc test identifies which specific group pairs differ significantly.

---

## Advanced Configuration

### Sample Key Format

The sample key CSV file should contain:

- A column with sample identifiers matching the file names in the input directories
- Additional columns for grouping and other metadata

**Example**:
```csv
SampleName,Group,Treatment,TimePoint
sample1,Control,Placebo,Day0
sample2,Treatment,Drug,Day0
sample3,Control,Placebo,Day7
```

### Filtering Low-abundance Taxa

You can filter low-abundance taxa based on abundance and prevalence thresholds:

```bash
# Filter taxa with < 1% mean abundance or present in < 10% of samples
kraken-tools --sample-key metadata.csv --bracken-dir bracken/ \
    --output-dir results/ --min-abundance 0.01 --min-prevalence 0.1
```

### Taxonomic Level Selection

Specify the taxonomic level for analysis (S=species, G=genus, etc.):

```bash
# Analyze at genus level
kraken-tools --sample-key metadata.csv --bracken-dir bracken/ \
    --output-dir results/ --taxonomic-level G
```

### Customizing Output

```bash
# Use custom output prefix and directory structure
kraken-tools --sample-key metadata.csv --bracken-dir bracken/ \
    --output-dir my_analysis --output-prefix MyPrefix
```

---

## Examples

### Example 1: Basic Analysis Workflow

```bash
# Process Bracken files and run standard downstream analysis
kraken-tools --sample-key metadata.csv \
    --bracken-dir bracken_files/ \
    --output-dir results/ \
    --group-col "DiseaseStatus"
```

### Example 2: End-to-End Pipeline

```bash
# Process raw sequence files through KneadData, Kraken2, Bracken
kraken-tools --run-preprocessing \
    --input-fastq sample1_R1.fastq.gz sample1_R2.fastq.gz sample2_R1.fastq.gz sample2_R2.fastq.gz \
    --paired \
    --kneaddata-dbs /path/to/human_db /path/to/contaminants_db \
    --kraken-db /path/to/kraken_db \
    --bracken-db /path/to/kraken_db/database150mers.kmer_distrib \
    --sample-key metadata.csv \
    --output-dir results/ \
    --group-col "DiseaseStatus" \
    --threads 8
```

### Example 3: Focused Differential Abundance Analysis

```bash
# Run only specific differential abundance methods, examining genus level
kraken-tools --sample-key metadata.csv \
    --bracken-dir bracken_files/ \
    --output-dir results/ \
    --taxonomic-level G \
    --group-col "DiseaseStatus" \
    --run-diff-abundance \
    --diff-methods aldex2,ancom-bc \
    --min-abundance 0.001 \
    --min-prevalence 0.2
```

### Example 4: Parallel Processing for Large Dataset

```bash
# Process large dataset with parallel execution
kraken-tools --run-preprocessing \
    --input-fastq sample*_R*.fastq.gz \
    --paired \
    --kneaddata-dbs /path/to/human_db \
    --kraken-db /path/to/kraken_db \
    --bracken-db /path/to/kraken_db/database150mers.kmer_distrib \
    --sample-key metadata.csv \
    --output-dir results/ \
    --use-parallel \
    --threads-per-sample 4 \
    --max-parallel 6
```

---

## Troubleshooting

### Common Issues

1. **File Format Issues**
   - **Problem**: "No valid Kraken/Bracken files found for any samples"
   - **Solution**: Check file naming patterns and directory structure. Use `--list-files` option to see detected files.

2. **Sample Key Issues**
   - **Problem**: Sample identifiers in metadata don't match file names
   - **Solution**: Make sure your sample key's sample identifier column matches the sample names in your files.

3. **Statistical Test Errors**
   - **Problem**: "No significant taxa found after FDR correction"
   - **Solution**: Try a less stringent significance threshold with `--min-abundance` and `--min-prevalence` options.

4. **Memory Issues**
   - **Problem**: Out of memory errors with large datasets
   - **Solution**: Use `--max-memory` option, or process in batches.

### Getting Help

If you encounter issues not covered in this documentation, please:

- Check the log file for detailed error messages
- Set `--log-level DEBUG` for more verbose output
- Open an issue on the GitHub repository with a description of the problem and relevant log entries

---

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use Kraken Tools in your research, please cite:

- The original Kraken2 paper: Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. Genome Biol. 2019;20:257.
- Bracken: Lu J, Breitwieser FP, Thielen P, Salzberg SL. Bracken: estimating species abundance in metagenomics data. PeerJ Comput Sci. 2017;3:e104.
- This tool: [Your Name]. (2025). Kraken Tools: A comprehensive framework for taxonomic analysis.
