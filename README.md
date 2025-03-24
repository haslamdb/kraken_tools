# Kraken Tools

A comprehensive Python package for analyzing taxonomic profiles from Kraken2 and Bracken, with tools for data processing, statistical analysis, and visualization.

## Overview

Kraken Tools provides an end-to-end solution for microbiome analysis:

1. **Preprocessing**: Run quality control and host depletion with KneadData
2. **Taxonomic Classification**: Process sequences through Kraken2 and estimate abundances with Bracken
3. **Processing**: Merge, normalize, and filter taxonomic abundance data
4. **Analysis**: Perform statistical tests, differential abundance analysis, and visualization

## Installation

### Conda Environment Setup (Recommended)

```bash
# Clone the repository
git clone https://github.com/haslamdb/kraken_tools.git
cd kraken_tools

# Create a Conda environment
conda create -n kraken_tools python=3.12 -y
conda activate kraken_tools

# Install dependencies
conda install -c conda-forge -c bioconda pandas numpy scipy scikit-bio scikit-learn scikit-posthocs statsmodels matplotlib seaborn matplotlib-venn tqdm psutil
conda install -c bioconda kraken2 bracken kneaddata

# Install the package
pip install -e .
```

### Standard Installation

```bash
# Install directly using pip (once published to PyPI)
pip install kraken-tools

# Note: External tools (Kraken2, Bracken, KneadData) must be installed separately
```

## Command-Line Usage

Kraken Tools provides a comprehensive set of commands for different stages of microbiome analysis:

### Basic Workflows

| Command | Description |
|---------|-------------|
| `full-pipeline` | Run the complete pipeline from raw reads to analysis |
| `preprocess` | Run only preprocessing (KneadData) for quality control and host depletion |
| `classify` | Run only taxonomic classification (Kraken2 + Bracken) |
| `process` | Process existing Kraken/Bracken files without re-running classification |
| `analyze` | Preliminary group-based analysis on processed abundance data |

### Advanced Analyses

| Command | Description |
|---------|-------------|
| `diff-abundance` | Differential abundance testing with multiple methods (ALDEx2, ANCOM, ANCOM-BC) |
| `glmm` | Generalized Linear Mixed Models for complex experimental designs |
| `permanova` | Permutational Multivariate Analysis of Variance for community-level differences |
| `feature-selection` | Random Forest feature importance analysis for microbiome drivers |
| `rf-shap` | Random Forest with SHAP (SHapley Additive exPlanations) values for interpretable ML |
| `tsne` | Run t-SNE dimensionality reduction and visualization for community structure analysis |

### Utilities

| Command | Description |
|---------|-------------|
| `list-files` | List discovered input files in Kraken/Bracken directories |

### Examples

#### 1. Full Pipeline (Raw Reads to Analysis)

```bash
kraken-tools full-pipeline \
    --input-fastq reads_1.fastq.gz reads_2.fastq.gz \
    --paired \
    --kneaddata-dbs /path/to/kneaddata_db \
    --kraken-db /path/to/kraken_db \
    --bracken-db /path/to/kraken_db/database150mers.kmer_distrib \
    --sample-key metadata.csv \
    --output-dir results/ \
    --group-col "Group" \
    --min-abundance 0.01 \
    --min-prevalence 0.1 \
    --threads 8
```

#### 2. Preprocessing Only (KneadData)

```bash
kraken-tools preprocess \
    --input-fastq reads_1.fastq.gz reads_2.fastq.gz \
    --paired \
    --kneaddata-dbs /path/to/kneaddata_db \
    --output-dir results/preprocessed/ \
    --threads 8
```

#### 3. Taxonomic Classification Only (Kraken2 + Bracken)

```bash
kraken-tools classify \
    --input-fastq clean_reads.fastq \
    --kraken-db /path/to/kraken_db \
    --bracken-db /path/to/kraken_db/database150mers.kmer_distrib \
    --output-dir results/taxonomy/ \
    --taxonomic-level S \
    --threads 8
```

#### 4. Process Existing Kraken/Bracken Files

```bash
kraken-tools process \
    --kreport-dir kraken_reports/ \
    --bracken-dir bracken_files/ \
    --sample-key metadata.csv \
    --output-dir results/processed/ \
    --min-abundance 0.01 \
    --min-prevalence 0.1
```

#### 5. Preliminary Group-Based Analysis
####    Heatmap, PCA, Diversity Metrics, and Statistical Tests


```bash
kraken-tools analyze \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/analysis/ \
    --group-col "Group"
```

#### 6. Differential Abundance Testing

```bash
kraken-tools diff-abundance \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/diff_abundance/ \
    --group-col "Group" \
    --methods aldex2,ancom,ancom-bc
```

#### 7. GLMM Analysis

```bash
kraken-tools glmm \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/glmm/ \
    --formula "Count ~ Group + (1|Subject)" \
    --model negbin
```

#### 8. PERMANOVA Analysis

```bash
kraken-tools permanova \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/permanova/ \
    --categorical-vars "Treatment,TimePoint,DiseaseStatus" \
    --distance-metric bray \
    --transform clr \
    --permutations 999
```

#### 9. Microbiome-Wide Feature Selection with Random Forest

```bash
kraken-tools feature-selection \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/feature_selection/ \
    --predictors "Treatment,TimePoint,Subject,Age,BMI" \
    --distance-metric bray \
    --transform clr
```

#### 10. Organism-Level Feature Selection with Random Forest with SHAP 

```bash
kraken-tools rf-shap \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/rf_shap/ \
    --target-taxa "Bacteroides.fragilis,Faecalibacterium.prausnitzii" \
    --predictors "Treatment,TimePoint,Age" \
    --random-effects "Subject" \
    --transform clr
```

#### 11. t-SNE Visualization

```bash
kraken-tools tsne \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/tsne/ \
    --target-taxa "Bacteroides.fragilis,Faecalibacterium.prausnitzii" \
    --categorical-vars "Treatment,TimePoint" \
    --transform clr
```

#### 12. Utility: List Files

```bash
kraken-tools list-files \
    --kreport-dir kraken_reports/ \
    --bracken-dir bracken_files/
```

### Common Options

| Option | Description |
|--------|-------------|
| `--log-file FILE` | Path to log file (default: stdout) |
| `--log-level LEVEL` | Logging level (DEBUG, INFO, WARNING, ERROR) |
| `--max-memory MB` | Maximum memory usage in MB |
| `--no-interactive` | Non-interactive mode for sample key selection |
| `--use-parallel` | Enable parallel processing |
| `--threads N` | Number of threads to use |
| `--threads-per-sample N` | Threads per sample in parallel mode |
| `--max-parallel N` | Maximum samples to process in parallel |


## Sample Key Format

The sample key CSV file should contain:

- A column with sample identifiers matching the file names in the input directories
- Additional columns for grouping and metadata

**Example**:
```csv
SampleName,Group,Treatment,TimePoint,Subject
sample1,Control,Placebo,Day0,Subject1
sample2,Treatment,Drug,Day0,Subject2
sample3,Control,Placebo,Day7,Subject1
sample4,Treatment,Drug,Day7,Subject2
```

## Output Structure

The output directory will contain the following structure:

```
output_dir/
├── PreprocessedData/           # If preprocessing was run
│   ├── kneaddata_output/       # Clean reads after host removal
│   └── ...
├── TaxonomyData/               # If classification was run
│   ├── kraken_reports/         # Kraken2 reports (.kreport)
│   └── bracken_output/         # Bracken abundance files
├── ProcessedData/              # Processed abundance data
│   ├── KrakenProcessed/        # Processed Kraken reports
│   ├── BrackenProcessed/       # Processed Bracken files
│   └── ...
├── DownstreamAnalysis/         # Analysis results
│   ├── taxonomy_heatmap.svg    # Taxonomic heatmap
│   ├── taxonomy_pca.svg        # PCA plot
│   ├── diversity_metrics.tsv   # Alpha diversity metrics
│   ├── diversity_boxplots.svg  # Diversity visualizations
│   ├── StatisticalTests/       # Statistical test results
│   └── ...
├── DifferentialAbundance/      # If diff-abundance was run
│   ├── aldex2_results.csv      # ALDEx2 results
│   ├── ancom_results.csv       # ANCOM results
│   ├── ancom_bc_results.csv    # ANCOM-BC results
│   ├── method_comparison.txt   # Comparison between methods
│   ├── venn_diagram.png        # Overlap of significant features
│   └── ...
└── GLMM/                       # If GLMM was run
    ├── glmm_results.csv        # Combined GLMM results
    ├── glmm_coefficients_*.png # Coefficient plots
    └── glmm_*.txt              # Individual taxon model results
```

## Python API Usage

You can also use the Python API for more flexibility:

```python
from kraken_tools import run_full_pipeline

# Run the complete pipeline
abundance_file, success = run_full_pipeline(
    sample_key="metadata.csv",
    kreport_dir="kraken_reports/",
    bracken_dir="bracken_files/",
    output_dir="results/",
    group_col="Group",
    min_abundance=0.01,
    min_prevalence=0.1,
    log_file="kraken_analysis.log"
)

# Process files only
from kraken_tools import process_kraken_files_only

abundance_file = process_kraken_files_only(
    sample_key="metadata.csv",
    kreport_dir="kraken_reports/",
    bracken_dir="bracken_files/",
    output_dir="results/",
    taxonomic_level="S",
    log_file="processing.log"
)

# Run differential abundance analysis
from kraken_tools import run_taxonomic_differential_abundance

results = run_taxonomic_differential_abundance(
    abundance_file="abundance.tsv",
    sample_key="metadata.csv",
    output_dir="diff_abundance/",
    group_col="Group",
    methods=["aldex2", "ancom-bc"],
    min_abundance=0.01,
    min_prevalence=0.1
)

# Access results
if 'aldex2' in results:
    significant = results['aldex2'][results['aldex2']['q_value'] < 0.05]
    print(f"Found {len(significant)} significant features with ALDEx2")
```

## Statistical Analysis Methods

### Alpha Diversity Metrics

The following alpha diversity metrics are calculated for each sample:

- **Richness**: Number of observed taxa
- **Shannon Index**: Measures both richness and evenness
- **Simpson Index**: Measures dominance (1 - sum of squared proportions)
- **Evenness**: Shannon diversity divided by log(richness)

### Statistical Tests

- **Kruskal-Wallis Test**: Non-parametric alternative to ANOVA for comparing across multiple groups
- **Dunn's Post-hoc Test**: For pairwise comparisons following significant Kruskal-Wallis results

### Differential Abundance Methods

- **ALDEx2**: Uses a Dirichlet-multinomial model to account for compositional data
  - Only works with two-group comparisons
  - Features with q-value < 0.05 are considered significant

- **ANCOM**: Analysis of Composition of Microbiomes, based on log-ratio testing
  - Works with multiple groups
  - Features with W-ratio > 0.7 are considered significant

- **ANCOM-BC**: ANCOM with bias correction
  - Works with multiple groups
  - Features with q-value < 0.05 are considered significant

### GLMM Analysis

The GLMM (Generalized Linear Mixed Model) module supports:

- **Zero-inflated Poisson models**: For count data with excess zeros
- **Zero-inflated Negative Binomial models**: For overdispersed count data
- **Mixed effects**: Inclusion of random effects for nested/repeated measures design
- **Multiple testing correction**: FDR correction for multiple taxa

## Visualization Types

- **Heatmaps**: Show taxonomic abundances across samples with clustering
- **PCA Plots**: Visualize relationships between samples based on taxonomic profiles
- **Bar Plots**: Display taxonomic composition by group
- **Diversity Plots**: Visualize alpha diversity metrics across samples
- **Volcano Plots**: Show effect size vs. significance for differentially abundant taxa
- **Venn Diagrams**: Show overlap of significant features identified by different methods
- **Coefficient Plots**: Visualize effect sizes from GLMM analysis

## Filtering Parameters

- **min_abundance**: Minimum relative abundance threshold (default: 0.01 = 1%)
  - Filters out taxa with mean abundance below this threshold
  - Helps remove rare taxa that may introduce noise

- **min_prevalence**: Minimum prevalence threshold (default: 0.1 = 10%)
  - Filters out taxa present in fewer than this proportion of samples
  - Helps focus on consistently detected taxa

## Examples

### Example 1: Basic Analysis Workflow

```bash
# Process Bracken files and run standard downstream analysis
kraken-tools process \
    --bracken-dir bracken_files/ \
    --sample-key metadata.csv \
    --output-dir results/ \
    --group-col "DiseaseStatus"
```

### Example 2: End-to-End Pipeline

```bash
# Process raw sequence files through KneadData, Kraken2, Bracken
kraken-tools full-pipeline \
    --input-fastq sample1_R1.fastq.gz sample1_R2.fastq.gz sample2_R1.fastq.gz sample2_R2.fastq.gz \
    --paired \
    --kneaddata-dbs /path/to/human_db \
    --kraken-db /path/to/kraken_db \
    --bracken-db /path/to/kraken_db/database150mers.kmer_distrib \
    --sample-key metadata.csv \
    --output-dir results/ \
    --group-col "DiseaseStatus" \
    --threads 8
```

### Example 3: Focused Differential Abundance Analysis

```bash
# Run only specific differential abundance methods at genus level
kraken-tools diff-abundance \
    --abundance-file abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/ \
    --group-col "DiseaseStatus" \
    --methods aldex2,ancom-bc \
    --min-abundance 0.001 \
    --min-prevalence 0.2
```

### Example 4: GLMM with Repeated Measures

```bash
# Run GLMM analysis with Subject as random effect
kraken-tools glmm \
    --abundance-file abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/glmm/ \
    --formula "Count ~ Group + TimePoint + (1|Subject)" \
    --model negbin
```

### Example 5: Parallel Processing for Large Dataset

```bash
# Process large dataset with parallel execution
kraken-tools full-pipeline \
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

## Citation

If you use Kraken Tools in your research, please cite:

- The original Kraken2 paper: Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. Genome Biol. 2019;20:257.
- Bracken: Lu J, Breitwieser FP, Thielen P, Salzberg SL. Bracken: estimating species abundance in metagenomics data. PeerJ Comput Sci. 2017;3:e104.
- KneadData: McIver LJ, et al. bioBakery: a meta'omic analysis environment. Bioinformatics. 2018;34(7):1235-1237.
- This tool: Haslam D. (2025). Kraken Tools: A comprehensive framework for taxonomic analysis.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
