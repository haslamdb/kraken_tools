# PERMANOVA Analysis

PERMANOVA (Permutational Multivariate Analysis of Variance) is a powerful method for analyzing differences in microbial community composition between groups. Kraken Tools implements PERMANOVA analysis to identify significant associations between metadata variables and overall microbiome structure.

## What is PERMANOVA?

PERMANOVA tests whether the centroids and dispersion of groups defined by categorical variables differ from each other. Unlike methods that examine individual taxa (such as differential abundance tests), PERMANOVA looks at the entire community structure simultaneously.

Key features:
- Non-parametric test using permutations instead of assuming a specific distribution
- Works with any distance metric (Bray-Curtis, UniFrac, Euclidean, etc.)
- Calculates R² values to quantify the percentage of variation explained by each variable
- Handles multiple categorical variables and their interactions

## Command-Line Usage

```bash
kraken-tools permanova \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/permanova/ \
    --group-col "Group" \
    --categorical-vars "Treatment,TimePoint,DiseaseStatus" \
    --distance-metric "bray" \
    --permutations 999 \
    --make-pcoa
```

### Required Parameters

- `--abundance-file`: Path to the processed abundance file (from the processing step)
- `--sample-key`: Path to sample metadata CSV file
- `--output-dir`: Directory for output files

### Optional Parameters

- `--group-col`: Main grouping variable to analyze (default: "Group")
- `--categorical-vars`: Comma-separated list of categorical metadata variables to analyze
- `--distance-metric`: Distance metric to use (options: "bray", "jaccard", "euclidean", "unifrac", default: "bray")
- `--permutations`: Number of permutations for significance testing (default: 999)
- `--make-pcoa`: Generate PCoA plots for significant variables (default: True)
- `--min-group-size`: Minimum number of samples required in each group (default: 3)
- `--transform`: Transformation to apply to abundance data ("clr", "hellinger", "none", default: "clr")

## Output Files

The PERMANOVA analysis produces the following output files:

### Main Results
- `permanova_results_by_r2.csv`: Results sorted by R² (variance explained)
- `permanova_results_by_pvalue.csv`: Results sorted by p-value
- `permanova_variance_explained.pdf`: Bar plot showing R² values for all tested variables

### PCoA Visualizations
- `pcoa_by_[variable].png`: PCoA plots colored by significant variables
- `pcoa_coordinates.tsv`: Raw PCoA coordinates for custom visualization

## Interpreting Results

The main PERMANOVA results file contains the following columns:

- **Feature**: Name of the metadata variable tested
- **p_value**: Statistical significance (p-value) 
- **R2**: Proportion of variance explained (0-1)
- **Significant**: Whether the variable is significant (Yes/No)

Key points for interpretation:
- **p-value < 0.05**: Indicates a statistically significant association
- **R²**: Represents the strength of the association; higher values indicate stronger effects
- PCoA plots help visualize how samples cluster by the significant variables

## Python API Usage

For more flexibility, you can use the Python API:

```python
from kraken_tools.analysis.permanova import run_permanova_analysis

results = run_permanova_analysis(
    abundance_file="processed_abundance.tsv",
    metadata_file="metadata.csv",
    output_dir="results/permanova/",
    categorical_vars=["Treatment", "TimePoint", "DiseaseStatus"],
    distance_metric="bray",
    transform="clr",
    permutations=999,
    make_pcoa=True
)

# Access results programmatically
for feature, result in results.items():
    print(f"{feature}: p-value = {result['p_value']:.4f}, R² = {result['R2']:.4f}")
```

## Example

Here's an example of PERMANOVA analysis on a human gut microbiome dataset with treatment groups and timepoints:

```bash
kraken-tools permanova \
    --abundance-file gut_microbiome_abundance.tsv \
    --sample-key patient_metadata.csv \
    --output-dir results/gut_study/ \
    --categorical-vars "Treatment,Timepoint,AntibioticHistory,DiseaseStatus,Gender" \
    --distance-metric "bray" \
    --transform "clr"
```

This might reveal that Treatment (R²=0.15, p=0.001) and AntibioticHistory (R²=0.08, p=0.003) are significantly associated with microbiome composition, while other variables are not significant.

## Tips for Effective PERMANOVA Analysis

1. **Sample size**: Ensure you have sufficient samples per group (at least 5-10 recommended)
2. **Group balance**: Groups with highly unbalanced sizes can affect results
3. **Multiple testing**: When testing many variables, consider multiple testing correction
4. **Distance metric**: Different distance metrics emphasize different aspects of community structure:
   - Bray-Curtis: Sensitive to both presence/absence and abundance
   - Jaccard: Only considers presence/absence
   - Euclidean: Works well with CLR-transformed data
   - UniFrac: Incorporates phylogenetic relationships (requires phylogenetic tree)
5. **Transform data**: Transform abundance data to account for compositional nature
   - CLR (centered log-ratio): Addresses the compositional nature of microbiome data
   - Hellinger: Less sensitive to rare taxa
