# t-SNE Analysis for Microbiome Visualization

t-SNE (t-Distributed Stochastic Neighbor Embedding) is a powerful dimensionality reduction technique that helps visualize complex, high-dimensional data like microbiome profiles. This module provides robust t-SNE visualization capabilities for taxonomic abundance data.

## Overview

The t-SNE implementation in Kraken Tools allows you to:

1. **Visualize Community Structure**: Reduce high-dimensional taxonomic profiles to 2D visualizations that capture sample relationships
2. **Map Taxonomic Abundance**: Color samples by the abundance of specific taxa of interest
3. **Overlay Metadata**: Visualize how experimental factors correlate with community structure
4. **Generate Publication-Quality Figures**: Create high-resolution plots suitable for publications

This approach is especially valuable for:
- Exploratory data analysis of complex microbiome datasets
- Identifying patterns and clusters in your samples
- Visualizing associations between community structure and metadata variables
- Generating figures for presentations and publications

## Command-Line Usage

```bash
kraken-tools tsne \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/tsne/ \
    --target-taxa "Bacteroides.fragilis,Faecalibacterium.prausnitzii" \
    --categorical-vars "Treatment,TimePoint,DiseaseStatus" \
    --transform clr
```

### Required Parameters

- `--abundance-file`: Path to the processed abundance file (from the processing step)
- `--sample-key`: Path to sample metadata CSV file
- `--output-dir`: Directory for output files

### Optional Parameters

- `--target-taxa`: Comma-separated list of taxa to visualize (default: top 10 by abundance)
- `--categorical-vars`: Comma-separated list of categorical metadata variables to use for coloring
- `--group-col`: Primary grouping variable (included if not in categorical_vars)
- `--transform`: Transformation to apply to abundance data ("clr", "hellinger", "log", "none", default: "clr")
- `--perplexity`: Perplexity parameter for t-SNE, controls the balance between local and global structure (default: 30)
- `--n-iter`: Number of iterations for t-SNE, higher values give more stable results (default: 1000)

## Output Files

The t-SNE analysis produces the following output files:

### Data Files
- `tsne_coordinates.tsv`: Raw t-SNE coordinates for all samples

### Individual Plots
- `tsne_[taxon].png`: Individual plots for each taxon, colored by abundance 
- `tsne_metadata_[variable].png`: Individual plots for each metadata variable

### Grid Plots
- `tsne_all_taxa_grid.png/pdf`: Grid visualization of all taxa in a single figure
- `tsne_all_metadata_grid.png/pdf`: Grid visualization of all metadata variables

### Combined Plots
- `tsne_combined_overview.png/pdf`: Combined visualization showing the most important taxa and metadata

### Multi-page PDFs
- `tsne_abundance_plots.pdf`: Multi-page PDF with all taxon plots
- `tsne_metadata_plots.pdf`: Multi-page PDF with all metadata variable plots

## Interpreting t-SNE Results

t-SNE visualizations require careful interpretation:

- **Sample Proximity**: Samples that are close together in the t-SNE plot have similar microbial profiles
- **Clusters**: Groups of samples that form distinct clusters may represent different community types
- **Local Structure**: t-SNE preserves local neighborhood relationships, not global distances
- **Color Patterns**: When colored by taxon abundance, look for gradients that suggest environmental factors
- **Metadata Correspondence**: When colored by metadata variables, look for separation between categories

Key considerations for interpretation:
- The absolute distances between clusters are not meaningful, only the presence of clustering
- t-SNE results can vary between runs due to its stochastic nature (use `--n-iter` to improve stability)
- The `perplexity` parameter affects clustering; try different values (10-50) if results are unclear

## Data Transformations

Different transformations can emphasize different aspects of your data:

- **CLR (default)**: Centered log-ratio transformation, addresses compositional nature of microbiome data
- **Hellinger**: Reduces impact of highly abundant taxa, highlighting rarer taxa
- **Log**: Simple log transformation that reduces dynamic range
- **None**: No transformation, uses raw abundance values

## Python API Usage

For more flexibility, you can use the Python API:

```python
from kraken_tools.analysis.tsne import run_tsne_analysis

output_dir = run_tsne_analysis(
    abundance_file="processed_abundance.tsv",
    metadata_file="metadata.csv",
    output_dir="results/tsne/",
    target_taxa=["Bacteroides.fragilis", "Faecalibacterium.prausnitzii"],
    categorical_vars=["Treatment", "TimePoint"],
    transform="clr",
    perplexity=30,
    n_iter=1000
)
```

## Example

Here's an example of t-SNE analysis on a gut microbiome dataset with treatment and time variables:

```bash
kraken-tools tsne \
    --abundance-file gut_microbiome.tsv \
    --sample-key patient_data.csv \
    --output-dir results/visualization/ \
    --target-taxa "Bacteroides.fragilis,Faecalibacterium.prausnitzii,Akkermansia.muciniphila,Escherichia.coli" \
    --categorical-vars "Treatment,TimePoint,AntibioticHistory" \
    --transform clr \
    --perplexity 40 \
    --n-iter 2000
```

This might reveal distinct clustering by treatment group and correlation between specific taxa and treatment outcomes.

## Tips for Effective t-SNE Visualization

1. **Data Preparation**: Always transform your data (CLR recommended) to account for compositional nature
2. **Perplexity Tuning**: Try different perplexity values (10-50) to find the best visualization
3. **Iteration Count**: Use at least 1000 iterations for stable results
4. **Sample Size**: Ensure sufficient samples for meaningful visualization (at least 20-30 recommended)
5. **Targeted Analysis**: Focus on specific taxa of interest rather than attempting to visualize all taxa
6. **Metadata Integration**: Include relevant metadata to help interpret clusters

## References

1. van der Maaten, L., & Hinton, G. (2008). Visualizing data using t-SNE. Journal of Machine Learning Research, 9, 2579-2605.
2. Friedman, J., & Alm, E. J. (2012). Inferring correlation networks from genomic survey data. PLoS Computational Biology, 8(9), e1002687.
