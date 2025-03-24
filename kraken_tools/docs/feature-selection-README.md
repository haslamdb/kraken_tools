# Feature Selection with Random Forest

Random Forest feature selection helps identify the most important variables that drive differences in microbiome composition. This approach complements PERMANOVA by not only detecting significant associations but also ranking their relative importance.

## Overview

The feature selection implementation in Kraken Tools uses Random Forest regression to:

1. Quantify the importance of metadata variables in explaining microbiome distances
2. Rank features by their impact on community structure
3. Visualize the relative importance of different clinical and environmental factors

This analysis is especially valuable for:
- Exploratory analysis to identify key drivers of microbiome variation
- Determining which variables to include in more complex statistical models
- Understanding the relative importance of different factors

## Command-Line Usage

```bash
kraken-tools feature-selection \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/feature_selection/ \
    --predictors "Treatment,TimePoint,Subject,Age,BMI" \
    --n-estimators 100 \
    --distance-metric "bray" \
    --transform "clr"
```

### Required Parameters

- `--abundance-file`: Path to the processed abundance file (from the processing step)
- `--sample-key`: Path to sample metadata CSV file
- `--output-dir`: Directory for output files

### Optional Parameters

- `--predictors`: Comma-separated list of metadata variables to analyze (default: all columns)
- `--n-estimators`: Number of trees in the Random Forest (default: 100)
- `--max-features`: Maximum number of features to consider (default: "sqrt")
- `--distance-metric`: Distance metric to use (options: "bray", "jaccard", "euclidean", default: "bray")
- `--transform`: Transformation to apply to abundance data ("clr", "hellinger", "none", default: "clr")
- `--test-size`: Proportion of data to use for testing (default: 0.2)
- `--random-state`: Random seed for reproducibility (default: 42)

## Output Files

The feature selection analysis produces the following output files:

### Results Files
- `feature_importance.csv`: Table of all features ranked by importance score
- `feature_importance_plot.pdf`: Visual representation of feature importance
- `rf_model_summary.txt`: Performance metrics and model details

### Optional Outputs
- `rf_predictions.csv`: Predicted vs. actual distance values for model assessment
- `partial_dependence_plots.pdf`: Visualization of how specific features affect distances

## Interpreting Results

The main feature importance file contains:

- **Feature**: Name of the metadata variable
- **Importance**: Random Forest feature importance score (higher = more important)
- **Rank**: Ranking of the feature (1 = most important)

Key considerations for interpretation:
- Feature importance is relative, not absolute
- The actual importance score depends on the specific Random Forest implementation
- Features with very low importance scores may not be meaningful predictors

## Python API Usage

For more flexibility, you can use the Python API:

```python
from kraken_tools.analysis.feature_selection import run_feature_selection

results = run_feature_selection(
    abundance_file="processed_abundance.tsv",
    metadata_file="metadata.csv",
    output_dir="results/feature_selection/",
    predictors=["Treatment", "TimePoint", "Subject", "Age", "BMI"],
    n_estimators=100,
    distance_metric="bray",
    transform="clr"
)

# Access results programmatically
for feature, importance in results.items():
    print(f"{feature}: {importance:.4f}")
```

## Example

Here's an example of feature selection on a soil microbiome dataset with environmental variables:

```bash
kraken-tools feature-selection \
    --abundance-file soil_microbiome_abundance.tsv \
    --sample-key environmental_data.csv \
    --output-dir results/soil_study/ \
    --predictors "pH,Temperature,Moisture,NitrogenContent,PhosphorusContent,Depth,Season" \
    --n-estimators 200 \
    --distance-metric "bray"
```

This might reveal that pH, Moisture, and Temperature are the most important predictors of soil microbiome composition.

## Methodology

The feature selection process works as follows:

1. **Preprocessing**:
   - The abundance data is transformed (CLR, Hellinger, or no transformation)
   - Pairwise distances between samples are calculated (Bray-Curtis, Jaccard, or Euclidean)
   - Metadata variables are encoded (numerical) or one-hot encoded (categorical)

2. **Creating Feature Vectors**:
   - For each pair of samples, the absolute difference in metadata values is calculated
   - These pairwise metadata differences become the predictors

3. **Random Forest Regression**:
   - The model tries to predict the microbiome distance between samples based on metadata differences
   - Feature importance is extracted from the trained model

4. **Evaluation**:
   - The model's accuracy is assessed on a test set
   - Performance metrics include R², RMSE, and MAE

## Tips for Effective Feature Selection

1. **Feature Selection**: Start with all available metadata and let the algorithm determine important features
2. **Sample Size**: Ensure you have a sufficiently large dataset (at least 30-50 samples recommended)
3. **Transformations**: Try different data transformations to see if they affect feature importance
4. **Categorical Variables**: Be cautious with categorical variables that have many levels
5. **Correlated Features**: Be aware that correlated features may share importance
6. **Validation**: Consider cross-validation for more robust feature importance estimates
