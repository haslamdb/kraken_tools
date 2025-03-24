# Random Forest with SHAP and Linear Mixed Models

This advanced analysis module combines Random Forest modeling with SHAP (SHapley Additive exPlanations) values and Linear Mixed Models (LMM) to provide detailed insights into factors that influence the abundance of specific taxa of interest.

## Overview

While PERMANOVA and general feature selection look at community-level differences, this analysis focuses on individual taxa. For each taxon of interest, the module:

1. Builds a Random Forest model to predict its abundance
2. Uses SHAP values to explain which factors affect that taxon's abundance
3. Applies Linear Mixed Models to account for random effects (e.g., repeated measures)

This approach is especially valuable for:
- Identifying clinical or environmental factors that influence specific taxa
- Understanding non-linear relationships between metadata and taxon abundance
- Accounting for hierarchical data structures (e.g., multiple samples per subject)

## Command-Line Usage

```bash
kraken-tools rf-shap \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/rf_shap/ \
    --target-taxa "Bacteroides.fragilis,Escherichia.coli" \
    --predictors "Treatment,TimePoint,Age" \
    --random-effects "Subject" \
    --transform "clr"
```

### Required Parameters

- `--abundance-file`: Path to the processed abundance file
- `--sample-key`: Path to sample metadata CSV file
- `--output-dir`: Directory for output files

### Optional Parameters

- `--target-taxa`: Comma-separated list of taxa to analyze (default: top 10 most abundant)
- `--predictors`: Comma-separated list of predictor variables (default: all columns)
- `--random-effects`: Variables to treat as random effects in mixed models (default: None)
- `--transform`: Transformation to apply to abundance data (default: "clr")
- `--n-estimators`: Number of trees in the Random Forest (default: 100)
- `--test-size`: Proportion of data to use for testing (default: 0.2)
- `--mixed-model`: Type of mixed model to use ("lmer" or "glmm", default: "lmer")

## Output Files

The RF-SHAP analysis produces the following output files for each taxon:

### SHAP Analysis
- `shap_summary_{taxon}.pdf`: Summary plot of SHAP values
- `shap_feature_importance_{taxon}.pdf`: Feature importance based on SHAP values
- `shap_dependence_{taxon}_{feature}.pdf`: Dependence plots for top features

### Mixed Models
- `mixed_model_results_{taxon}.txt`: Statistical results from the mixed model
- `coefficient_plot_{taxon}.pdf`: Visualization of model coefficients

### Summary Files
- `shap_feature_importance_all_microbes.csv`: Combined SHAP results across all analyzed taxa
- `model_summary_all_microbes.csv`: Combined mixed model results for all taxa

## Interpreting Results

### SHAP Summary Plots
- Show the impact of each feature on the model output
- Features are ordered by importance (top = most important)
- Each point represents a sample
- Position on x-axis shows the impact on model prediction
- Color shows the feature value (red = high, blue = low)

### Feature Importance Plots
- Rank features by their overall impact on predictions
- Based on the mean absolute SHAP value across all samples
- Provides a clearer ranking than traditional feature importance

### Mixed Model Results
- Estimate the effect size of each predictor
- Account for random effects (e.g., within-subject correlations)
- Provide statistical significance (p-values) for each effect
- Validate the relationships identified by SHAP analysis

## Python API Usage

For more flexibility, you can use the Python API:

```python
from kraken_tools.analysis.rf_shap import run_rf_shap_analysis

results = run_rf_shap_analysis(
    abundance_file="processed_abundance.tsv",
    metadata_file="metadata.csv",
    output_dir="results/rf_shap/",
    target_taxa=["Bacteroides.fragilis", "Escherichia.coli"],
    predictors=["Treatment", "TimePoint", "Age"],
    random_effects=["Subject"],
    transform="clr"
)

# Access results programmatically
for taxon, taxon_results in results.items():
    print(f"\nResults for {taxon}:")
    print(f"Top predictors (SHAP): {taxon_results['top_predictors']}")
    print(f"Mixed model p-values: {taxon_results['pvalues']}")
```

## Example

Here's an example analyzing the factors affecting *Clostridioides difficile* abundance in a hospital microbiome study:

```bash
kraken-tools rf-shap \
    --abundance-file hospital_microbiome.tsv \
    --sample-key patient_data.csv \
    --output-dir results/hospital_study/ \
    --target-taxa "Clostridioides.difficile" \
    --predictors "AntibioticUse,Age,LengthOfStay,DietType,RoomType" \
    --random-effects "PatientID,HospitalWard" \
    --transform "clr"
```

This might reveal that AntibioticUse and LengthOfStay are the strongest predictors of *C. difficile* abundance, with particular antibiotics showing the most pronounced effects.

## Methodology

The RF-SHAP analysis process works as follows:

1. **Data Preparation**:
   - Transform abundance data (CLR, TSS, or other method)
   - Extract target taxa
   - Encode categorical variables
   - Split data into training and testing sets

2. **Random Forest Modeling**:
   - Train a Random Forest regressor for each taxon
   - Use the model to predict taxon abundance from metadata
   - Evaluate model performance on test data

3. **SHAP Analysis**:
   - Calculate SHAP values for each feature and sample
   - Generate summary plots showing feature importance
   - Create dependence plots for top features

4. **Mixed Effects Modeling**:
   - Formulate mixed models with fixed and random effects
   - Fit models using different optimization methods
   - Calculate p-values and confidence intervals for coefficients

## Tips for Effective RF-SHAP Analysis

1. **Choosing Target Taxa**: Focus on taxa of biological interest or those with sufficient abundance
2. **Sample Size**: Ensure adequate samples for reliable model building (~50+ recommended)
3. **Predictors**: Include relevant variables but avoid too many predictors relative to sample size
4. **Random Effects**: Include subject IDs or other grouping variables when you have repeated measures
5. **Transformations**: Consider the appropriate transformation for your data type
6. **Model Validation**: Check model performance metrics before interpreting SHAP values
7. **Mixed Models**: If mixed models fail to converge, try simpler models or different optimizers
