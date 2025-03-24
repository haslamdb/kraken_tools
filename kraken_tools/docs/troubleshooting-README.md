# Troubleshooting Guide

This guide addresses common issues you might encounter when using Kraken Tools and provides solutions to help you resolve them.

## Installation Issues

### Missing Dependencies

**Problem**: Error messages about missing Python packages when running `kraken-tools`.

**Solution**:
```bash
# Make sure your environment is activated
conda activate kraken_tools

# Install missing dependencies
conda install -c conda-forge -c bioconda PACKAGE_NAME
# or
pip install PACKAGE_NAME
```

### External Tool Installation

**Problem**: "Command not found" errors for Kraken2, Bracken, or KneadData.

**Solution**:
```bash
# Install external tools with conda
conda install -c bioconda kraken2 bracken kneaddata

# Verify installation
kraken2 --version
est_abundance.py --help
kneaddata --version
```

## Input File Issues

### File Format Problems

**Problem**: "No valid Kraken/Bracken files found for any samples."

**Solution**:
1. Check file names and extensions
2. Use the `--list-files` command to see what files are detected:
   ```bash
   kraken-tools list-files --kreport-dir kraken_reports/ --bracken-dir bracken_files/
   ```
3. Ensure files follow naming conventions:
   - Kraken reports: `SAMPLE_ID.kreport`, `SAMPLE_ID.kreport.txt`
   - Bracken files: `SAMPLE_ID.bracken`, `SAMPLE_ID_abundance.txt`

### Sample Key Issues

**Problem**: Sample identifiers in metadata don't match file names.

**Solution**:
1. Make sure your sample key's sample identifier column (usually "SampleName") matches the sample IDs used in your Kraken/Bracken files
2. If needed, rename files or update the sample key to ensure consistency
3. Verify column names in the sample key match those expected by the tools

## Processing Issues

### Memory Errors

**Problem**: Out of memory errors when processing large datasets.

**Solution**:
1. Use the `--max-memory` option to limit memory usage
2. Process data in smaller batches
3. Increase system swap space
4. Use `--use-parallel` with `--max-parallel` to control resource usage:
   ```bash
   kraken-tools process --kreport-dir kraken_reports/ --use-parallel --max-parallel 4
   ```

### Long Processing Times

**Problem**: Analysis takes too long to complete.

**Solution**:
1. Use parallel processing with `--use-parallel` and `--threads-per-sample` options
2. Filter data earlier in the pipeline with `--min-abundance` and `--min-prevalence` options
3. Process only a subset of samples for initial testing

## Analysis Issues

### No Significant Results

**Problem**: "No significant taxa found after FDR correction."

**Solution**:
1. Try a less stringent significance threshold:
   ```bash
   kraken-tools diff-abundance --abundance-file abundance.tsv --sample-key metadata.csv --output-dir results/ --fdr-cutoff 0.1
   ```
2. Adjust filtering parameters to include more taxa:
   ```bash
   kraken-tools process --kreport-dir kraken_reports/ --bracken-dir bracken_files/ --min-abundance 0.001 --min-prevalence 0.05
   ```
3. Check sample size - insufficient samples may limit statistical power

### PERMANOVA Errors

**Problem**: "Error in PERMANOVA analysis: singular matrix."

**Solution**:
1. Ensure you have enough samples in each group (minimum 3-5 recommended)
2. Remove groups with too few samples using `--min-group-size` option
3. Check for high correlation between variables

### Mixed Model Convergence Issues

**Problem**: "Mixed model failed to converge."

**Solution**:
1. Try a different optimizer:
   ```bash
   kraken-tools glmm --abundance-file abundance.tsv --formula "Count ~ Group + (1|Subject)" --optimizer "lbfgs"
   ```
2. Simplify the model by removing complex terms
3. Ensure you have sufficient samples for the model complexity
4. Check for collinearity among predictors

## Transformation Issues

**Problem**: "Error in CLR transformation: log(0) is undefined."

**Solution**:
1. Use a different transformation method:
   ```bash
   kraken-tools analyze --abundance-file abundance.tsv --transform "hellinger"
   ```
2. Add a pseudocount to avoid zeros:
   ```bash
   kraken-tools process --kreport-dir kraken_reports/ --bracken-dir bracken_files/ --pseudocount 0.5
   ```

## Visualization Issues

**Problem**: Missing or empty plots in output directory.

**Solution**:
1. Check log file for specific error messages
2. Ensure you have the required matplotlib dependencies:
   ```bash
   conda install -c conda-forge matplotlib seaborn matplotlib-venn
   ```
3. Check that your data has sufficient variation to be visualized

## Command-Specific Troubleshooting

### KneadData Issues

**Problem**: "KneadData run failed for sample."

**Solution**:
1. Check input file quality and format
2. Verify reference database paths
3. Check disk space for temporary files
4. Run with `--log-level DEBUG` for more detailed output

### Kraken2/Bracken Issues

**Problem**: "Kraken2/Bracken run failed for sample."

**Solution**:
1. Verify database paths are correct
2. Check that reference databases are built properly
3. Ensure enough memory is available for database loading
4. Run commands directly to see detailed error messages:
   ```bash
   kraken2 --db /path/to/db --output output.kraken --report output.kreport input.fastq
   ```

### Feature Selection Issues

**Problem**: "Error in feature selection: insufficient number of samples."

**Solution**:
1. Increase the number of samples in your dataset
2. Reduce the number of features by filtering or selecting key variables
3. Use a simpler model with fewer parameters

### RF-SHAP Issues

**Problem**: "Error calculating SHAP values."

**Solution**:
1. Update the SHAP library to the latest version:
   ```bash
   pip install --upgrade shap
   ```
2. Try a different machine learning model type
3. Reduce the dataset size if memory is an issue
4. Check for NaN or infinite values in your data

## Getting More Help

If you encounter issues not covered in this troubleshooting guide:

1. Check the log file for detailed error messages:
   ```bash
   cat kraken_analysis.log
   ```

2. Run with increased verbosity:
   ```bash
   kraken-tools COMMAND --log-level DEBUG
   ```

3. Open an issue on the GitHub repository with:
   - A description of the problem
   - Command used
   - Log file contents
   - Input data description (file sizes, number of samples, etc.)
   - System information (OS, memory, CPU)
