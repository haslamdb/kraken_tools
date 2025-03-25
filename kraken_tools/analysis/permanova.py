# kraken_tools/analysis/permanova.py
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import squareform, pdist
import skbio.stats.ordination as ordination
from skbio.stats.distance import DistanceMatrix
from matplotlib.lines import Line2D
import logging
import traceback

from kraken_tools.utils.file_utils import check_file_exists_with_logger

def fix_sample_ids(df1, df2, df1_name="DataFrame 1", df2_name="DataFrame 2", logger=None):
    """
    Check and report sample ID overlap issues between two dataframes
    
    Args:
        df1: First DataFrame
        df2: Second DataFrame
        df1_name: Name for the first DataFrame
        df2_name: Name for the second DataFrame
        logger: Logger instance
        
    Returns:
        Set of common samples
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    logger.info(f"\nSample ID consistency check between {df1_name} and {df2_name}:")
    logger.info(f"- {df1_name} samples: {len(df1.index)} ({len(df1.index.unique())} unique)")
    logger.info(f"- {df2_name} samples: {len(df2.index)} ({len(df2.index.unique())} unique)")
    
    common = set(df1.index) & set(df2.index)
    logger.info(f"- Common samples: {len(common)}")
    
    only_in_df1 = set(df1.index) - set(df2.index)
    only_in_df2 = set(df2.index) - set(df1.index)
    
    if only_in_df1:
        logger.info(f"- Samples only in {df1_name}: {len(only_in_df1)}")
        if len(only_in_df1) < 10:
            logger.info(f"  IDs: {sorted(only_in_df1)}")
    
    if only_in_df2:
        logger.info(f"- Samples only in {df2_name}: {len(only_in_df2)}")
        if len(only_in_df2) < 10:
            logger.info(f"  IDs: {sorted(only_in_df2)}")
    
    # Return common samples
    return common

def transform_abundance_data(abundance_df, transform_method="clr", logger=None):
    """
    Apply various transformations to abundance data
    
    Args:
        abundance_df: DataFrame with abundance data
        transform_method: Transformation method ('clr', 'hellinger', 'log', 'none')
        logger: Logger instance
        
    Returns:
        Transformed DataFrame
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    if transform_method.lower() == "none":
        logger.info("No transformation applied to abundance data")
        return abundance_df
    
    try:
        # Try to import microbiome_transform module if available
        try:
            import microbiome_transform as mt
            
            if transform_method.lower() == "clr":
                logger.info("Applying CLR transformation using microbiome_transform")
                return mt.clr_transform(abundance_df)
            elif transform_method.lower() == "hellinger":
                logger.info("Applying Hellinger transformation using microbiome_transform")
                return mt.hellinger_transform(abundance_df)
            
        except ImportError:
            logger.warning("microbiome_transform module not found, using built-in transformations")
            pass
        
        # Built-in implementations if module not available
        if transform_method.lower() == "clr":
            logger.info("Applying built-in CLR transformation")
            # Add small pseudocount to zeros
            df_pseudo = abundance_df.replace(0, np.nextafter(0, 1))
            # Log transform
            df_log = np.log(df_pseudo)
            # Subtract sample-wise mean (CLR transformation)
            clr_data = df_log.subtract(df_log.mean(axis=1), axis=0)
            return clr_data
            
        elif transform_method.lower() == "hellinger":
            logger.info("Applying built-in Hellinger transformation")
            # Calculate sample sums
            sample_sums = abundance_df.sum(axis=1)
            # Divide each value by its sample sum
            df_rel = abundance_df.div(sample_sums, axis=0)
            # Apply square root
            return np.sqrt(df_rel)
            
        elif transform_method.lower() == "log":
            logger.info("Applying log transformation")
            # Add pseudocount
            return np.log(abundance_df + 1)
            
        else:
            logger.warning(f"Unknown transformation method: {transform_method}. Using raw data.")
            return abundance_df
            
    except Exception as e:
        logger.error(f"Error during data transformation: {str(e)}")
        logger.error(traceback.format_exc())
        logger.warning("Using untransformed data due to transformation error")
        return abundance_df

def calculate_distance_matrix(abundance_df, distance_metric="bray", logger=None):
    """
    Calculate a distance matrix from abundance data
    
    Args:
        abundance_df: DataFrame with abundance data
        distance_metric: Distance metric to use ('bray', 'jaccard', 'euclidean')
        logger: Logger instance
        
    Returns:
        skbio.DistanceMatrix object
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    logger.info(f"Calculating distance matrix using {distance_metric} distance")
    
    try:
        if distance_metric.lower() == "bray":
            # For Bray-Curtis, use scipy's implementation
            distances = pdist(abundance_df.values, metric="braycurtis")
        elif distance_metric.lower() == "jaccard":
            # For Jaccard distance
            distances = pdist(abundance_df.values, metric="jaccard")
        elif distance_metric.lower() == "euclidean":
            # For Euclidean distance
            distances = pdist(abundance_df.values, metric="euclidean")
        else:
            logger.warning(f"Unknown distance metric: {distance_metric}. Using Euclidean distance.")
            distances = pdist(abundance_df.values, metric="euclidean")
        
        # Create a DistanceMatrix object
        distance_matrix = DistanceMatrix(distances, ids=abundance_df.index)
        
        logger.info(f"Distance matrix calculated: {distance_matrix.shape[0]} x {distance_matrix.shape[0]}")
        return distance_matrix
        
    except Exception as e:
        logger.error(f"Error calculating distance matrix: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def run_permanova_tests(distance_matrix, metadata_df, categorical_vars, 
                      min_group_size=3, permutations=999, logger=None):
    """
    Run PERMANOVA tests for each categorical variable
    
    Args:
        distance_matrix: DistanceMatrix object
        metadata_df: DataFrame with metadata
        categorical_vars: List of categorical variables to test
        min_group_size: Minimum group size to include in analysis
        permutations: Number of permutations for PERMANOVA
        logger: Logger instance
        
    Returns:
        Dictionary of PERMANOVA results
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    logger.info(f"Performing PERMANOVA analysis for {len(categorical_vars)} variables")
    permanova_results = {}
    
    for feature in categorical_vars:
        if feature not in metadata_df.columns:
            logger.warning(f"Variable '{feature}' not found in metadata. Skipping.")
            continue
        
        try:
            # Drop NAs for the feature being tested
            valid_samples = metadata_df[feature].dropna().index
            if len(valid_samples) < 5:
                logger.warning(f"Skipping {feature} - not enough valid samples (n={len(valid_samples)})")
                continue
            
            # Make sure valid_samples are all in the distance matrix
            valid_samples = [sample for sample in valid_samples if sample in distance_matrix.ids]
            if len(valid_samples) < 5:
                logger.warning(f"Skipping {feature} - not enough samples after filtering for distance matrix (n={len(valid_samples)})")
                continue
                
            # Subset distance matrix and grouping variable
            feature_dm = distance_matrix.filter(valid_samples)
            grouping = metadata_df.loc[valid_samples, feature].astype(str)
            
            # Skip if only one unique value
            if len(grouping.unique()) < 2:
                logger.warning(f"Skipping {feature} - only one unique value")
                continue
            
            # Print group information for debugging
            group_counts = grouping.value_counts()
            logger.info(f"{feature} groups: {dict(group_counts)}")
            
            # Check if we have enough groups and samples
            if len(group_counts) < 2:
                logger.warning(f"Skipping {feature} - need at least 2 groups, found {len(group_counts)}")
                continue
                
            # Check for groups with too few samples
            small_groups = [group for group, count in group_counts.items() if count < min_group_size]
            if small_groups:
                logger.warning(f"Some groups in {feature} have fewer than {min_group_size} samples: {small_groups}")
                # Filter out small groups
                keep_groups = [group for group, count in group_counts.items() if count >= min_group_size]
                valid_samples = [s for s in valid_samples if grouping[s] in keep_groups]
                
                # Re-subset distance matrix and grouping variable
                feature_dm = distance_matrix.filter(valid_samples)
                grouping = metadata_df.loc[valid_samples, feature].astype(str)
                
                # Re-check if we have enough groups
                if len(grouping.unique()) < 2:
                    logger.warning(f"Skipping {feature} - not enough valid groups after filtering small groups")
                    continue
            
            # Run PERMANOVA
            logger.info(f"Running PERMANOVA for {feature} with {permutations} permutations")
            result = skbio.stats.distance.permanova(feature_dm, grouping, permutations=permutations)
            
            # Extract results based on scikit-bio version
            test_stat = result.get('test statistic', result.get('F', 0.0))
            p_value = result.get('p-value', result.get('p', 1.0))
            
            # Calculate R² (proportion of variance explained)
            if 'test statistic' in result and 'denominator' in result:
                r_squared = result['test statistic'] / (result['test statistic'] + result['denominator'])
            elif 'R2' in result:
                r_squared = result['R2']
            else:
                # Simple approximation if other methods fail
                r_squared = test_stat / (test_stat + 1.0)
            
            # Store results
            permanova_results[feature] = {
                'test_statistic': test_stat,
                'p_value': p_value,
                'R2': r_squared,
                'sample_size': len(grouping),
                'groups': dict(group_counts)
            }
            
            logger.info(f"{feature}: p-value = {p_value:.4f}, R² = {r_squared:.4f}, n={len(grouping)}")
            
        except Exception as e:
            logger.error(f"Error in PERMANOVA analysis for {feature}: {str(e)}")
            logger.error(traceback.format_exc())
    
    return permanova_results

def run_pcoa_analysis(distance_matrix, metadata_df, permanova_results, output_dir, logger=None):
    """
    Run PCoA analysis and create visualization plots
    
    Args:
        distance_matrix: DistanceMatrix object
        metadata_df: DataFrame with metadata
        permanova_results: Dictionary of PERMANOVA results
        output_dir: Directory to save output files
        logger: Logger instance
        
    Returns:
        DataFrame with PCoA results
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    logger.info("Running PCoA analysis")
    try:
        # Run PCoA analysis
        pcoa_result = ordination.pcoa(distance_matrix, number_of_dimensions=5)
        
        # Create a DataFrame with PCoA results
        pc_cols = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5']
        n_dimensions = pcoa_result.samples.shape[1]
        actual_pc_cols = pc_cols[:min(5, n_dimensions)]
        
        pcoa_df = pd.DataFrame(
            data=pcoa_result.samples.values,
            columns=actual_pc_cols,
            index=distance_matrix.ids
        )
        
        # Save PCoA coordinates
        pcoa_file = os.path.join(output_dir, "pcoa_coordinates.tsv")
        pcoa_df.to_csv(pcoa_file, sep='\t')
        logger.info(f"PCoA coordinates saved to {pcoa_file}")
        
        # Report variance explained
        if hasattr(pcoa_result, 'proportion_explained'):
            logger.info("Proportion of variance explained by each PC:")
            for i, prop in enumerate(pcoa_result.proportion_explained[:n_dimensions]):
                logger.info(f"  PC{i+1}: {prop:.2%}")
        
        # Merge with metadata
        pcoa_with_metadata = pd.merge(
            pcoa_df, 
            metadata_df, 
            left_index=True, 
            right_index=True, 
            how='inner'
        )
        
        # Create PCoA plots for significant variables
        significant_features = [(feature, results['p_value'], results['R2']) 
                              for feature, results in permanova_results.items() 
                              if results['p_value'] < 0.05]
        
        for feature, p_value, r_squared in significant_features:
            if feature in pcoa_with_metadata.columns:
                # Create PCoA plot
                plt.figure(figsize=(12, 10))
                
                # Get unique categories
                categories = pcoa_with_metadata[feature].astype(str).unique()
                
                if len(categories) <= 10:  # Only create if we don't have too many categories
                    # Define a color palette
                    from matplotlib.cm import get_cmap
                    cmap = get_cmap('tab10' if len(categories) <= 10 else 'tab20')
                    colors = [cmap(i) for i in range(len(categories))]
                    
                    # Create scatter plot with categories
                    for i, category in enumerate(sorted(categories)):
                        mask = pcoa_with_metadata[feature].astype(str) == category
                        if sum(mask) > 0:
                            plt.scatter(
                                pcoa_with_metadata.loc[mask, 'PC1'],
                                pcoa_with_metadata.loc[mask, 'PC2'],
                                s=50, 
                                alpha=0.7,
                                color=colors[i],
                                label=f"{category} (n={sum(mask)})"
                            )
                    
                    # Add axis labels with variance explained
                    if hasattr(pcoa_result, 'proportion_explained') and len(pcoa_result.proportion_explained) >= 2:
                        plt.xlabel(f"PC1 ({pcoa_result.proportion_explained[0]:.2%} explained var.)")
                        plt.ylabel(f"PC2 ({pcoa_result.proportion_explained[1]:.2%} explained var.)")
                    else:
                        plt.xlabel("PC1")
                        plt.ylabel("PC2")
                    
                    plt.title(f"PCoA colored by {feature} (p={p_value:.4f}, R²={r_squared:.4f})")
                    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                    plt.grid(linestyle='--', alpha=0.7)
                    
                    # Save plot
                    plot_file = os.path.join(output_dir, f"pcoa_by_{feature}.png")
                    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
                    logger.info(f"PCoA plot saved to {plot_file}")
                    plt.close()
        
        return pcoa_df
    
    except Exception as e:
        logger.error(f"Error in PCoA analysis: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def create_summary_plot(permanova_results, output_dir, logger=None):
    """
    Create summary plot of PERMANOVA results
    
    Args:
        permanova_results: Dictionary of PERMANOVA results
        output_dir: Directory to save output files
        logger: Logger instance
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    try:
        # Sort features by R² (highest to lowest)
        features_sorted = sorted(
            [(feature, results['p_value'], results['R2']) 
             for feature, results in permanova_results.items()],
            key=lambda x: x[2],
            reverse=True
        )
        
        if not features_sorted:
            logger.warning("No PERMANOVA results to plot")
            return
        
        # Create summary plot
        plt.figure(figsize=(10, 6))
        
        feature_names = [f[0] for f in features_sorted]
        r2_values = [f[2] for f in features_sorted]
        p_values = [f[1] for f in features_sorted]
        
        # Create horizontal barplot
        plt.barh(
            feature_names, r2_values, 
            color=['#2c7bb6' if p < 0.05 else '#d7191c' for p in p_values]
        )
        
        plt.xlabel('R² (Proportion of Variance Explained)')
        plt.ylabel('Metadata Features')
        plt.title('PERMANOVA Results: Variance Explained by Metadata Features')
        
        # Add significance annotations
        max_r2 = max(r2_values) if r2_values else 0.1
        offset = max_r2 * 0.03
        
        for i, p in enumerate(p_values):
            significance = "**" if p < 0.01 else ("*" if p < 0.05 else "ns")
            plt.text(r2_values[i] + offset, i, significance, va='center')
        
        # Add R² values inside bars if there's enough space
        for i, r2 in enumerate(r2_values):
            if r2 > 0.03:
                plt.text(r2/2, i, f"R²={r2:.3f}", ha='center', va='center', color='white')
        
        # Add legend
        custom_lines = [
            Line2D([0], [0], color='#2c7bb6', lw=4),
            Line2D([0], [0], color='#d7191c', lw=4),
        ]
        plt.legend(custom_lines, ['p < 0.05', 'p ≥ 0.05'], loc='upper right')
        
        plt.tight_layout()
        
        # Save plot as both PNG and PDF
        summary_png = os.path.join(output_dir, "permanova_variance_explained.png")
        summary_pdf = os.path.join(output_dir, "permanova_variance_explained.pdf")
        
        plt.savefig(summary_png, dpi=300, bbox_inches='tight')
        plt.savefig(summary_pdf, bbox_inches='tight')
        
        logger.info(f"PERMANOVA summary plot saved to {summary_png} and {summary_pdf}")
        plt.close()
        
    except Exception as e:
        logger.error(f"Error creating summary plot: {str(e)}")
        logger.error(traceback.format_exc())

def save_permanova_results(permanova_results, output_dir, logger=None):
    """
    Save PERMANOVA results to CSV files
    
    Args:
        permanova_results: Dictionary of PERMANOVA results
        output_dir: Directory to save output files
        logger: Logger instance
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    try:
        # Convert results to DataFrames
        results_list = []
        for feature, results in permanova_results.items():
            results_list.append({
                'Feature': feature,
                'p_value': results['p_value'],
                'R2': results['R2'],
                'test_statistic': results['test_statistic'],
                'sample_size': results['sample_size'],
                'Significant': 'Yes' if results['p_value'] < 0.05 else 'No'
            })
        
        if not results_list:
            logger.warning("No PERMANOVA results to save")
            return
        
        # Create DataFrames sorted by p-value and R²
        results_df = pd.DataFrame(results_list)
        
        # By p-value
        results_pval = results_df.sort_values('p_value')
        results_pval_file = os.path.join(output_dir, "permanova_results_by_pvalue.csv")
        results_pval.to_csv(results_pval_file, index=False)
        
        # By R²
        results_r2 = results_df.sort_values('R2', ascending=False)
        results_r2_file = os.path.join(output_dir, "permanova_results_by_r2.csv")
        results_r2.to_csv(results_r2_file, index=False)
        
        logger.info(f"PERMANOVA results saved to {results_pval_file} and {results_r2_file}")
        
    except Exception as e:
        logger.error(f"Error saving PERMANOVA results: {str(e)}")
        logger.error(traceback.format_exc())

def run_permanova_analysis(abundance_file, metadata_file, output_dir, categorical_vars=None,
                         group_col=None, distance_metric="bray", transform="clr",
                         permutations=999, min_group_size=3, make_pcoa=True,
                         log_file=None):
    """
    Run PERMANOVA analysis on microbiome data
    
    Args:
        abundance_file: Path to abundance file
        metadata_file: Path to metadata file
        output_dir: Directory to save output files
        categorical_vars: List of categorical variables to test
        group_col: Primary grouping variable (included if not in categorical_vars)
        distance_metric: Distance metric to use
        transform: Transformation to apply to abundance data
        permutations: Number of permutations for PERMANOVA
        min_group_size: Minimum group size to include in analysis
        make_pcoa: Whether to generate PCoA plots
        log_file: Path to log file
        
    Returns:
        Dictionary of PERMANOVA results
    """
    # Setup logging
    if log_file:
        import logging.handlers
        logger = logging.getLogger('kraken_analysis')
        handler = logging.handlers.RotatingFileHandler(log_file, maxBytes=10_485_760, backupCount=5)
        handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(handler)
    else:
        logger = logging.getLogger('kraken_analysis')
    
    logger.info("Starting PERMANOVA analysis")
    
    # Validate input files
    if not check_file_exists_with_logger(abundance_file, "Abundance file", logger):
        return None
    
    if not check_file_exists_with_logger(metadata_file, "Metadata file", logger):
        return None
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Read abundance data
        abundance_df = pd.read_csv(abundance_file, sep='\t', index_col=0)
        logger.info(f"Loaded abundance data: {abundance_df.shape[0]} taxa, {abundance_df.shape[1]} samples")
        
        # Read metadata
        metadata_df = pd.read_csv(metadata_file, index_col=0)
        logger.info(f"Loaded metadata: {metadata_df.shape[0]} samples, {metadata_df.shape[1]} variables")
        
        # Ensure indices are strings for consistent joining
        abundance_df.index = abundance_df.index.astype(str)
        metadata_df.index = metadata_df.index.astype(str)
        
        # Identify categorical variables
        if categorical_vars:
            if isinstance(categorical_vars, str):
                categorical_vars = [var.strip() for var in categorical_vars.split(',')]
            logger.info(f"Using specified categorical variables: {categorical_vars}")
        else:
            # Use all object/category columns in metadata
            categorical_vars = metadata_df.select_dtypes(include=['object', 'category']).columns.tolist()
            logger.info(f"Using all categorical variables from metadata: {categorical_vars}")
        
        # Ensure group_col is in categorical_vars
        if group_col and group_col not in categorical_vars:
            categorical_vars.append(group_col)
            logger.info(f"Added group_col '{group_col}' to categorical variables")
        
        # Transform abundance data
        abundance_transformed = transform_abundance_data(abundance_df, transform, logger)
        
        # Calculate distance matrix
        distance_matrix = calculate_distance_matrix(abundance_transformed, distance_metric, logger)
        if distance_matrix is None:
            logger.error("Failed to calculate distance matrix")
            return None
        
        # Run PERMANOVA tests
        permanova_results = run_permanova_tests(
            distance_matrix, 
            metadata_df, 
            categorical_vars, 
            min_group_size, 
            permutations, 
            logger
        )
        
        if not permanova_results:
            logger.warning("No valid PERMANOVA results")
            return {}
        
        # Save results
        save_permanova_results(permanova_results, output_dir, logger)
        
        # Create summary plot
        create_summary_plot(permanova_results, output_dir, logger)
        
        # Run PCoA analysis if requested
        if make_pcoa:
            pcoa_df = run_pcoa_analysis(distance_matrix, metadata_df, permanova_results, output_dir, logger)
        
        logger.info("PERMANOVA analysis completed successfully")
        return permanova_results
        
    except Exception as e:
        logger.error(f"Error in PERMANOVA analysis: {str(e)}")
        logger.error(traceback.format_exc())
        return None
