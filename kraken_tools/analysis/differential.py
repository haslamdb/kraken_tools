# kraken_tools/analysis/differential.py
import os
import numpy as np
import pandas as pd
import logging
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests
import traceback

# Create our own CLR implementation to avoid skbio dependency
def clr_transform(data_matrix):
    """
    Compute the centered log-ratio (CLR) transformation.
    
    Parameters:
    -----------
    data_matrix : numpy.ndarray or pandas.DataFrame
        The data matrix to transform, with features as columns
        
    Returns:
    --------
    numpy.ndarray or pandas.DataFrame
        The CLR-transformed data matrix
    """
    if isinstance(data_matrix, pd.DataFrame):
        # For pandas DataFrame
        # Add pseudocount to avoid log(0)
        pseudo_df = data_matrix.replace(0, np.nextafter(0, 1))
        log_data = np.log(pseudo_df)
        geometric_means = log_data.mean(axis=1)
        clr_data = log_data.subtract(geometric_means, axis=0)
        return clr_data
    else:
        # For numpy array
        # Add pseudocount to avoid log(0)
        pseudo_array = np.where(data_matrix == 0, np.nextafter(0, 1), data_matrix)
        log_data = np.log(pseudo_array)
        geometric_means = np.mean(log_data, axis=1, keepdims=True)
        return log_data - geometric_means

def aldex2_like(abundance_df, metadata_df, group_col, mc_samples=128, denom="all", filter_groups=None):
    """
    A Python implementation similar to ALDEx2 for differential abundance testing
    
    Parameters:
    -----------
    abundance_df : pandas DataFrame
        Feature table with samples as columns and features as rows
    metadata_df : pandas DataFrame
        Metadata with sample IDs as index and metadata as columns
    group_col : str
        Column name in metadata_df that contains the grouping variable
    mc_samples : int
        Number of Monte Carlo samples to generate
    denom : str
        Features to use as denominator: "all" for all features, 
        "unmapped_excluded" to exclude unmapped features
    filter_groups : list or None
        List of group names to include in the analysis. If provided, only these groups will be used.
        Must contain exactly 2 groups for ALDEx2.
        
    Returns:
    --------
    pandas DataFrame with test results
    """
    logger = logging.getLogger('kraken_analysis')
    
    # Make sure metadata and abundance data have matching samples
    shared_samples = list(set(abundance_df.columns) & set(metadata_df.index))
    if len(shared_samples) == 0:
        logger.error("No shared samples between abundance data and metadata")
        raise ValueError("No shared samples between abundance data and metadata")
    
    abundance = abundance_df[shared_samples].copy()
    metadata = metadata_df.loc[shared_samples].copy()
    
    # Handle unmapped reads based on denom parameter
    if denom == "unmapped_excluded" and "UNMAPPED" in abundance.index:
        logger.info("Excluding unmapped reads from denominator")
        abundance = abundance.drop("UNMAPPED", axis=0)
    
    # Replace zeros with small value (pseudocount)
    min_val = abundance[abundance > 0].min().min() / 2
    abundance = abundance.replace(0, min_val)
    
    # Get group information
    groups = metadata[group_col]
    unique_groups = groups.unique()
    
    # Apply group filtering if specified
    if filter_groups is not None:
        if not isinstance(filter_groups, list):
            filter_groups = [filter_groups]
        
        # Verify the specified groups exist in the data
        missing_groups = [g for g in filter_groups if g not in unique_groups]
        if missing_groups:
            logger.error(f"The following specified groups don't exist in the data: {missing_groups}")
            raise ValueError(f"Groups not found in data: {missing_groups}")
        
        # Filter metadata and abundance to only include samples from specified groups
        valid_samples = metadata[metadata[group_col].isin(filter_groups)].index
        if len(valid_samples) == 0:
            logger.error(f"No samples found for groups: {filter_groups}")
            raise ValueError(f"No samples found for groups: {filter_groups}")
        
        metadata = metadata.loc[valid_samples]
        abundance = abundance[valid_samples]
        groups = metadata[group_col]
        unique_groups = groups.unique()
        
        logger.info(f"Filtered to {len(unique_groups)} groups: {unique_groups}")
    
    # ALDEx2 requires exactly 2 groups for comparison
    if len(unique_groups) != 2:
        logger.error(f"ALDEx2 implementation requires exactly 2 groups for comparison, found {len(unique_groups)}: {unique_groups}")
        if filter_groups:
            logger.error(f"Please filter to exactly 2 groups using --filter-groups option")
        else:
            logger.error(f"Please use --filter-groups option to select 2 groups from: {unique_groups}")
        raise ValueError("ALDEx2 implementation requires exactly 2 groups for comparison")
        
    group1_samples = groups[groups == unique_groups[0]].index
    group2_samples = groups[groups == unique_groups[1]].index
    
    logger.info(f"Running ALDEx2 analysis with {len(group1_samples)} samples in group '{unique_groups[0]}' and "
                f"{len(group2_samples)} samples in group '{unique_groups[1]}'")
    
    # Results dataframe
    results = pd.DataFrame(index=abundance.index)
    results['feature'] = abundance.index
    
    # Monte Carlo sampling of CLR transformations
    all_clrs = []
    
    for i in range(mc_samples):
        # Generate Dirichlet Monte Carlo instance
        mc_instance = pd.DataFrame(index=abundance.index, columns=abundance.columns)
        
        for col in abundance.columns:
            # Add random noise according to Dirichlet distribution
            mc_instance[col] = np.random.dirichlet(abundance[col], 1)[0] * abundance[col].sum()
        
        # CLR transformation
        clr_data = clr_transform(mc_instance.T + 0.5).T
        all_clrs.append(clr_data)
    
    # Calculate effect sizes and p-values across MC instances
    effect_sizes = []
    pvals = []
    
    for clr_data in all_clrs:
        # For each feature, compare between groups
        for feature in clr_data.index:
            group1_values = clr_data.loc[feature, group1_samples]
            group2_values = clr_data.loc[feature, group2_samples]
            
            # Calculate effect size (difference of means)
            effect = group1_values.mean() - group2_values.mean()
            effect_sizes.append((feature, effect))
            
            # Welch's t-test
            t_stat, p_val = stats.ttest_ind(
                group1_values, group2_values, equal_var=False
            )
            pvals.append((feature, p_val))
    
    # Aggregate results
    effect_df = pd.DataFrame(effect_sizes, columns=['feature', 'effect'])
    pval_df = pd.DataFrame(pvals, columns=['feature', 'pval'])
    
    # Group by feature and calculate median effect and p-value
    median_effects = effect_df.groupby('feature')['effect'].median()
    median_pvals = pval_df.groupby('feature')['pval'].median()
    
    # Add to results
    results['effect_size'] = median_effects
    results['p_value'] = median_pvals
    
    # Multiple testing correction
    results['q_value'] = multipletests(results['p_value'], method='fdr_bh')[1]
    
    # Add mean abundance information
    results['mean_abundance_group1'] = abundance[group1_samples].mean(axis=1)
    results['mean_abundance_group2'] = abundance[group2_samples].mean(axis=1)
    
    return results.sort_values('q_value')

def ancom(abundance_df, metadata_df, group_col, alpha=0.05, denom="all", filter_groups=None):
    """
    ANCOM for differential abundance testing
    
    Parameters:
    -----------
    abundance_df : pandas DataFrame
        Feature table with samples as columns and features as rows
    metadata_df : pandas DataFrame
        Metadata with sample IDs as index and metadata as columns
    group_col : str
        Column name in metadata_df that contains the grouping variable
    alpha : float
        Significance level for tests
    denom : str
        Features to use as denominator: "all" for all features, 
        "unmapped_excluded" to exclude unmapped features
    filter_groups : list or None
        List of group names to include in the analysis. If provided, only these groups will be used.
        
    Returns:
    --------
    pandas DataFrame with test results
    """
    logger = logging.getLogger('kraken_analysis')
    
    # Make sure metadata and abundance data have matching samples
    shared_samples = list(set(abundance_df.columns) & set(metadata_df.index))
    if len(shared_samples) == 0:
        logger.error("No shared samples between abundance data and metadata")
        raise ValueError("No shared samples between abundance data and metadata")
    
    abundance = abundance_df[shared_samples].copy()
    metadata = metadata_df.loc[shared_samples].copy()
    
    # Handle unmapped reads based on denom parameter
    if denom == "unmapped_excluded" and "UNMAPPED" in abundance.index:
        logger.info("Excluding unmapped reads from denominator")
        abundance = abundance.drop("UNMAPPED", axis=0)
    
    # Replace zeros with small value (pseudocount)
    min_val = 0.5
    abundance = abundance.replace(0, min_val)
    
    # Get group information
    groups = metadata[group_col]
    unique_groups = groups.unique()
    
    # Apply group filtering if specified
    if filter_groups is not None:
        if not isinstance(filter_groups, list):
            filter_groups = [filter_groups]
        
        # Verify the specified groups exist in the data
        missing_groups = [g for g in filter_groups if g not in unique_groups]
        if missing_groups:
            logger.error(f"The following specified groups don't exist in the data: {missing_groups}")
            raise ValueError(f"Groups not found in data: {missing_groups}")
        
        # Filter metadata and abundance to only include samples from specified groups
        valid_samples = metadata[metadata[group_col].isin(filter_groups)].index
        if len(valid_samples) == 0:
            logger.error(f"No samples found for groups: {filter_groups}")
            raise ValueError(f"No samples found for groups: {filter_groups}")
        
        metadata = metadata.loc[valid_samples]
        abundance = abundance[valid_samples]
        groups = metadata[group_col]
        unique_groups = groups.unique()
        
        logger.info(f"Filtered to {len(unique_groups)} groups: {unique_groups}")
    
    # ANCOM needs at least 2 groups
    if len(unique_groups) < 2:
        logger.error(f"ANCOM requires at least 2 groups for comparison, found {len(unique_groups)}")
        raise ValueError("ANCOM requires at least 2 groups for comparison")
    
    logger.info(f"Running ANCOM analysis with {len(unique_groups)} groups: {unique_groups}")
    
    # Results dataframe
    results = pd.DataFrame(index=abundance.index)
    results['feature'] = abundance.index
    
    # For each feature, compare against all other features
    n_features = len(abundance.index)
    W = {feature: 0 for feature in abundance.index}
    
    # For each feature (i)
    for i, feature_i in enumerate(abundance.index):
        if i % 50 == 0:  # Progress update for large datasets
            logger.info(f"ANCOM: Processing feature {i+1}/{n_features}")
            
        # Compare with all other features (j)
        for j, feature_j in enumerate(abundance.index):
            if i != j:
                # Log ratio
                log_ratio = np.log(abundance.loc[feature_i] / abundance.loc[feature_j])
                
                # Perform statistical test on log ratio between groups
                group_values = {}
                for group in groups.unique():
                    group_values[group] = log_ratio[groups == group]
                
                if len(groups.unique()) == 2:
                    # Two groups: t-test
                    group_list = list(group_values.values())
                    t_stat, p_val = stats.ttest_ind(
                        group_list[0], group_list[1], equal_var=False
                    )
                else:
                    # Multiple groups: ANOVA
                    anova_groups = []
                    for group, values in group_values.items():
                        anova_groups.append(values)
                    f_stat, p_val = stats.f_oneway(*anova_groups)
                
                # Count significant tests
                if p_val < alpha:
                    W[feature_i] += 1
    
    # Calculate W statistic and add to results
    results['W'] = [W[feature] for feature in results['feature']]
    results['W_ratio'] = results['W'] / (n_features - 1)
    
    # Add detection threshold (cutoff)
    cutoff = 0.7  # Can be adjusted based on desired sensitivity
    results['significant'] = results['W_ratio'] > cutoff
    
    # Add mean abundance information
    for group in groups.unique():
        group_samples = groups[groups == group].index
        results[f'mean_abundance_{group}'] = abundance[group_samples].mean(axis=1)
    
    return results.sort_values('W', ascending=False)

def ancom_bc(abundance_df, metadata_df, group_col, formula=None, denom="all", filter_groups=None):
    """
    ANCOM-BC for differential abundance testing
    
    Parameters:
    -----------
    abundance_df : pandas DataFrame
        Feature table with samples as columns and features as rows
    metadata_df : pandas DataFrame
        Metadata with sample IDs as index and metadata as columns
    group_col : str
        Column name in metadata_df that contains the grouping variable
    formula : str
        R-style formula for the model (e.g., "~ Group + Covariate")
        If None, will use simple one-way formula with group_col
    denom : str
        Features to use as denominator: "all" for all features, 
        "unmapped_excluded" to exclude unmapped features
    filter_groups : list or None
        List of group names to include in the analysis. If provided, only these groups will be used.
        
    Returns:
    --------
    pandas DataFrame with test results
    """
    logger = logging.getLogger('kraken_analysis')
    
    try:
        import statsmodels.api as sm
        from statsmodels.formula.api import ols
    except ImportError:
        logger.error("statsmodels is required for ANCOM-BC")
        raise ImportError("statsmodels is required for ANCOM-BC")
    
    # Make sure metadata and abundance data have matching samples
    shared_samples = list(set(abundance_df.columns) & set(metadata_df.index))
    if len(shared_samples) == 0:
        logger.error("No shared samples between abundance data and metadata")
        raise ValueError("No shared samples between abundance data and metadata")
    
    abundance = abundance_df[shared_samples].copy()
    metadata = metadata_df.loc[shared_samples].copy()
    
    # Handle unmapped reads based on denom parameter
    if denom == "unmapped_excluded" and "UNMAPPED" in abundance.index:
        logger.info("Excluding unmapped reads from denominator")
        abundance = abundance.drop("UNMAPPED", axis=0)
    
    # Replace zeros with small value (pseudocount)
    min_val = 0.5
    abundance = abundance.replace(0, min_val)
    
    # Get group information
    groups = metadata[group_col]
    unique_groups = groups.unique()
    
    # Apply group filtering if specified
    if filter_groups is not None:
        if not isinstance(filter_groups, list):
            filter_groups = [filter_groups]
        
        # Verify the specified groups exist in the data
        missing_groups = [g for g in filter_groups if g not in unique_groups]
        if missing_groups:
            logger.error(f"The following specified groups don't exist in the data: {missing_groups}")
            raise ValueError(f"Groups not found in data: {missing_groups}")
        
        # Filter metadata and abundance to only include samples from specified groups
        valid_samples = metadata[metadata[group_col].isin(filter_groups)].index
        if len(valid_samples) == 0:
            logger.error(f"No samples found for groups: {filter_groups}")
            raise ValueError(f"No samples found for groups: {filter_groups}")
        
        metadata = metadata.loc[valid_samples]
        abundance = abundance[valid_samples]
        groups = metadata[group_col]
        unique_groups = groups.unique()
        
        logger.info(f"Filtered to {len(unique_groups)} groups: {unique_groups}")
    
    # ANCOM-BC needs at least 2 groups
    if len(unique_groups) < 2:
        logger.error(f"ANCOM-BC requires at least 2 groups for comparison, found {len(unique_groups)}")
        raise ValueError("ANCOM-BC requires at least 2 groups for comparison")
    
    logger.info(f"Running ANCOM-BC analysis with {len(unique_groups)} groups: {unique_groups}")
    
    # Results dataframe
    results = pd.DataFrame(index=abundance.index)
    results['feature'] = abundance.index
    
    # Estimate sequencing depth using ANCOM-BC procedure
    # (This is a simplified version of the bias correction step)
    
    # 1. Log transform the data
    log_abundance = np.log(abundance)
    
    # 2. Get sample-wise means (log geometric means)
    sample_means = log_abundance.mean(axis=0)
    
    # 3. Center log-ratio transform to remove compositional bias
    clr_abundance = log_abundance.sub(sample_means, axis=1)
    
    # 4. Transpose for regression (samples as rows, features as columns)
    clr_abundance_t = clr_abundance.T
    
    # Set up the model formula if not provided
    if formula is None:
        formula = f"feature ~ C({group_col})"
        logger.info(f"Using formula: {formula}")
    
    # Run models for each feature
    pvals = []
    effects = []
    
    for i, feature in enumerate(clr_abundance_t.columns):
        if i % 100 == 0:  # Progress update for large datasets
            logger.info(f"ANCOM-BC: Processing feature {i+1}/{len(clr_abundance_t.columns)}")
            
        # Create a temporary dataframe for regression
        temp_df = pd.DataFrame({
            'feature': clr_abundance_t[feature],
            **metadata
        })
        
        try:
            # Fit the model
            model = ols(formula, data=temp_df).fit()
            
            # Extract p-values for the group effect
            for term in model.pvalues.index:
                if group_col in term and term != 'Intercept':
                    pvals.append((feature, model.pvalues[term]))
                    effects.append((feature, model.params[term]))
                    break
        except Exception as e:
            logger.warning(f"Error in ANCOM-BC regression for feature {feature}: {str(e)}")
    
    # Compile results
    pval_df = pd.DataFrame(pvals, columns=['feature', 'p_value'])
    effect_df = pd.DataFrame(effects, columns=['feature', 'effect_size'])
    
    # Add to results
    results = results.merge(pval_df, on='feature', how='left')
    results = results.merge(effect_df, on='feature', how='left')
    
    # Multiple testing correction
    results['q_value'] = multipletests(results['p_value'], method='fdr_bh')[1]
    
    # Add mean abundance information
    for group in groups.unique():
        group_samples = groups[groups == group].index
        results[f'mean_abundance_{group}'] = abundance[group_samples].mean(axis=1)
    
    return results.sort_values('q_value')

def run_differential_abundance_analysis(abundance_df, metadata_df, output_dir, group_col="Group", 
                                      methods=["aldex2", "ancom", "ancom-bc"], denom="all",
                                      filter_groups=None, logger=None):
    """
    Run multiple differential abundance testing methods and compare results
    
    Parameters:
    -----------
    abundance_df : pandas DataFrame
        Feature table with samples as columns and features as rows
    metadata_df : pandas DataFrame
        Metadata with sample IDs as index and metadata as columns
    output_dir : str
        Directory to save output files
    group_col : str
        Column name in metadata_df that contains the grouping variable
    methods : list
        List of methods to run. Options: "aldex2", "ancom", "ancom-bc"
    denom : str
        Features to use as denominator: "all" for all features,
        "unmapped_excluded" to exclude unmapped features
    filter_groups : list or None
        List of group names to include in the analysis. If provided, only these groups will be used.
    logger : logging.Logger
        Logger for output
        
    Returns:
    --------
    dict with results from each method
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    # Make output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    results = {}
    
    # Get unique groups in the full dataset
    unique_groups = metadata_df[group_col].unique()
    logger.info(f"Full dataset contains {len(unique_groups)} groups: {unique_groups}")
    
    # Handle group filtering
    if filter_groups is not None:
        # Convert single string to list if needed
        if isinstance(filter_groups, str):
            filter_groups = [g.strip() for g in filter_groups.split(',')]
        
        logger.info(f"Filtering to specified groups: {filter_groups}")
        
        # Validate that all specified groups exist in the data
        missing_groups = [g for g in filter_groups if g not in unique_groups]
        if missing_groups:
            logger.error(f"The following specified groups don't exist in the data: {missing_groups}")
            logger.error(f"Available groups: {unique_groups}")
            return {}
        
        # For ALDEx2, check if we have exactly 2 groups after filtering
        if "aldex2" in methods and len(filter_groups) != 2:
            logger.warning(f"ALDEx2 requires exactly 2 groups, but {len(filter_groups)} were specified.")
            if len(unique_groups) == 2:
                logger.info(f"Consider using all groups for ALDEx2 as there are exactly 2 in the data: {unique_groups}")
            else:
                logger.info(f"Available groups: {unique_groups}")
            if "aldex2" in methods and len(methods) == 1:
                logger.error("Cannot proceed with ALDEx2 analysis without exactly 2 groups")
                return {}
    
    # Check if we have exactly two groups (required for ALDEx2)
    n_unique_groups = len(unique_groups)
    if filter_groups:
        n_filtered_groups = len(filter_groups)
    else:
        n_filtered_groups = n_unique_groups
    
    # ALDEx2-like analysis
    if "aldex2" in methods:
        if n_filtered_groups != 2 and "aldex2" in methods:
            logger.warning(f"Skipping ALDEx2 analysis: found {n_filtered_groups} groups after filtering, but ALDEx2 requires exactly 2")
        else:
            logger.info("Running ALDEx2-like analysis...")
            try:
                aldex2_results = aldex2_like(
                    abundance_df, metadata_df, group_col=group_col, denom=denom,
                    filter_groups=filter_groups
                )
                aldex2_results.to_csv(os.path.join(output_dir, "aldex2_results.csv"))
                logger.info(f"  Significant features (q < 0.05): {sum(aldex2_results['q_value'] < 0.05)}")
                
                # Volcano plot
                plt.figure(figsize=(10, 6))
                plt.scatter(
                    aldex2_results['effect_size'], 
                    -np.log10(aldex2_results['p_value']),
                    alpha=0.7
                )
                # Highlight significant features
                sig_features = aldex2_results[aldex2_results['q_value'] < 0.05]
                plt.scatter(
                    sig_features['effect_size'], 
                    -np.log10(sig_features['p_value']),
                    color='red',
                    alpha=0.7
                )
                plt.axhline(-np.log10(0.05), linestyle='--', color='gray')
                plt.axvline(0, linestyle='--', color='gray')
                plt.xlabel('Effect Size')
                plt.ylabel('-log10(p-value)')
                plt.title('ALDEx2 Volcano Plot')
                plt.savefig(os.path.join(output_dir, "aldex2_volcano.png"), dpi=300, bbox_inches='tight')
                plt.close()
                
                results['aldex2'] = aldex2_results
            except Exception as e:
                logger.error(f"Error in ALDEx2 analysis: {str(e)}")
                logger.error(traceback.format_exc())
    
    # ANCOM analysis
    if "ancom" in methods:
        logger.info("Running ANCOM analysis...")
        try:
            ancom_results = ancom(
                abundance_df, metadata_df, group_col=group_col, denom=denom,
                filter_groups=filter_groups
            )
            ancom_results.to_csv(os.path.join(output_dir, "ancom_results.csv"))
            logger.info(f"  Significant features: {sum(ancom_results['significant'])}")
            
            # Bar plot for top ANCOM features
            top_ancom = ancom_results.head(20)
            plt.figure(figsize=(12, 8))
            plt.barh(top_ancom['feature'], top_ancom['W_ratio'])
            plt.axvline(0.7, linestyle='--', color='red', label='Significance threshold')
            plt.xlabel('W ratio')
            plt.ylabel('Feature')
            plt.title('Top 20 Features by ANCOM W-ratio')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "ancom_top_features.png"), dpi=300, bbox_inches='tight')
            plt.close()
            
            results['ancom'] = ancom_results
        except Exception as e:
            logger.error(f"Error in ANCOM analysis: {str(e)}")
            logger.error(traceback.format_exc())
    
    # ANCOM-BC analysis
    if "ancom-bc" in methods:
        logger.info("Running ANCOM-BC analysis...")
        try:
            ancom_bc_results = ancom_bc(
                abundance_df, metadata_df, group_col=group_col, denom=denom,
                filter_groups=filter_groups
            )
            ancom_bc_results.to_csv(os.path.join(output_dir, "ancom_bc_results.csv"))
            logger.info(f"  Significant features (q < 0.05): {sum(ancom_bc_results['q_value'] < 0.05)}")
            results['ancom_bc'] = ancom_bc_results
        except Exception as e:
            logger.error(f"Error in ANCOM-BC analysis: {str(e)}")
            logger.error(traceback.format_exc())
    
    # Compare methods if we have more than one
    if len(results) > 1:
        logger.info("Comparing results across methods...")
        significant_features = {}
        
        if 'aldex2' in results:
            significant_features['aldex2'] = set(results['aldex2'][results['aldex2']['q_value'] < 0.05]['feature'])
            
        if 'ancom' in results:
            significant_features['ancom'] = set(results['ancom'][results['ancom']['significant']]['feature'])
            
        if 'ancom_bc' in results:
            significant_features['ancom_bc'] = set(results['ancom_bc'][results['ancom_bc']['q_value'] < 0.05]['feature'])
        
        # Log the comparison information
        comparison_log = ["Overlap between significant features:"]
        for method, features in significant_features.items():
            comparison_log.append(f"{method.upper()} significant features: {len(features)}")
        
        # Pairwise comparisons
        methods = list(significant_features.keys())
        for i in range(len(methods)):
            for j in range(i+1, len(methods)):
                method1, method2 = methods[i], methods[j]
                overlap = len(significant_features[method1].intersection(significant_features[method2]))
                comparison_log.append(f"Overlap between {method1.upper()} and {method2.upper()}: {overlap}")
        
        # Three-way comparison if applicable
        if len(methods) >= 3:
            overlap = len(significant_features[methods[0]].intersection(
                significant_features[methods[1]]).intersection(
                significant_features[methods[2]]))
            comparison_log.append(f"Overlap between all three methods: {overlap}")
        
        # Log the comparison and save to file
        for line in comparison_log:
            logger.info(line)
            
        with open(os.path.join(output_dir, "method_comparison.txt"), "w") as f:
            f.write("\n".join(comparison_log))
        
        # Optionally, create a Venn diagram if matplotlib-venn is available
        try:
            from matplotlib_venn import venn2, venn3
            
            if len(methods) == 2:
                plt.figure(figsize=(8, 6))
                venn2([significant_features[methods[0]], significant_features[methods[1]]],
                      [methods[0].upper(), methods[1].upper()])
                plt.title("Overlap of Significant Features")
                plt.savefig(os.path.join(output_dir, "venn_diagram.png"), dpi=300, bbox_inches='tight')
                plt.close()
            elif len(methods) == 3:
                plt.figure(figsize=(8, 6))
                venn3([significant_features[methods[0]], significant_features[methods[1]], significant_features[methods[2]]],
                      [methods[0].upper(), methods[1].upper(), methods[2].upper()])
                plt.title("Overlap of Significant Features")
                plt.savefig(os.path.join(output_dir, "venn_diagram.png"), dpi=300, bbox_inches='tight')
                plt.close()
        except ImportError:
            logger.info("matplotlib-venn not available; skipping Venn diagram")
    
    logger.info(f"Differential abundance analysis complete. Results saved to {output_dir}")
    return results