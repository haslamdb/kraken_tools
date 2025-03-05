# kraken_tools/analysis/differential.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import skbio.stats.composition as composition
from statsmodels.stats.multitest import multipletests
import os
import logging
import warnings

warnings.filterwarnings('ignore')

def aldex2_like(abundance_df, metadata_df, group_col, mc_samples=128):
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
    
    # Replace zeros with small value (pseudocount)
    min_val = abundance[abundance > 0].min().min() / 2
    abundance = abundance.replace(0, min_val)
    
    # Get group information
    groups = metadata[group_col]
    unique_groups = groups.unique()
    if len(unique_groups) != 2:
        logger.error("ALDEx2 implementation only supports two groups for comparison")
        raise ValueError("This implementation only supports two groups for comparison")
        
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
        clr_data = composition.clr(mc_instance.T + 0.5).T
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

def ancom(abundance_df, metadata_df, group_col, alpha=0.05):
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
    
    # Replace zeros with small value (pseudocount)
    min_val = 0.5
    abundance = abundance.replace(0, min_val)
    
    # Get group information
    groups = metadata[group_col]
    unique_groups = groups.unique()
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
                        anova_groups.extend([values])
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

def ancom_bc(abundance_df, metadata_df, group_col, formula=None):
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
    
    # Replace zeros with small value (pseudocount)
    min_val = 0.5
    abundance = abundance.replace(0, min_val)
    
    # Get group information
    groups = metadata[group_col]
    unique_groups = groups.unique()
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
                                      methods=["aldex2", "ancom", "ancom-bc"],
                                      logger=None):
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
    
    # Check if we have exactly two groups (required for ALDEx2)
    unique_groups = metadata_df[group_col].unique()
    
    # ALDEx2-like analysis
    if "aldex2" in methods:
        if len(unique_groups) != 2 and "aldex2" in methods:
            logger.warning(f"Skipping ALDEx2 analysis: found {len(unique_groups)} groups, but ALDEx2 requires exactly 2")
        else:
            logger.info("Running ALDEx2-like analysis...")
            try:
                aldex2_results = aldex2_like(
                    abundance_df, metadata_df, group_col=group_col
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
    
    # ANCOM analysis
    if "ancom" in methods:
        logger.info("Running ANCOM analysis...")
        try:
            ancom_results = ancom(
                abundance_df, metadata_df, group_col=group_col
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
    
    # ANCOM-BC analysis
    if "ancom-bc" in methods:
        logger.info("Running ANCOM-BC analysis...")
        try:
            ancom_bc_results = ancom_bc(
                abundance_df, metadata_df, group_col=group_col
            )
            ancom_bc_results.to_csv(os.path.join(output_dir, "ancom_bc_results.csv"))
            logger.info(f"  Significant features (q < 0.05): {sum(ancom_bc_results['q_value'] < 0.05)}")
            results['ancom_bc'] = ancom_bc_results
        except Exception as e:
            logger.error(f"Error in ANCOM-BC analysis: {str(e)}")
    
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
