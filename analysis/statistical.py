# kraken_tools/analysis/statistical.py
import pandas as pd
import numpy as np
import logging
import os
import traceback

from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests
import scikit_posthocs as sp

from kraken_tools.utils.file_utils import sanitize_filename

def kruskal_wallis_dunn(df_long, group_col="Group", feature_col="Taxon", 
                        abundance_col="Abundance", alpha=0.05, logger=None):
    """
    1) Kruskal-Wallis across multiple groups
    2) Adjust p-values (Benjamini–Hochberg)
    3) Dunn's post-hoc for significant features
    
    Args:
        df_long: Long-format DataFrame with samples, features, and abundances
        group_col: Column name for grouping variable
        feature_col: Column name for feature (taxon)
        abundance_col: Column name for abundance values
        alpha: Significance threshold
        logger: Logger instance for logging
        
    Returns:
        Tuple of (kw_results_df, dict_of_posthoc_dfs)
    """
    if logger:
        logger.info(f"Running Kruskal-Wallis and Dunn's (group={group_col}, feature={feature_col})")
    
    features = df_long[feature_col].unique()
    results = []
    
    for i, feat in enumerate(features):
        sub = df_long[df_long[feature_col] == feat]
        groups = sub[group_col].unique()
        
        # Skip taxa with too few groups
        if len(groups) < 2:
            continue
            
        # Get data for each group
        group_data = [sub.loc[sub[group_col] == g, abundance_col] for g in groups]
        
        # Skip if any group has less than 2 samples
        if any(len(x) < 2 for x in group_data):
            continue
            
        try:
            # Run Kruskal-Wallis test
            stat, pval = kruskal(*group_data)
            
            # Add to results
            results.append({
                feature_col: feat,
                "KW_stat": stat,
                "KW_pvalue": pval
            })
        except Exception as e:
            if logger:
                logger.warning(f"Error in Kruskal-Wallis on {feat}: {str(e)}")
    
    if not results:
        return pd.DataFrame(), {}
        
    # Create DataFrame with results
    kw_df = pd.DataFrame(results)
    
    # FDR correction
    reject, pvals_corrected, _, _ = multipletests(kw_df["KW_pvalue"], alpha=alpha, method="fdr_bh")
    kw_df["KW_padj"] = pvals_corrected
    kw_df["Reject_H0"] = reject
    
    # Dunn's post-hoc for those with Reject_H0 = True
    posthoc_results = {}
    sig_features = kw_df[kw_df["Reject_H0"]][feature_col].tolist()
    
    for feat in sig_features:
        sub = df_long[df_long[feature_col] == feat]
        try:
            # Run Dunn's test
            posthoc_df = sp.posthoc_dunn(sub, val_col=abundance_col, group_col=group_col, p_adjust="holm")
            posthoc_results[feat] = posthoc_df
        except Exception as e:
            if logger:
                logger.warning(f"Error in Dunn's post-hoc on {feat}: {str(e)}")
    
    return kw_df, posthoc_results


def run_statistical_tests(data_df, output_dir, logger, group_col="Group"):
    """
    Run statistical tests on taxonomic data and save results.
    
    Args:
        data_df: DataFrame with taxonomic data in long format
        output_dir: Directory to save results
        logger: Logger instance
        group_col: Column name for grouping variable
    """
    logger.info(f"Running statistical tests on taxonomic data (Kruskal-Wallis + Dunn) with grouping variable '{group_col}'.")
    
    try:
        # Create output subdirectory for statistical tests
        stats_dir = os.path.join(output_dir, "StatisticalTests")
        os.makedirs(stats_dir, exist_ok=True)
        
        # Check if the group column exists
        if group_col not in data_df.columns:
            logger.error(f"Group column '{group_col}' not found in data. Available columns: {data_df.columns.tolist()}")
            return
            
        # Run Kruskal-Wallis and Dunn's tests
        kw_results, dunn_results = kruskal_wallis_dunn(
            data_df,
            group_col=group_col,  # Use the user-specified grouping column
            feature_col="Taxon",
            abundance_col="Abundance",
            alpha=0.05,
            logger=logger
        )
        
        if kw_results.empty:
            logger.warning("No valid Kruskal-Wallis results. No file saved.")
            return
            
        # Save Kruskal-Wallis results
        kw_path = os.path.join(stats_dir, "kruskal_wallis_results.csv")
        kw_results.to_csv(kw_path, index=False)
        logger.info(f"Saved Kruskal-Wallis results: {kw_path}")
        
        # Report significant results
        sig_count = sum(kw_results["Reject_H0"])
        logger.info(f"{sig_count} significant taxa found after FDR correction")
        
        # Save Dunn's post-hoc test results
        if dunn_results:
            dunn_dir = os.path.join(stats_dir, "dunn_posthoc_tests")
            os.makedirs(dunn_dir, exist_ok=True)
            
            for feat, pdf in dunn_results.items():
                safe_feat = sanitize_filename(feat)
                path = os.path.join(dunn_dir, f"dunn_{safe_feat}.csv")
                pdf.to_csv(path)
                
            logger.info(f"Saved Dunn's post-hoc results for {len(dunn_results)} features to {dunn_dir}")
            
            # Create summary file with significant comparisons
            summary_rows = []
            
            for taxon, dunn_df in dunn_results.items():
                # Get all significant comparisons (p < 0.05)
                for col in dunn_df.columns:
                    for idx in dunn_df.index:
                        if col != idx and dunn_df.loc[idx, col] < 0.05:
                            summary_rows.append({
                                'Taxon': taxon,
                                'Group1': idx,
                                'Group2': col,
                                'P_value': dunn_df.loc[idx, col]
                            })
            
            if summary_rows:
                summary_df = pd.DataFrame(summary_rows)
                summary_df = summary_df.sort_values(['Taxon', 'P_value'])
                summary_path = os.path.join(stats_dir, "significant_comparisons.csv")
                summary_df.to_csv(summary_path, index=False)
                logger.info(f"Saved summary of significant comparisons: {summary_path}")
    
    except Exception as e:
        logger.error(f"Error in statistical tests: {str(e)}")
        logger.error(traceback.format_exc())
