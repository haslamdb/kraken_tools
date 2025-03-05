# kraken_tools/analysis/taxonomy.py
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import logging

def read_and_process_taxonomy(abundance_file, sample_key_df, output_dir, min_abundance=0.01, min_prevalence=0.1, logger=None):
    """
    Read taxonomic abundance data, merge with metadata, and prepare for visualization.
    
    Args:
        abundance_file: Path to abundance file
        sample_key_df: DataFrame with sample metadata
        output_dir: Directory to save outputs
        min_abundance: Minimum relative abundance threshold for inclusion (default: 0.01 = 1%)
        min_prevalence: Minimum prevalence threshold (proportion of samples) for inclusion (default: 0.1 = 10%)
        logger: Logger instance
        
    Returns:
        DataFrame with processed taxonomic data in long format
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    try:
        # Read abundance data
        df = pd.read_csv(abundance_file, sep='\t', index_col=0)
        logger.info(f"Loaded abundance data with {df.shape[0]} taxa and {df.shape[1]} samples")
        
        # Check if the data needs normalization (values > 1 suggest percentages)
        if df.values.max() > 1:
            logger.info("Converting percentages to relative abundance (0-1)")
            df = df / 100
        
        # Filter by abundance and prevalence
        if min_abundance > 0 or min_prevalence > 0:
            logger.info(f"Filtering taxa with minimum abundance {min_abundance} and minimum prevalence {min_prevalence}")
            
            # Calculate mean abundance per taxon
            mean_abundance = df.mean(axis=1)
            
            # Calculate prevalence (proportion of samples where taxon is present)
            prevalence = (df > 0).mean(axis=1)
            
            # Apply filters
            before_count = df.shape[0]
            df = df[(mean_abundance >= min_abundance) & (prevalence >= min_prevalence)]
            after_count = df.shape[0]
            
            logger.info(f"Filtered from {before_count} to {after_count} taxa")
        
        # Melt the DataFrame for easier visualization
        long_df = df.reset_index().melt(id_vars="index", var_name="SampleName", value_name="Abundance")
        long_df = long_df.rename(columns={"index": "Taxon"})
        
        # Merge with sample metadata
        merged = pd.merge(long_df, sample_key_df, on="SampleName", how="inner")
        
        if merged.empty:
            logger.error("No samples matched between abundance data and metadata")
            return None
        
        # Create basic plots
        
        # 1. Heatmap of top taxa
        plt.figure(figsize=(12, 10))
        
        # Get top N taxa by mean abundance
        top_taxa = df.mean(axis=1).sort_values(ascending=False).head(30).index
        heatmap_data = df.loc[top_taxa]
        
        # Plot heatmap with clustering
        sns.clustermap(
            heatmap_data, 
            cmap="viridis", 
            figsize=(14, 10), 
            z_score=0,  # Standardize rows
            xticklabels=True,
            yticklabels=True
        )
        
        heatmap_path = os.path.join(output_dir, "taxonomy_heatmap.svg")
        plt.savefig(heatmap_path, format="svg", dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved taxonomy heatmap: {heatmap_path}")
        
        # 2. PCA plot
        plt.figure(figsize=(10, 8))
        
        # Prepare for PCA
        X = df.T  # Transpose to get samples as rows
        X_scaled = StandardScaler().fit_transform(X)
        
        # Run PCA
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(X_scaled)
        
        # Prepare PCA DataFrame
        pca_df = pd.DataFrame({
            'PC1': pca_result[:, 0],
            'PC2': pca_result[:, 1],
            'SampleName': X.index
        })
        
        # Merge with metadata for grouping
        pca_merged = pd.merge(pca_df, sample_key_df, on="SampleName", how="inner")
        
        # Plot PCA with group coloring
        if 'Group' in pca_merged.columns:
            sns.scatterplot(
                data=pca_merged,
                x="PC1",
                y="PC2",
                hue="Group",
                palette="Set1",
                s=100
            )
        else:
            # If no Group column, just use the first categorical column
            categorical_cols = pca_merged.select_dtypes(include=['object']).columns
            if len(categorical_cols) > 1:  # First one is SampleName
                group_col = categorical_cols[1]
                sns.scatterplot(
                    data=pca_merged,
                    x="PC1",
                    y="PC2",
                    hue=group_col,
                    palette="Set1",
                    s=100
                )
            else:
                # No grouping available
                sns.scatterplot(
                    data=pca_merged,
                    x="PC1",
                    y="PC2",
                    s=100
                )
        
        plt.title(f"PCA of Taxonomic Profiles")
        plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)")
        plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)")
        plt.tight_layout()
        
        pca_path = os.path.join(output_dir, "taxonomy_pca.svg")
        plt.savefig(pca_path, format="svg", dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved PCA plot: {pca_path}")
        
        # 3. Stacked bar plot for top taxa
        plt.figure(figsize=(14, 8))
        
        # Prepare data for stacked bars
        # First, get top taxa to show (too many make the plot unreadable)
        top_taxa = df.mean(axis=1).sort_values(ascending=False).head(15).index
        
        # Subset the long dataframe to include only top taxa
        top_taxa_df = merged[merged['Taxon'].isin(top_taxa)]
        
        # Group by sample and calculate the remaining abundance for "Other"
        if 'Group' in merged.columns:
            # First calculate total abundance per sample
            total_by_sample = merged.groupby('SampleName')['Abundance'].sum().reset_index()
            total_by_sample.rename(columns={'Abundance': 'Total'}, inplace=True)
            
            # Calculate sum of top taxa per sample
            top_sum_by_sample = top_taxa_df.groupby('SampleName')['Abundance'].sum().reset_index()
            top_sum_by_sample.rename(columns={'Abundance': 'TopSum'}, inplace=True)
            
            # Merge to calculate "Other"
            sample_totals = pd.merge(total_by_sample, top_sum_by_sample, on='SampleName', how='left')
            sample_totals['Other'] = sample_totals['Total'] - sample_totals['TopSum']
            
            # Create a dataframe for "Other" category
            other_df = sample_totals[['SampleName', 'Other']].copy()
            other_df['Taxon'] = 'Other'
            other_df.rename(columns={'Other': 'Abundance'}, inplace=True)
            
            # Merge with metadata to get grouping
            other_df = pd.merge(other_df, sample_key_df, on='SampleName', how='inner')
            
            # Combine top taxa and "Other"
            plot_df = pd.concat([top_taxa_df, other_df])
            
            # Calculate mean abundance by group and taxon
            group_taxa_means = plot_df.groupby(['Group', 'Taxon'])['Abundance'].mean().reset_index()
            
            # Pivot for the stacked bar plot
            pivot_df = group_taxa_means.pivot(index='Group', columns='Taxon', values='Abundance')
            
            # Plot
            pivot_df.plot(kind='bar', stacked=True, colormap='viridis', figsize=(14, 8))
            plt.title('Mean Taxonomic Composition by Group')
            plt.xlabel('Group')
            plt.ylabel('Relative Abundance')
            plt.xticks(rotation=45)
            plt.legend(title='Taxon', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            
            bar_path = os.path.join(output_dir, "taxonomy_bars_by_group.svg")
            plt.savefig(bar_path, format="svg", dpi=300, bbox_inches="tight")
            plt.close()
            logger.info(f"Saved stacked bar plot: {bar_path}")
        
        return merged
    
    except Exception as e:
        logger.error(f"Error in taxonomic analysis: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return None

def calculate_diversity_metrics(abundance_df, output_dir, logger=None):
    """
    Calculate diversity metrics for taxonomic abundance data.
    
    Args:
        abundance_df: DataFrame with abundance data (taxa as rows, samples as columns)
        output_dir: Directory to save outputs
        logger: Logger instance
        
    Returns:
        DataFrame with diversity metrics
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    try:
        # Ensure the data is in relative abundance format (0-1)
        if abundance_df.values.max() > 1:
            abundance_df = abundance_df / 100
        
        # Create a DataFrame to store diversity metrics
        diversity_df = pd.DataFrame(index=abundance_df.columns)
        
        # 1. Richness (number of taxa present)
        diversity_df['Richness'] = (abundance_df > 0).sum(axis=0)
        
        # 2. Shannon diversity
        shannon = []
        for sample in abundance_df.columns:
            p = abundance_df[sample]
            p = p[p > 0]  # Remove zeros to avoid log(0)
            h = -np.sum(p * np.log(p))
            shannon.append(h)
        diversity_df['Shannon'] = shannon
        
        # 3. Simpson diversity (1 - sum of squared probabilities)
        simpson = []
        for sample in abundance_df.columns:
            p = abundance_df[sample]
            d = 1 - np.sum(p**2)
            simpson.append(d)
        diversity_df['Simpson'] = simpson
        
        # 4. Pielou's evenness (Shannon diversity / log(richness))
        diversity_df['Evenness'] = diversity_df['Shannon'] / np.log(diversity_df['Richness'])
        
        # Fill NaN values (e.g., if richness is 0)
        diversity_df = diversity_df.fillna(0)
        
        # Save to file
        diversity_file = os.path.join(output_dir, "diversity_metrics.tsv")
        diversity_df.to_csv(diversity_file, sep='\t')
        logger.info(f"Saved diversity metrics to {diversity_file}")
        
        # Create visualizations
        
        # Box plots of diversity metrics
        plt.figure(figsize=(14, 10))
        
        # Reshape for seaborn
        diversity_long = diversity_df.reset_index().melt(
            id_vars="index", 
            var_name="Metric", 
            value_name="Value"
        )
        diversity_long = diversity_long.rename(columns={"index": "SampleName"})
        
        # Create a 2x2 grid for the four metrics
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        metrics = ['Richness', 'Shannon', 'Simpson', 'Evenness']
        
        for i, metric in enumerate(metrics):
            ax = axes[i // 2, i % 2]
            sns.boxplot(
                data=diversity_long[diversity_long['Metric'] == metric],
                y='Value',
                ax=ax
            )
            ax.set_title(f"{metric}")
            ax.set_ylabel(metric)
        
        plt.tight_layout()
        box_path = os.path.join(output_dir, "diversity_boxplots.svg")
        plt.savefig(box_path, format="svg", dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved diversity boxplots: {box_path}")
        
        return diversity_df
    
    except Exception as e:
        logger.error(f"Error calculating diversity metrics: {str(e)}")
        return None
