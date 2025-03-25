# kraken_tools/analysis/tsne.py
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors
import logging
import traceback

from kraken_tools.utils.file_utils import check_file_exists_with_logger
from kraken_tools.analysis.permanova import transform_abundance_data

def get_organism_colors(taxa_list):
    """
    Generate a distinct color for each taxon.
    
    Args:
        taxa_list: List of taxa names
        
    Returns:
        Dictionary mapping taxa to colors
    """
    # Define a set of distinct colors
    base_colors = [
        '#e41a1c',  # red
        '#377eb8',  # blue
        '#4daf4a',  # green
        '#984ea3',  # purple
        '#ff7f00',  # orange
        '#ffff33',  # yellow
        '#a65628',  # brown
        '#f781bf',  # pink
        '#999999',  # grey
        '#8dd3c7',  # teal
        '#bebada'   # light purple
    ]
    
    # If we have more taxa than colors, we'll cycle through the colors
    taxa_colors = {}
    for i, taxon in enumerate(taxa_list):
        taxa_colors[taxon] = base_colors[i % len(base_colors)]
    
    return taxa_colors

def create_custom_cmap(base_color):
    """
    Create a custom colormap that goes from white to the base color.
    
    Args:
        base_color: Base color (hex code)
        
    Returns:
        Custom colormap
    """
    return mcolors.LinearSegmentedColormap.from_list('custom_cmap', ['#ffffff', base_color])

def generate_tsne_embeddings(abundance_df, perplexity=30, n_iter=1000, random_state=42, logger=None):
    """
    Generate t-SNE embeddings from abundance data.
    
    Args:
        abundance_df: DataFrame with abundance data
        perplexity: Perplexity parameter for t-SNE
        n_iter: Number of iterations for t-SNE
        random_state: Random seed for reproducibility
        logger: Logger instance
        
    Returns:
        DataFrame with t-SNE coordinates
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    logger.info(f"Generating t-SNE embeddings with perplexity={perplexity}, n_iter={n_iter}")
    
    try:
        # Transpose to get samples as rows
        abundance_matrix = abundance_df.T
        
        # Perform t-SNE
        tsne = TSNE(
            n_components=2,
            perplexity=perplexity,
            n_iter=n_iter,
            random_state=random_state
        )
        
        tsne_results = tsne.fit_transform(abundance_matrix)
        
        # Create a DataFrame with t-SNE results
        tsne_df = pd.DataFrame(
            data=tsne_results,
            columns=['t-SNE1', 't-SNE2'],
            index=abundance_matrix.index
        )
        
        logger.info(f"t-SNE embedding completed successfully for {len(tsne_df)} samples")
        return tsne_df
        
    except Exception as e:
        logger.error(f"Error generating t-SNE embeddings: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def create_taxa_plots(tsne_df, abundance_df, target_taxa, output_dir, logger=None):
    """
    Create t-SNE plots colored by taxa abundance.
    
    Args:
        tsne_df: DataFrame with t-SNE coordinates
        abundance_df: DataFrame with abundance data
        target_taxa: List of taxa to visualize
        output_dir: Directory to save outputs
        logger: Logger instance
        
    Returns:
        List of output files
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    # Transpose abundance data for sample-wise access
    abundance_matrix = abundance_df.T
    
    # Get colors for each taxon
    organism_colors = get_organism_colors(target_taxa)
    
    # Create taxa plots directory
    taxa_dir = os.path.join(output_dir, "taxa_plots")
    os.makedirs(taxa_dir, exist_ok=True)
    
    output_files = []
    
    # Create a multi-page PDF
    pdf_path = os.path.join(output_dir, "tsne_abundance_plots.pdf")
    with PdfPages(pdf_path) as pdf:
        # Filter to include only taxa that exist in the dataset
        valid_taxa = [taxon for taxon in target_taxa if taxon in abundance_df.index]
        
        if not valid_taxa:
            logger.warning("None of the specified target taxa found in the dataset")
            return output_files
        
        # Create plots for each taxon
        for taxon in valid_taxa:
            try:
                logger.info(f"Creating t-SNE plot for {taxon}")
                
                # Create figure
                plt.figure(figsize=(10, 8))
                
                # Get base color for this taxon
                base_color = organism_colors[taxon]
                
                # Create custom colormap from white to the base color
                custom_cmap = create_custom_cmap(base_color)
                
                # Create a scatter plot with t-SNE coordinates
                scatter = plt.scatter(
                    tsne_df['t-SNE1'], 
                    tsne_df['t-SNE2'],
                    c=abundance_matrix[taxon],
                    cmap=custom_cmap,
                    alpha=0.9,
                    s=50
                )
                
                # Add colorbar
                cbar = plt.colorbar(scatter)
                cbar.set_label(f'{taxon} abundance')
                
                # Add labels and title
                plt.xlabel('t-SNE Component 1')
                plt.ylabel('t-SNE Component 2')
                plt.title(f't-SNE plot colored by {taxon} abundance')
                
                # Add grid
                plt.grid(alpha=0.3)
                
                # Save as PNG
                taxon_safe = taxon.replace(".", "_")
                png_path = os.path.join(taxa_dir, f"tsne_{taxon_safe}.png")
                plt.savefig(png_path, dpi=300, bbox_inches="tight")
                output_files.append(png_path)
                
                # Add to PDF
                pdf.savefig()
                
                plt.close()
                
            except Exception as e:
                logger.error(f"Error creating plot for {taxon}: {str(e)}")
        
        # Create grid visualization
        try:
            logger.info("Creating grid visualization of taxa plots")
            create_taxa_grid(tsne_df, abundance_matrix, valid_taxa, organism_colors, pdf, output_dir, logger)
        except Exception as e:
            logger.error(f"Error creating grid visualization: {str(e)}")
    
    output_files.append(pdf_path)
    logger.info(f"Taxa plots saved to {pdf_path}")
    return output_files

def create_taxa_grid(tsne_df, abundance_matrix, valid_taxa, organism_colors, pdf, output_dir, logger=None):
    """
    Create grid visualization of multiple taxa.
    
    Args:
        tsne_df: DataFrame with t-SNE coordinates
        abundance_matrix: Transposed abundance data
        valid_taxa: List of valid taxa to visualize
        organism_colors: Dictionary mapping taxa to colors
        pdf: PdfPages object for multi-page PDF
        output_dir: Directory to save outputs
        logger: Logger instance
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    # Calculate grid dimensions
    max_cols = 3
    n_taxa = len(valid_taxa)
    n_cols = min(max_cols, n_taxa)
    n_rows = (n_taxa + n_cols - 1) // n_cols
    
    # Create figure and axes
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*5, n_rows*4))
    if n_rows == 1 and n_cols == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # Plot each taxon
    for i, taxon in enumerate(valid_taxa):
        # Get base color for this taxon
        base_color = organism_colors[taxon]
        
        # Create custom colormap
        custom_cmap = create_custom_cmap(base_color)
        
        # Create scatter plot
        scatter = axes[i].scatter(
            tsne_df['t-SNE1'], 
            tsne_df['t-SNE2'],
            c=abundance_matrix[taxon],
            cmap=custom_cmap,
            alpha=0.9,
            s=30
        )
        
        # Add colorbar
        cbar = fig.colorbar(scatter, ax=axes[i])
        cbar.set_label('Abundance')
        
        # Add labels and title
        axes[i].set_xlabel('t-SNE1')
        axes[i].set_ylabel('t-SNE2')
        axes[i].set_title(taxon)
        axes[i].grid(alpha=0.3)
    
    # Hide any unused subplots
    for j in range(n_taxa, len(axes)):
        axes[j].set_visible(False)
    
    plt.tight_layout()
    
    # Save as PNG
    grid_path = os.path.join(output_dir, "tsne_all_taxa_grid.png")
    plt.savefig(grid_path, dpi=300, bbox_inches="tight")
    
    # Save as PDF
    pdf_path = os.path.join(output_dir, "tsne_all_taxa_grid.pdf")
    plt.savefig(pdf_path, bbox_inches="tight")
    
    # Add to multi-page PDF
    pdf.savefig(fig)
    
    plt.close(fig)
    
    logger.info(f"Grid visualization saved to {grid_path} and {pdf_path}")

def create_metadata_plots(tsne_df, metadata_df, categorical_vars, output_dir, logger=None):
    """
    Create t-SNE plots colored by metadata variables.
    
    Args:
        tsne_df: DataFrame with t-SNE coordinates
        metadata_df: DataFrame with metadata
        categorical_vars: List of categorical variables to visualize
        output_dir: Directory to save outputs
        logger: Logger instance
        
    Returns:
        List of output files
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    # Create metadata plots directory
    metadata_dir = os.path.join(output_dir, "metadata_plots")
    os.makedirs(metadata_dir, exist_ok=True)
    
    output_files = []
    
    # Merge t-SNE coordinates with metadata
    try:
        logger.info("Merging t-SNE coordinates with metadata")
        tsne_with_meta = pd.merge(
            tsne_df, 
            metadata_df, 
            left_index=True, 
            right_index=True,
            how='inner'
        )
        
        logger.info(f"Found {len(tsne_with_meta)} samples with both t-SNE coordinates and metadata")
        
        # Create metadata PDF
        pdf_path = os.path.join(output_dir, "tsne_metadata_plots.pdf")
        with PdfPages(pdf_path) as pdf:
            # Filter to valid categorical variables
            valid_vars = [var for var in categorical_vars if var in tsne_with_meta.columns]
            
            if not valid_vars:
                logger.warning("None of the specified categorical variables found in the metadata")
                return output_files
            
            # Create plots for each variable
            for var in valid_vars:
                try:
                    logger.info(f"Creating t-SNE plot for metadata variable {var}")
                    
                    # Create figure
                    plt.figure(figsize=(12, 10))
                    
                    # Create categorical plot
                    sns.scatterplot(
                        data=tsne_with_meta,
                        x='t-SNE1',
                        y='t-SNE2',
                        hue=var,
                        palette='Set1',
                        s=100,
                        alpha=0.8
                    )
                    
                    plt.title(f't-SNE plot colored by {var}')
                    plt.xlabel('t-SNE Component 1')
                    plt.ylabel('t-SNE Component 2')
                    plt.grid(alpha=0.3)
                    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                    
                    plt.tight_layout()
                    
                    # Save as PNG
                    var_safe = var.replace(" ", "_")
                    png_path = os.path.join(metadata_dir, f"tsne_meta_{var_safe}.png")
                    plt.savefig(png_path, dpi=300, bbox_inches="tight")
                    output_files.append(png_path)
                    
                    # Add to PDF
                    pdf.savefig()
                    
                    plt.close()
                    
                except Exception as e:
                    logger.error(f"Error creating plot for {var}: {str(e)}")
            
            # Create grid visualization if multiple variables
            if len(valid_vars) > 1:
                try:
                    logger.info("Creating grid visualization of metadata variables")
                    create_metadata_grid(tsne_with_meta, valid_vars, pdf, output_dir, logger)
                except Exception as e:
                    logger.error(f"Error creating metadata grid: {str(e)}")
        
        output_files.append(pdf_path)
        logger.info(f"Metadata plots saved to {pdf_path}")
        
    except Exception as e:
        logger.error(f"Error creating metadata plots: {str(e)}")
        logger.error(traceback.format_exc())
    
    return output_files

def create_metadata_grid(tsne_with_meta, valid_vars, pdf, output_dir, logger=None):
    """
    Create grid visualization of multiple metadata variables.
    
    Args:
        tsne_with_meta: DataFrame with t-SNE coordinates and metadata
        valid_vars: List of valid metadata variables
        pdf: PdfPages object for multi-page PDF
        output_dir: Directory to save outputs
        logger: Logger instance
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    # Calculate grid dimensions
    max_cols = 2
    n_vars = len(valid_vars)
    n_cols = min(max_cols, n_vars)
    n_rows = (n_vars + n_cols - 1) // n_cols
    
    # Create figure and axes
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*6, n_rows*5))
    if n_rows == 1 and n_cols == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # Plot each variable
    for i, var in enumerate(valid_vars):
        # Use Seaborn for categorical plots
        sns.scatterplot(
            data=tsne_with_meta,
            x='t-SNE1',
            y='t-SNE2',
            hue=var,
            palette='Set1',
            s=80,
            alpha=0.8,
            ax=axes[i]
        )
        
        axes[i].set_title(f'{var}')
        axes[i].set_xlabel('t-SNE1')
        axes[i].set_ylabel('t-SNE2')
        axes[i].grid(alpha=0.3)
        
        # Move legend outside for better visibility
        axes[i].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Hide any unused subplots
    for j in range(n_vars, len(axes)):
        axes[j].set_visible(False)
    
    plt.tight_layout()
    
    # Save as PNG
    grid_path = os.path.join(output_dir, "tsne_all_metadata_grid.png")
    plt.savefig(grid_path, dpi=300, bbox_inches="tight")
    
    # Save as PDF
    pdf_path = os.path.join(output_dir, "tsne_all_metadata_grid.pdf")
    plt.savefig(pdf_path, bbox_inches="tight")
    
    # Add to multi-page PDF
    pdf.savefig(fig)
    
    plt.close(fig)
    
    logger.info(f"Metadata grid visualization saved to {grid_path} and {pdf_path}")

def create_combined_overview(tsne_df, abundance_df, metadata_df, target_taxa, categorical_vars, output_dir, logger=None):
    """
    Create a combined overview visualization with key taxa and metadata.
    
    Args:
        tsne_df: DataFrame with t-SNE coordinates
        abundance_df: DataFrame with abundance data
        metadata_df: DataFrame with metadata
        target_taxa: List of taxa to visualize
        categorical_vars: List of categorical variables to visualize
        output_dir: Directory to save outputs
        logger: Logger instance
        
    Returns:
        Path to the output file
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    try:
        logger.info("Creating combined overview visualization")
        
        # Transpose abundance data
        abundance_matrix = abundance_df.T
        
        # Merge with metadata
        merged_data = pd.merge(
            tsne_df, 
            metadata_df, 
            left_index=True, 
            right_index=True,
            how='inner'
        )
        
        # Filter to valid taxa and variables
        valid_taxa = [taxon for taxon in target_taxa if taxon in abundance_df.index]
        valid_vars = [var for var in categorical_vars if var in merged_data.columns]
        
        if not valid_taxa or not valid_vars:
            logger.warning("Not enough valid taxa or metadata variables for combined visualization")
            return None
        
        # Select top taxa by abundance
        top_taxa = valid_taxa[:min(3, len(valid_taxa))]
        
        # Select key metadata variables
        key_vars = valid_vars[:min(3, len(valid_vars))]
        
        # Calculate grid dimensions
        n_plots = len(top_taxa) + len(key_vars)
        n_cols = min(3, n_plots)
        n_rows = (n_plots + n_cols - 1) // n_cols
        
        # Create figure and axes
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*5, n_rows*4))
        if n_rows == 1 and n_cols == 1:
            axes = np.array([axes])
        axes = axes.flatten()
        
        # Track current axis index
        ax_idx = 0
        
        # Plot top taxa
        organism_colors = get_organism_colors(valid_taxa)
        for taxon in top_taxa:
            # Get base color for this taxon
            base_color = organism_colors[taxon]
            
            # Create custom colormap
            custom_cmap = create_custom_cmap(base_color)
            
            # Create scatter plot
            scatter = axes[ax_idx].scatter(
                tsne_df['t-SNE1'], 
                tsne_df['t-SNE2'],
                c=abundance_matrix[taxon],
                cmap=custom_cmap,
                alpha=0.9,
                s=30
            )
            
            # Add colorbar
            cbar = fig.colorbar(scatter, ax=axes[ax_idx])
            cbar.set_label('Abundance')
            
            # Add labels and title
            axes[ax_idx].set_xlabel('t-SNE1')
            axes[ax_idx].set_ylabel('t-SNE2')
            axes[ax_idx].set_title(f"Taxa: {taxon}")
            axes[ax_idx].grid(alpha=0.3)
            
            ax_idx += 1
        
        # Plot key metadata variables
        for var in key_vars:
            # Use Seaborn for categorical plots
            sns.scatterplot(
                data=merged_data,
                x='t-SNE1',
                y='t-SNE2',
                hue=var,
                palette='Set1',
                s=50,
                alpha=0.8,
                ax=axes[ax_idx]
            )
            
            axes[ax_idx].set_title(f"Metadata: {var}")
            axes[ax_idx].set_xlabel('t-SNE1')
            axes[ax_idx].set_ylabel('t-SNE2')
            axes[ax_idx].grid(alpha=0.3)
            axes[ax_idx].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            
            ax_idx += 1
        
        # Hide any unused subplots
        for j in range(ax_idx, len(axes)):
            axes[j].set_visible(False)
        
        plt.suptitle("t-SNE Overview: Key Taxa and Metadata Variables", fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.97])  # Make room for suptitle
        
        # Save as PNG
        overview_path = os.path.join(output_dir, "tsne_combined_overview.png")
        plt.savefig(overview_path, dpi=300, bbox_inches="tight")
        
        # Save as PDF
        pdf_path = os.path.join(output_dir, "tsne_combined_overview.pdf")
        plt.savefig(pdf_path, bbox_inches="tight")
        
        plt.close(fig)
        
        logger.info(f"Combined overview visualization saved to {overview_path}")
        return overview_path
        
    except Exception as e:
        logger.error(f"Error creating combined overview: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def run_tsne_analysis(abundance_file, metadata_file, output_dir, target_taxa=None, 
                    categorical_vars=None, group_col=None, transform="clr",
                    perplexity=30, n_iter=1000, random_state=42, log_file=None):
    """
    Run t-SNE analysis on microbiome abundance data.
    
    Args:
        abundance_file: Path to abundance file
        metadata_file: Path to metadata file
        output_dir: Directory to save outputs
        target_taxa: List of taxa to visualize (comma-separated string or list)
        categorical_vars: List of categorical variables to visualize (comma-separated string or list)
        group_col: Primary grouping variable
        transform: Transformation to apply to abundance data
        perplexity: Perplexity parameter for t-SNE
        n_iter: Number of iterations for t-SNE
        random_state: Random seed for reproducibility
        log_file: Path to log file
        
    Returns:
        Directory with output files
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
    
    logger.info("Starting t-SNE analysis")
    
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
        
        # Process target taxa
        if target_taxa:
            if isinstance(target_taxa, str):
                target_taxa = [t.strip() for t in target_taxa.split(',')]
            
            # Verify taxa exist in the data
            existing_taxa = [t for t in target_taxa if t in abundance_df.index]
            if not existing_taxa:
                logger.warning(f"None of the specified target taxa found in data: {target_taxa}")
                # Use top taxa by abundance instead
                mean_abundance = abundance_df.mean(axis=1)
                target_taxa = mean_abundance.sort_values(ascending=False).head(10).index.tolist()
                logger.info(f"Using top 10 taxa by abundance instead: {target_taxa}")
            else:
                target_taxa = existing_taxa
                logger.info(f"Found {len(target_taxa)} of the specified taxa in the data")
        else:
            # If no taxa specified, use top taxa by abundance
            mean_abundance = abundance_df.mean(axis=1)
            target_taxa = mean_abundance.sort_values(ascending=False).head(10).index.tolist()
            logger.info(f"Using top 10 taxa by abundance: {target_taxa}")
        
        # Process categorical variables
        if categorical_vars:
            if isinstance(categorical_vars, str):
                categorical_vars = [v.strip() for v in categorical_vars.split(',')]
            
            # Verify variables exist in the data
            existing_vars = [v for v in categorical_vars if v in metadata_df.columns]
            if not existing_vars:
                logger.warning(f"None of the specified categorical variables found in metadata: {categorical_vars}")
                # Find categorical columns
                categorical_cols = metadata_df.select_dtypes(include=['object', 'category']).columns.tolist()
                if categorical_cols:
                    categorical_vars = categorical_cols[:min(5, len(categorical_cols))]
                    logger.info(f"Using {len(categorical_vars)} categorical variables instead: {categorical_vars}")
                else:
                    categorical_vars = []
            else:
                categorical_vars = existing_vars
                logger.info(f"Found {len(categorical_vars)} of the specified variables in the metadata")
        else:
            # Find categorical columns
            categorical_cols = metadata_df.select_dtypes(include=['object', 'category']).columns.tolist()
            if categorical_cols:
                categorical_vars = categorical_cols[:min(5, len(categorical_cols))]
                logger.info(f"Using {len(categorical_vars)} categorical variables: {categorical_vars}")
            else:
                categorical_vars = []
        
        # Include group_col if specified and not already in categorical_vars
        if group_col and group_col not in categorical_vars and group_col in metadata_df.columns:
            categorical_vars.append(group_col)
            logger.info(f"Added group column '{group_col}' to categorical variables")
        
        # Transform abundance data
        abundance_transformed = transform_abundance_data(abundance_df, transform, logger)
        
        # Generate t-SNE embeddings
        tsne_df = generate_tsne_embeddings(
            abundance_transformed,
            perplexity=perplexity,
            n_iter=n_iter,
            random_state=random_state,
            logger=logger
        )
        
        if tsne_df is None:
            logger.error("Failed to generate t-SNE embeddings")
            return None
        
        # Save t-SNE coordinates
        tsne_coords_file = os.path.join(output_dir, "tsne_coordinates.tsv")
        tsne_df.to_csv(tsne_coords_file, sep='\t')
        logger.info(f"t-SNE coordinates saved to {tsne_coords_file}")
        
        # Create taxa plots
        taxa_files = create_taxa_plots(
            tsne_df,
            abundance_transformed,
            target_taxa,
            output_dir,
            logger
        )
        
        # Create metadata plots
        if categorical_vars:
            metadata_files = create_metadata_plots(
                tsne_df,
                metadata_df,
                categorical_vars,
                output_dir,
                logger
            )
        
        # Create combined overview
        overview_file = create_combined_overview(
            tsne_df,
            abundance_transformed,
            metadata_df,
            target_taxa,
            categorical_vars,
            output_dir,
            logger
        )
        
        logger.info("t-SNE analysis completed successfully")
        return output_dir
        
    except Exception as e:
        logger.error(f"Error in t-SNE analysis: {str(e)}")
        logger.error(traceback.format_exc())
        return None

# If run standalone
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Run t-SNE analysis on microbiome data")
    parser.add_argument("--abundance-file", required=True, help="Path to abundance file")
    parser.add_argument("--metadata-file", required=True, help="Path to metadata file")
    parser.add_argument("--output-dir", required=True, help="Directory to save outputs")
    parser.add_argument("--target-taxa", help="Comma-separated list of taxa to visualize")
    parser.add_argument("--categorical-vars", help="Comma-separated list of categorical variables to visualize")
    parser.add_argument("--group-col", help="Primary grouping variable")
    parser.add_argument("--transform", default="clr", choices=["clr", "hellinger", "log", "none"], help="Transformation to apply")
    parser.add_argument("--perplexity", type=int, default=30, help="Perplexity parameter for t-SNE")
    parser.add_argument("--n-iter", type=int, default=1000, help="Number of iterations for t-SNE")
    parser.add_argument("--log-file", help="Path to log file")
    
    args = parser.parse_args()
    
    # Run analysis
    run_tsne_analysis(
        abundance_file=args.abundance_file,
        metadata_file=args.metadata_file,
        output_dir=args.output_dir,
        target_taxa=args.target_taxa,
        categorical_vars=args.categorical_vars,
        group_col=args.group_col,
        transform=args.transform,
        perplexity=args.perplexity,
        n_iter=args.n_iter,
        log_file=args.log_file
    )