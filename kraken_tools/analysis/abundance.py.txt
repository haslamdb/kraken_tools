# kraken_tools/analysis/abundance.py
import os
import pandas as pd
import numpy as np
import logging

def read_bracken_abundance_file(file_path, logger=None):
    """
    Read a Bracken abundance file and extract relevant data.
    
    Args:
        file_path: Path to Bracken abundance file
        logger: Logger instance
        
    Returns:
        DataFrame with taxonomic abundance data
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    try:
        # Bracken files have a specific format with columns for:
        # - name (taxon name)
        # - taxonomy_id (NCBI taxonomy ID)
        # - taxonomy_lvl (taxonomic level code)
        # - kreads (Kraken assigned reads)
        # - kreads_percentage (Kraken percentage)
        # - breads (Bracken re-estimated reads)
        # - breads_percentage (Bracken percentage)
        
        df = pd.read_csv(file_path, sep='\t', header=0)
        
        # Log the column names for debugging
        logger.debug(f"Column names in Bracken file: {list(df.columns)}")
        
        # Standardize column names (dealing with different versions)
        column_map = {
            'name': 'taxon_name',
            'taxonomy_id': 'taxonomy_id',
            'taxonomy_lvl': 'taxonomy_level',
            'new_est_reads': 'breads', 
            'est_reads': 'breads',
            'new_est_fraction': 'breads_percentage',
            'fraction_total_reads': 'breads_percentage'
        }
        
        # Rename only the columns that exist
        rename_dict = {old: new for old, new in column_map.items() if old in df.columns}
        df = df.rename(columns=rename_dict)
        
        # If we're missing any of the essential columns, log a warning
        essential_cols = ['taxon_name', 'taxonomy_id', 'breads', 'breads_percentage']
        missing_cols = [col for col in essential_cols if col not in df.columns]
        if missing_cols:
            logger.warning(f"Missing essential columns in Bracken file: {missing_cols}")
            return None
        
        return df
    
    except Exception as e:
        logger.error(f"Error reading Bracken file {file_path}: {str(e)}")
        return None

def merge_bracken_files(bracken_files, output_file, logger=None):
    """
    Merge multiple Bracken abundance files into a single table.
    
    Args:
        bracken_files: Dict mapping sample IDs to file paths
        output_file: Path to output merged file
        logger: Logger instance
        
    Returns:
        Path to merged file, or None if merging failed
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    # Create a dict to store abundance data by taxon
    abundance_data = {}
    
    # Process each sample
    for sample_id, file_path in bracken_files.items():
        logger.info(f"Processing Bracken abundance file for sample {sample_id}")
        
        # Read the Bracken file
        bracken_df = read_bracken_abundance_file(file_path, logger)
        if bracken_df is None:
            logger.warning(f"Skipping sample {sample_id} due to file reading error")
            continue
        
        # Extract relevant columns
        taxa = bracken_df['taxon_name'].values
        abundances = bracken_df['breads_percentage'].values
        
        # Add to our abundance dict
        for taxon, abundance in zip(taxa, abundances):
            if taxon not in abundance_data:
                abundance_data[taxon] = {}
            abundance_data[taxon][sample_id] = abundance
    
    if not abundance_data:
        logger.error("No valid data found in any Bracken file")
        return None
    
    # Convert to DataFrame
    merged_df = pd.DataFrame(abundance_data).T
    
    # Fill missing values with 0
    merged_df = merged_df.fillna(0)
    
    # Save to file
    try:
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Write to TSV file
        merged_df.to_csv(output_file, sep='\t')
        logger.info(f"Merged abundance data saved to {output_file}")
        
        return output_file
    except Exception as e:
        logger.error(f"Error saving merged abundance data: {str(e)}")
        return None

def normalize_abundance(abundance_file, output_file=None, method='relabundance', logger=None):
    """
    Normalize taxonomic abundance data.
    
    Args:
        abundance_file: Path to abundance file
        output_file: Path to save normalized data (if None, modify in place)
        method: Normalization method ('relabundance', 'cpm', 'log10', etc.)
        logger: Logger instance
        
    Returns:
        Path to normalized file
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    try:
        # Read abundance data
        df = pd.read_csv(abundance_file, sep='\t', index_col=0)
        
        # Apply normalization
        if method == 'relabundance':
            # Relative abundance (0-1) - typically Bracken files are already percentages
            # Convert from percentage to relative abundance
            df = df / 100
        elif method == 'cpm':
            # Counts per million
            df = df * 10000  # Since Bracken is already in percentage (0-100)
        elif method == 'log10':
            # Log10 transformation (adding pseudocount to avoid log(0))
            df = np.log10(df / 100 + 0.000001)
        elif method == 'clr':
            # Centered log-ratio transformation
            # Add small pseudocount to avoid log(0)
            pseudo_df = df / 100 + 0.000001
            # Calculate geometric mean of each sample
            geo_mean = np.exp(np.mean(np.log(pseudo_df), axis=0))
            # CLR transformation
            df = np.log(pseudo_df / geo_mean)
        else:
            logger.warning(f"Unknown normalization method: {method}. Using raw values.")
        
        # Determine output path
        if output_file is None:
            # Modify in place
            output_file = abundance_file
        
        # Save normalized data
        df.to_csv(output_file, sep='\t')
        logger.info(f"Normalized abundance data saved to {output_file}")
        
        return output_file
    
    except Exception as e:
        logger.error(f"Error normalizing abundance data: {str(e)}")
        return abundance_file  # Return the original file on error

def process_bracken_abundance(bracken_files, output_dir, output_prefix, taxonomic_level, logger=None):
    """
    Process Bracken abundance files: merge, normalize, and prepare for analysis.
    
    Args:
        bracken_files: Dict mapping sample IDs to file paths
        output_dir: Directory for output files
        output_prefix: Prefix for output files
        taxonomic_level: Taxonomic level (e.g., 'S' for species)
        logger: Logger instance
        
    Returns:
        Path to processed abundance file
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Merge Bracken files
    merged_file = os.path.join(output_dir, f"{output_prefix}_bracken_{taxonomic_level.lower()}_merged.tsv")
    merged_result = merge_bracken_files(bracken_files, merged_file, logger)
    
    if not merged_result:
        logger.error("Failed to merge Bracken files")
        return None
    
    # Normalize to relative abundance (Bracken already gives percentages, so just convert to 0-1)
    normalized_file = os.path.join(output_dir, f"{output_prefix}_bracken_{taxonomic_level.lower()}_relabundance.tsv")
    normalized_result = normalize_abundance(merged_file, normalized_file, method='relabundance', logger=logger)
    
    if not normalized_result:
        logger.error("Failed to normalize abundance data")
        return merged_result  # Return the merged file if normalization fails
    
    logger.info(f"Bracken abundance processing completed. Final file: {normalized_result}")
    return normalized_result

def read_kraken_report(file_path, taxonomic_level='S', logger=None):
    """
    Read a Kraken2 report file (.kreport) and extract abundance for a specific taxonomic level.
    
    Args:
        file_path: Path to Kraken report file
        taxonomic_level: Taxonomic level to extract (D, P, C, O, F, G, S)
        logger: Logger instance
        
    Returns:
        DataFrame with taxonomic abundance data
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    try:
        # Define the column names for Kraken report format
        cols = ['percentage', 'clade_reads', 'taxon_reads', 'rank_code', 'taxonomy_id', 'taxon_name']
        
        # Read the file with fixed width format
        df = pd.read_csv(file_path, sep='\t', header=None, names=cols)
        
        # Map taxonomic level codes to their standard codes
        level_map = {
            'D': 'D',      # Domain
            'P': 'P',      # Phylum
            'C': 'C',      # Class
            'O': 'O',      # Order
            'F': 'F',      # Family
            'G': 'G',      # Genus
            'S': 'S'       # Species
        }
        
        # Filter by taxonomic level
        if taxonomic_level in level_map:
            level_code = level_map[taxonomic_level]
            filtered_df = df[df['rank_code'] == level_code].copy()
            
            # Clean up taxon names (Kraken uses indentation with spaces)
            filtered_df['taxon_name'] = filtered_df['taxon_name'].str.strip()
            
            return filtered_df
        else:
            logger.error(f"Invalid taxonomic level: {taxonomic_level}")
            return None
    
    except Exception as e:
        logger.error(f"Error reading Kraken report {file_path}: {str(e)}")
        return None

def merge_kraken_reports(report_files, output_file, taxonomic_level='S', logger=None):
    """
    Merge multiple Kraken report files into a single abundance table.
    
    Args:
        report_files: Dict mapping sample IDs to report file paths
        output_file: Path to output merged file
        taxonomic_level: Taxonomic level to extract
        logger: Logger instance
        
    Returns:
        Path to merged file, or None if merging failed
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    # Create a dict to store abundance data by taxon
    abundance_data = {}
    
    # Process each sample
    for sample_id, file_path in report_files.items():
        logger.info(f"Processing Kraken report for sample {sample_id}")
        
        # Read the Kraken report
        kraken_df = read_kraken_report(file_path, taxonomic_level, logger)
        if kraken_df is None:
            logger.warning(f"Skipping sample {sample_id} due to file reading error")
            continue
        
        # Extract relevant columns
        taxa = kraken_df['taxon_name'].values
        percentages = kraken_df['percentage'].values
        
        # Add to our abundance dict
        for taxon, percentage in zip(taxa, percentages):
            if taxon not in abundance_data:
                abundance_data[taxon] = {}
            abundance_data[taxon][sample_id] = percentage
    
    if not abundance_data:
        logger.error("No valid data found in any Kraken report")
        return None
    
    # Convert to DataFrame
    merged_df = pd.DataFrame(abundance_data).T
    
    # Fill missing values with 0
    merged_df = merged_df.fillna(0)
    
    # Save to file
    try:
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Write to TSV file
        merged_df.to_csv(output_file, sep='\t')
        logger.info(f"Merged Kraken report data saved to {output_file}")
        
        return output_file
    except Exception as e:
        logger.error(f"Error saving merged Kraken data: {str(e)}")
        return None

def process_kraken_reports(report_files, output_dir, output_prefix, taxonomic_level, logger=None):
    """
    Process Kraken report files: merge, normalize, and prepare for analysis.
    
    Args:
        report_files: Dict mapping sample IDs to report file paths
        output_dir: Directory for output files
        output_prefix: Prefix for output files
        taxonomic_level: Taxonomic level (e.g., 'S' for species)
        logger: Logger instance
        
    Returns:
        Path to processed abundance file
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Merge Kraken reports
    merged_file = os.path.join(output_dir, f"{output_prefix}_kraken_{taxonomic_level.lower()}_merged.tsv")
    merged_result = merge_kraken_reports(report_files, merged_file, taxonomic_level, logger)
    
    if not merged_result:
        logger.error("Failed to merge Kraken reports")
        return None
    
    # Normalize to relative abundance (Kraken percentages to 0-1)
    normalized_file = os.path.join(output_dir, f"{output_prefix}_kraken_{taxonomic_level.lower()}_relabundance.tsv")
    normalized_result = normalize_abundance(merged_file, normalized_file, method='relabundance', logger=logger)
    
    if not normalized_result:
        logger.error("Failed to normalize abundance data")
        return merged_result  # Return the merged file if normalization fails
    
    logger.info(f"Kraken report processing completed. Final file: {normalized_result}")
    return normalized_result
