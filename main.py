# kraken_tools/main.py
"""
Main module for kraken_tools package.

This module provides functions to run the full analysis pipeline
or individual components as needed.
"""
import os
import logging
import time
import traceback
import pandas as pd

from kraken_tools.logger import setup_logger, log_print
from kraken_tools.utils.sample_utils import validate_sample_key, check_input_files_exist
from kraken_tools.utils.file_utils import check_file_exists_with_logger
from kraken_tools.preprocessing.kraken_run import run_kraken
from kraken_tools.preprocessing.bracken_run import run_bracken
from kraken_tools.analysis.abundance import process_bracken_abundance, normalize_abundance
from kraken_tools.analysis.metadata import read_and_process_metadata
from kraken_tools.analysis.taxonomy import read_and_process_taxonomy
from kraken_tools.analysis.statistical import run_statistical_tests
from kraken_tools.analysis.differential import run_differential_abundance_analysis


def run_full_pipeline(
    sample_key,
    kreport_dir,
    bracken_dir,
    output_dir,
    output_prefix="ProcessedFiles",
    taxonomic_level="S",
    group_col="Group",
    skip_kraken=False,
    skip_bracken=False,
    skip_downstream=False,
    min_abundance=0.01,
    min_prevalence=0.1,
    no_interactive=False,
    log_file=None,
):
    """
    Run the full Kraken2/Bracken processing and analysis pipeline.
    
    Args:
        sample_key: Path to CSV file with sample metadata
        kreport_dir: Directory containing Kraken2 report files
        bracken_dir: Directory containing Bracken abundance files
        output_dir: Directory for output files
        output_prefix: Prefix for output files
        taxonomic_level: Taxonomic level for analysis (D, P, C, O, F, G, S)
        group_col: Column name to use for statistical grouping
        skip_kraken: Skip Kraken2 report processing
        skip_bracken: Skip Bracken abundance processing
        skip_downstream: Skip downstream analysis
        min_abundance: Minimum relative abundance threshold for inclusion
        min_prevalence: Minimum prevalence threshold for inclusion
        no_interactive: Disable interactive prompts
        log_file: Path to log file
        
    Returns:
        Tuple of (abundance_file, success_flag)
    """
    # Setup logging
    logger = setup_logger(log_file=log_file)
    log_print("Starting Kraken Analysis Pipeline", level="info")
    start_time = time.time()

    # Process sample metadata
    samples, selected_columns = validate_sample_key(sample_key, no_interactive=no_interactive)

    # Check input directories
    if not skip_kraken and not os.path.isdir(kreport_dir):
        log_print(f"ERROR: Kraken report directory not found: {kreport_dir}", level="error")
        return None, False
    
    if not skip_bracken and not os.path.isdir(bracken_dir):
        log_print(f"ERROR: Bracken directory not found: {bracken_dir}", level="error")
        return None, False

    # Process Kraken/Bracken files
    abundance_file = None
    if not skip_bracken:
        # Look for existing Bracken files
        bracken_files = {}
        for sample in samples:
            # Look for Bracken abundance file patterns
            potential_files = [
                os.path.join(bracken_dir, f"{sample}_{taxonomic_level.lower()}_abundance.txt"),
                os.path.join(bracken_dir, f"{sample}.{taxonomic_level.lower()}.bracken"),
                os.path.join(bracken_dir, f"{sample}_{taxonomic_level.lower()}.bracken"),
                os.path.join(bracken_dir, f"{sample}_bracken_{taxonomic_level.lower()}_abundance.txt")
            ]
            
            for file_path in potential_files:
                if os.path.isfile(file_path):
                    bracken_files[sample] = file_path
                    break
        
        if bracken_files:
            log_print(f"Found {len(bracken_files)} Bracken abundance files at {taxonomic_level} level", level="info")
            # Process Bracken files
            bracken_output_dir = os.path.join(output_dir, "BrackenProcessed")
            abundance_file = process_bracken_abundance(
                bracken_files,
                bracken_output_dir,
                output_prefix,
                taxonomic_level,
                logger=logger
            )
        else:
            log_print(f"No Bracken abundance files found for taxonomic level {taxonomic_level}", level="warning")
    
    # If Bracken processing didn't produce a file, try with Kraken
    if not abundance_file and not skip_kraken:
        log_print("Processing Kraken reports directly...", level="info")
        # TODO: Implement Kraken report processing
        # This would involve parsing .kreport files and converting to a format for analysis
        
    # If we still don't have an abundance file, we can't continue
    if not abundance_file:
        log_print("ERROR: Could not produce an abundance file for analysis", level="error")
        return None, False

    # Downstream analysis
    success = True
    if skip_downstream:
        log_print("Skipping downstream analysis stage", level="info")
    else:
        try:
            # Create analysis output directory
            downstream_out = os.path.join(output_dir, "DownstreamAnalysis")
            os.makedirs(downstream_out, exist_ok=True)
            logger.info(f"Downstream analysis output will be in: {downstream_out}")

            # Read sample metadata
            if not check_file_exists_with_logger(sample_key, "Sample key", logger):
                log_print("Cannot proceed with downstream analysis; missing sample key", level="error")
                success = False
            else:
                sample_key_df = read_and_process_metadata(sample_key, logger)

                # Process abundance data
                abundance_data = read_and_process_taxonomy(
                    abundance_file,
                    sample_key_df,
                    downstream_out,
                    min_abundance=min_abundance,
                    min_prevalence=min_prevalence,
                    logger=logger
                )

                # Run statistical tests
                run_statistical_tests(
                    abundance_data, 
                    downstream_out, 
                    logger, 
                    group_col=group_col
                )
                
                # Generate additional plots
                # TODO: Implement specific taxonomy plots like:
                # - Stacked bar charts by taxonomic level
                # - Heatmap of top taxa
                # - Diversity metrics
                
        except Exception as e:
            logger.error(f"Downstream analysis failed: {e}")
            logger.error(traceback.format_exc())
            success = False

    elapsed = time.time() - start_time
    hh, rr = divmod(elapsed, 3600)
    mm, ss = divmod(rr, 60)
    log_print(f"Pipeline finished in {int(hh)}h {int(mm)}m {int(ss)}s", level="info")

    return abundance_file, success


def process_kraken_files_only(
    sample_key,
    kreport_dir=None,
    bracken_dir=None,
    output_dir="./KrakenProcessed",
    output_prefix="ProcessedFiles",
    taxonomic_level="S",
    skip_kraken=False,
    skip_bracken=False,
    no_interactive=False,
    log_file=None,
):
    """
    Run only the Kraken/Bracken processing stage (no downstream analysis).
    
    Args:
        Same as run_full_pipeline but without downstream parameters
        
    Returns:
        Path to processed abundance file
    """
    # Run full pipeline with skip_downstream=True
    abundance_file, _ = run_full_pipeline(
        sample_key=sample_key,
        kreport_dir=kreport_dir,
        bracken_dir=bracken_dir,
        output_dir=output_dir,
        output_prefix=output_prefix,
        taxonomic_level=taxonomic_level,
        skip_kraken=skip_kraken,
        skip_bracken=skip_bracken,
        skip_downstream=True,
        no_interactive=no_interactive,
        log_file=log_file
    )
    
    return abundance_file


def analyze_existing_kraken_files(
    abundance_file,
    sample_key,
    output_dir,
    group_col="Group",
    min_abundance=0.01,
    min_prevalence=0.1,
    log_file=None,
):
    """
    Run only the downstream analysis on existing Kraken/Bracken abundance files.
    
    Args:
        abundance_file: Path to abundance file
        sample_key: Path to sample metadata CSV
        output_dir: Directory for output files
        group_col: Column name to use for statistical grouping
        min_abundance: Minimum relative abundance threshold for inclusion
        min_prevalence: Minimum prevalence threshold for inclusion
        log_file: Path to log file
        
    Returns:
        Boolean success flag
    """
    # Setup logging
    logger = setup_logger(log_file=log_file)
    log_print("Starting downstream analysis of existing Kraken/Bracken files", level="info")
    start_time = time.time()

    if not abundance_file:
        log_print("Error: abundance_file must be provided", level="error")
        return False

    try:
        # Create analysis output directory
        downstream_out = os.path.join(output_dir, "DownstreamAnalysis")
        os.makedirs(downstream_out, exist_ok=True)
        logger.info(f"Downstream analysis output will be in: {downstream_out}")

        # Read sample metadata
        if not check_file_exists_with_logger(sample_key, "Sample key", logger):
            log_print("Cannot proceed with analysis; missing sample key", level="error")
            return False

        sample_key_df = read_and_process_metadata(sample_key, logger)

        # Process abundance data
        abundance_data = read_and_process_taxonomy(
            abundance_file,
            sample_key_df,
            downstream_out,
            min_abundance=min_abundance,
            min_prevalence=min_prevalence,
            logger=logger
        )

        # Run statistical tests
        run_statistical_tests(abundance_data, downstream_out, logger, group_col=group_col)

        elapsed = time.time() - start_time
        hh, rr = divmod(elapsed, 3600)
        mm, ss = divmod(rr, 60)
        log_print(f"Analysis finished in {int(hh)}h {int(mm)}m {int(ss)}s", level="info")

        return True

    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        logger.error(traceback.format_exc())
        return False


def run_taxonomic_differential_abundance(
    abundance_file,
    sample_key,
    output_dir,
    group_col="Group",
    methods=["aldex2", "ancom", "ancom-bc"],
    min_abundance=0.01,
    min_prevalence=0.1,
    log_file=None,
):
    """
    Run differential abundance tests on taxonomic data.
    
    Args:
        abundance_file: Path to abundance file
        sample_key: Path to sample metadata
        output_dir: Directory for output files
        group_col: Column name for grouping
        methods: List of methods to run
        min_abundance: Minimum relative abundance threshold for inclusion
        min_prevalence: Minimum prevalence threshold for inclusion
        log_file: Path to log file
        
    Returns:
        Dictionary with results from each method
    """
    logger = logging.getLogger("kraken_analysis")

    if not check_file_exists_with_logger(abundance_file, "Abundance file", logger):
        logger.error("Cannot run differential abundance analysis: abundance file not found")
        return None

    if not check_file_exists_with_logger(sample_key, "Sample key", logger):
        logger.error("Cannot run differential abundance analysis: sample key not found")
        return None

    try:
        # Make output directory
        diff_abund_dir = os.path.join(output_dir, "DifferentialAbundance")
        os.makedirs(diff_abund_dir, exist_ok=True)

        # Read data
        abundance_df = pd.read_csv(abundance_file, sep="\t", index_col=0)
        metadata_df = pd.read_csv(sample_key, index_col=None)

        # Get sample ID column (attempt common naming)
        sample_id_col = None
        common_id_names = ["SampleName", "Sample", "SampleID", "Sample_ID", "sample_name", "sample_id"]
        for col in common_id_names:
            if col in metadata_df.columns:
                sample_id_col = col
                break

        if sample_id_col is None:
            logger.error("Could not find sample ID column in metadata")
            return None

        # Set the sample ID as index
        metadata_df = metadata_df.set_index(sample_id_col)

        # Filter low-abundance taxa if needed
        if min_abundance > 0 or min_prevalence > 0:
            before_count = abundance_df.shape[0]
            
            # Calculate mean abundance per taxon
            mean_abundance = abundance_df.mean(axis=1)
            
            # Calculate prevalence (proportion of samples where taxon is present)
            prevalence = (abundance_df > 0).mean(axis=1)
            
            # Apply filters
            abundance_df = abundance_df[(mean_abundance >= min_abundance) & 
                                      (prevalence >= min_prevalence)]
            
            after_count = abundance_df.shape[0]
            logger.info(f"Filtered from {before_count} to {after_count} taxa based on abundance/prevalence thresholds")

        # Run differential abundance analysis
        results = run_differential_abundance_analysis(
            abundance_df,
            metadata_df,
            diff_abund_dir,
            group_col=group_col,
            methods=methods,
            logger=logger
        )

        return results

    except Exception as e:
        logger.error(f"Error in differential abundance analysis: {str(e)}")
        logger.error(traceback.format_exc())
        return None
