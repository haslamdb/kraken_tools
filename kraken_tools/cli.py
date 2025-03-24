# kraken_tools/cli.py
import os
import sys
import argparse
import time
import logging
import traceback

from kraken_tools.logger import setup_logger, log_print
from kraken_tools.utils.sample_utils import validate_sample_key, check_input_files_exist
from kraken_tools.utils.file_utils import check_file_exists_with_logger
from kraken_tools.utils.resource_utils import limit_memory_usage
from kraken_tools.preprocessing.kneaddata import check_kneaddata_installation
from kraken_tools.preprocessing.kraken_run import check_kraken_installation
from kraken_tools.preprocessing.bracken_run import check_bracken_installation
from kraken_tools.preprocessing.pipeline import run_preprocessing_pipeline, run_preprocessing_pipeline_parallel
from kraken_tools.main import (
    run_full_pipeline,
    process_kraken_files_only,
    analyze_existing_kraken_files,
    run_taxonomic_differential_abundance
)
from kraken_tools.analysis.permanova import run_permanova_analysis
from kraken_tools.analysis.feature_selection import run_feature_selection
from kraken_tools.analysis.rf_shap import run_rf_shap_analysis


def setup_common_args(parser):
    """Add common arguments to a parser."""
    parser.add_argument("--log-file", default=None, help="Path to log file")
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default=INFO)",
    )
    parser.add_argument(
        "--max-memory", type=int, default=None, help="Maximum memory usage in MB (default: unlimited)"
    )
    parser.add_argument(
        "--no-interactive", action="store_true", help="Non-interactive mode for sample key column selection"
    )
    return parser


def setup_input_output_args(parser):
    """Add input/output arguments to a parser."""
    parser.add_argument("--output-dir", required=True, help="Directory for output files")
    parser.add_argument(
        "--output-prefix",
        default="ProcessedFiles",
        help="Prefix for intermediate output files",
    )
    return parser


def setup_preprocessing_args(parser):
    """Add preprocessing arguments to a parser."""
    parser.add_argument("--input-fastq", nargs="+", help="Input FASTQ file(s)")
    parser.add_argument(
        "--paired", action="store_true", help="Input files are paired-end reads"
    )
    parser.add_argument(
        "--kneaddata-dbs", nargs="+", help="Path(s) to KneadData reference database(s)"
    )
    parser.add_argument(
        "--threads", type=int, default=1, help="Number of threads to use"
    )
    return parser


def setup_kraken_args(parser):
    """Add Kraken/Bracken arguments to a parser."""
    parser.add_argument(
        "--kraken-db", help="Path to Kraken2 database"
    )
    parser.add_argument(
        "--bracken-db", help="Path to Bracken database file"
    )
    parser.add_argument(
        "--taxonomic-level", 
        default="S", 
        choices=["D", "P", "C", "O", "F", "G", "S"],
        help="Taxonomic level for analysis (D=domain, P=phylum, C=class, O=order, F=family, G=genus, S=species)"
    )
    parser.add_argument(
        "--threshold", 
        type=float, 
        default=10.0,
        help="Threshold for Bracken (minimum number of reads required)"
    )
    return parser


def setup_analysis_args(parser):
    """Add analysis arguments to a parser."""
    parser.add_argument("--sample-key", required=True, help="CSV file with sample metadata")
    parser.add_argument(
        "--group-col", default="Group", help="Column name for grouping in stats"
    )
    parser.add_argument(
        "--min-abundance", type=float, default=0.01, 
        help="Minimum relative abundance threshold (default: 0.01 = 1%%)"
    )
    parser.add_argument(
        "--min-prevalence", type=float, default=0.1, 
        help="Minimum prevalence threshold (default: 0.1 = 10%%)"
    )
    return parser


def setup_diff_abundance_args(parser):
    """Add differential abundance arguments to a parser."""
    parser.add_argument(
        "--methods",
        default="aldex2,ancom,ancom-bc",
        help="Comma-separated list of methods (default: aldex2,ancom,ancom-bc)",
    )
    return parser


def setup_parallel_args(parser):
    """Add parallel processing arguments to a parser."""
    parser.add_argument(
        "--threads-per-sample", type=int, default=1, help="Number of threads per sample"
    )
    parser.add_argument(
        "--max-parallel", type=int, default=None, help="Maximum samples to process in parallel"
    )
    parser.add_argument(
        "--use-parallel", action="store_true", help="Use parallel processing"
    )
    return parser


def main():
    """Main entry point for the kraken_tools CLI."""
    
    # Create the top-level parser
    parser = argparse.ArgumentParser(
        description="Kraken Tools: Process and analyze Kraken2/Bracken taxonomic output",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example commands:
  # Run full pipeline (raw reads to analysis)
  kraken-tools full-pipeline --input-fastq reads_1.fastq reads_2.fastq --paired --sample-key metadata.csv --output-dir results/
  
  # Process existing Kraken/Bracken files
  kraken-tools process --kreport-dir kraken_reports/ --bracken-dir bracken_files/ --sample-key metadata.csv --output-dir results/
  
  # Run only differential abundance analysis
  kraken-tools diff-abundance --abundance-file abundance.tsv --sample-key metadata.csv --output-dir results/diff_abundance/
  
See documentation for more examples and detailed parameter descriptions.
        """
    )
    
    # Create subparsers for different commands
    subparsers = parser.add_subparsers(dest="command", help="Command to run")
    
    # 1. Full pipeline command
    full_pipeline_parser = subparsers.add_parser("full-pipeline", help="Run the complete pipeline from raw reads to analysis")
    full_pipeline_parser = setup_common_args(full_pipeline_parser)
    full_pipeline_parser = setup_input_output_args(full_pipeline_parser)
    full_pipeline_parser = setup_preprocessing_args(full_pipeline_parser)
    full_pipeline_parser = setup_kraken_args(full_pipeline_parser)
    full_pipeline_parser = setup_analysis_args(full_pipeline_parser)
    full_pipeline_parser = setup_parallel_args(full_pipeline_parser)
    full_pipeline_parser.add_argument(
        "--skip-preprocessing", action="store_true", help="Skip preprocessing (KneadData) step"
    )
    full_pipeline_parser.add_argument(
        "--skip-classification", action="store_true", help="Skip classification (Kraken/Bracken) step"
    )
    full_pipeline_parser.add_argument(
        "--skip-downstream", action="store_true", help="Skip downstream analysis"
    )
    full_pipeline_parser.add_argument(
        "--run-diff-abundance", action="store_true", help="Run differential abundance analysis"
    )
    
    # 2. Preprocess command
    preprocess_parser = subparsers.add_parser("preprocess", help="Run preprocessing only (KneadData)")
    preprocess_parser = setup_common_args(preprocess_parser)
    preprocess_parser = setup_input_output_args(preprocess_parser)
    preprocess_parser = setup_preprocessing_args(preprocess_parser)
    preprocess_parser = setup_parallel_args(preprocess_parser)
    
    # 3. Classify command
    classify_parser = subparsers.add_parser("classify", help="Run taxonomic classification only (Kraken2 + Bracken)")
    classify_parser = setup_common_args(classify_parser)
    classify_parser = setup_input_output_args(classify_parser)
    classify_parser = setup_preprocessing_args(classify_parser)
    classify_parser = setup_kraken_args(classify_parser)
    classify_parser = setup_parallel_args(classify_parser)
    classify_parser.add_argument(
        "--skip-kraken", action="store_true", help="Skip Kraken2 step (use existing reports)"
    )
    classify_parser.add_argument(
        "--skip-bracken", action="store_true", help="Skip Bracken step"
    )
    
    # 4. Process command
    process_parser = subparsers.add_parser("process", help="Process existing Kraken/Bracken files")
    process_parser = setup_common_args(process_parser)
    process_parser = setup_input_output_args(process_parser)
    process_parser = setup_analysis_args(process_parser)
    process_parser.add_argument(
        "--kreport-dir", help="Directory containing Kraken2 report files"
    )
    process_parser.add_argument(
        "--bracken-dir", help="Directory containing Bracken abundance files"
    )
    process_parser.add_argument(
        "--skip-kraken", action="store_true", help="Skip Kraken2 report processing"
    )
    process_parser.add_argument(
        "--skip-bracken", action="store_true", help="Skip Bracken abundance processing"
    )
    
    # 5. Analyze command
    analyze_parser = subparsers.add_parser("analyze", help="Run downstream analysis on processed abundance data")
    analyze_parser = setup_common_args(analyze_parser)
    analyze_parser = setup_input_output_args(analyze_parser)
    analyze_parser = setup_analysis_args(analyze_parser)
    analyze_parser.add_argument(
        "--abundance-file", required=True, help="Path to abundance file"
    )
    
    # 6. Diff-abundance command
    diff_abundance_parser = subparsers.add_parser("diff-abundance", help="Run differential abundance analysis")
    diff_abundance_parser = setup_common_args(diff_abundance_parser)
    diff_abundance_parser = setup_input_output_args(diff_abundance_parser)
    diff_abundance_parser = setup_analysis_args(diff_abundance_parser)
    diff_abundance_parser = setup_diff_abundance_args(diff_abundance_parser)
    diff_abundance_parser.add_argument(
        "--abundance-file", required=True, help="Path to abundance file"
    )
    
    # 7. GLMM command
    glmm_parser = subparsers.add_parser("glmm", help="Run GLMM analysis")
    glmm_parser = setup_common_args(glmm_parser)
    glmm_parser = setup_input_output_args(glmm_parser)
    glmm_parser = setup_analysis_args(glmm_parser)
    glmm_parser.add_argument(
        "--abundance-file", required=True, help="Path to abundance file"
    )
    glmm_parser.add_argument(
        "--formula", required=True, help="R-style formula for GLMM (e.g., 'Count ~ Group + (1|Subject)')"
    )
    glmm_parser.add_argument(
        "--model", default="negbin", choices=["poisson", "negbin"], 
        help="Model family for GLMM"
    )
    
    # 8. List-files command (utility)
    list_files_parser = subparsers.add_parser("list-files", help="List files in specified directories")
    list_files_parser = setup_common_args(list_files_parser)
    list_files_parser.add_argument("--kreport-dir", help="Directory containing Kraken2 report files")
    list_files_parser.add_argument("--bracken-dir", help="Directory containing Bracken abundance files")
    
    # 9. PERMANOVA command
    permanova_parser = subparsers.add_parser("permanova", help="Run PERMANOVA analysis")
    permanova_parser = setup_common_args(permanova_parser)
    permanova_parser = setup_input_output_args(permanova_parser)
    permanova_parser = setup_analysis_args(permanova_parser)
    permanova_parser.add_argument(
        "--abundance-file", required=True, help="Path to abundance file"
    )
    permanova_parser.add_argument(
        "--categorical-vars", help="Comma-separated list of categorical variables to test"
    )
    permanova_parser.add_argument(
        "--distance-metric", default="bray", choices=["bray", "jaccard", "euclidean"],
        help="Distance metric to use"
    )
    permanova_parser.add_argument(
        "--transform", default="clr", choices=["clr", "hellinger", "log", "none"],
        help="Transformation to apply to abundance data"
    )
    permanova_parser.add_argument(
        "--permutations", type=int, default=999, help="Number of permutations"
    )
    permanova_parser.add_argument(
        "--min-group-size", type=int, default=3, help="Minimum number of samples per group"
    )
    permanova_parser.add_argument(
        "--make-pcoa", action="store_true", default=True, help="Generate PCoA plots"
    )

    # 10. Feature Selection command
    feature_selection_parser = subparsers.add_parser("feature-selection", help="Run Random Forest feature selection")
    feature_selection_parser = setup_common_args(feature_selection_parser)
    feature_selection_parser = setup_input_output_args(feature_selection_parser)
    feature_selection_parser = setup_analysis_args(feature_selection_parser)
    feature_selection_parser.add_argument(
        "--abundance-file", required=True, help="Path to abundance file"
    )
    feature_selection_parser.add_argument(
        "--predictors", help="Comma-separated list of predictor variables"
    )
    feature_selection_parser.add_argument(
        "--n-estimators", type=int, default=100, help="Number of trees in Random Forest"
    )
    feature_selection_parser.add_argument(
        "--distance-metric", default="bray", choices=["bray", "jaccard", "euclidean"],
        help="Distance metric to use"
    )
    feature_selection_parser.add_argument(
        "--transform", default="clr", choices=["clr", "hellinger", "log", "none"],
        help="Transformation to apply to abundance data"
    )
    feature_selection_parser.add_argument(
        "--test-size", type=float, default=0.2, help="Proportion of data for testing"
    )
    feature_selection_parser.add_argument(
        "--random-state", type=int, default=42, help="Random seed for reproducibility"
    )

    # 11. RF-SHAP command
    rf_shap_parser = subparsers.add_parser("rf-shap", help="Run Random Forest with SHAP analysis")
    rf_shap_parser = setup_common_args(rf_shap_parser)
    rf_shap_parser = setup_input_output_args(rf_shap_parser)
    rf_shap_parser = setup_analysis_args(rf_shap_parser)
    rf_shap_parser.add_argument(
        "--abundance-file", required=True, help="Path to abundance file"
    )
    rf_shap_parser.add_argument(
        "--target-taxa", help="Comma-separated list of taxa to analyze"
    )
    rf_shap_parser.add_argument(
        "--predictors", help="Comma-separated list of predictor variables"
    )
    rf_shap_parser.add_argument(
        "--random-effects", help="Comma-separated list of random effect variables"
    )
    rf_shap_parser.add_argument(
        "--transform", default="clr", choices=["clr", "hellinger", "log", "none"],
        help="Transformation to apply to abundance data"
    )
    rf_shap_parser.add_argument(
        "--n-estimators", type=int, default=100, help="Number of trees in Random Forest"
    )
    rf_shap_parser.add_argument(
        "--test-size", type=float, default=0.2, help="Proportion of data for testing"
    )
    rf_shap_parser.add_argument(
        "--top-n", type=int, default=10, help="Number of top taxa to analyze if target-taxa not specified"
    )
    rf_shap_parser.add_argument(
        "--mixed-model", default="lmer", choices=["lmer", "glmm"],
        help="Type of mixed model to use"
    )

    # Parse arguments
    args = parser.parse_args()
    
    # If no command specified, show help
    if args.command is None:
        parser.print_help()
        sys.exit(0)
    
    # Setup logging
    logger = setup_logger(log_file=args.log_file, log_level=getattr(logging, args.log_level))
    log_print(f"Starting Kraken Tools - {args.command}", level="info")
    
    start_time = time.time()
    
    # Apply memory limit if specified
    if hasattr(args, "max_memory") and args.max_memory:
        success = limit_memory_usage(args.max_memory)
        if success:
            log_print(f"Set memory limit to {args.max_memory} MB", level="info")
        else:
            log_print("Failed to set memory limit, proceeding with unlimited memory", level="warning")
    
    # Handle list-files command
    if args.command == "list-files":
        if args.kreport_dir:
            log_print(f"Files in Kraken report dir: {args.kreport_dir}", level="info")
            if os.path.isdir(args.kreport_dir):
                for f in sorted(os.listdir(args.kreport_dir)):
                    log_print("  " + f, level="info")
            else:
                log_print("  Kraken report dir not found.", level="warning")
        
        if args.bracken_dir:
            log_print(f"Files in Bracken dir: {args.bracken_dir}", level="info")
            if os.path.isdir(args.bracken_dir):
                for f in sorted(os.listdir(args.bracken_dir)):
                    log_print("  " + f, level="info")
            else:
                log_print("  Bracken dir not found.", level="warning")
        
        sys.exit(0)
    
    # Handle full-pipeline command
    if args.command == "full-pipeline":
        # Validate sample key
        samples, selected_columns = validate_sample_key(args.sample_key, no_interactive=args.no_interactive)
        
        # Setup output directories
        os.makedirs(args.output_dir, exist_ok=True)
        preproc_dir = os.path.join(args.output_dir, "PreprocessedData")
        taxonomy_dir = os.path.join(args.output_dir, "TaxonomyData")
        processed_dir = os.path.join(args.output_dir, "ProcessedData")
        analysis_dir = os.path.join(args.output_dir, "AnalysisResults")
        
        # Preprocessing step
        if not args.skip_preprocessing:
            if not args.input_fastq:
                log_print("ERROR: --input-fastq is required for preprocessing", level="error")
                sys.exit(1)
            
            # Check KneadData installation
            kneaddata_ok, kneaddata_msg = check_kneaddata_installation()
            if not kneaddata_ok:
                log_print(f"ERROR: KneadData not properly installed: {kneaddata_msg}", level="error")
                sys.exit(1)
            
            # Create output directory
            os.makedirs(preproc_dir, exist_ok=True)
            
            # Run preprocessing
            log_print("Starting preprocessing (KneadData) step...", level="info")
            
            if args.use_parallel:
                preprocessing_results = run_preprocessing_pipeline_parallel(
                    input_files=args.input_fastq,
                    output_dir=preproc_dir,
                    threads_per_sample=args.threads_per_sample,
                    max_parallel=args.max_parallel,
                    kneaddata_dbs=args.kneaddata_dbs,
                    paired=args.paired,
                    logger=logger
                )
            else:
                preprocessing_results = run_preprocessing_pipeline(
                    input_files=args.input_fastq,
                    output_dir=preproc_dir,
                    threads=args.threads,
                    kneaddata_dbs=args.kneaddata_dbs,
                    paired=args.paired,
                    logger=logger
                )
            
            if not preprocessing_results:
                log_print("ERROR: Preprocessing failed", level="error")
                sys.exit(1)
            
            # Get clean reads for next step
            kneaddata_files = preprocessing_results.get('kneaddata_files', [])
            log_print(f"Preprocessing completed with {len(kneaddata_files)} clean read files", level="info")
            
            # Use these as input for classification
            input_for_classification = kneaddata_files
        else:
            input_for_classification = args.input_fastq
        
        # Classification step
        kreport_dir = None
        bracken_dir = None
        
        if not args.skip_classification:
            if not input_for_classification:
                log_print("ERROR: No input files for classification", level="error")
                sys.exit(1)
            
            # Check installations
            kraken_ok, kraken_msg = check_kraken_installation()
            if not kraken_ok:
                log_print(f"ERROR: Kraken2 not properly installed: {kraken_msg}", level="error")
                sys.exit(1)
            
            bracken_ok, bracken_msg = check_bracken_installation()
            if not bracken_ok:
                log_print(f"ERROR: Bracken not properly installed: {bracken_msg}", level="error")
                sys.exit(1)
            
            # Create output directories
            os.makedirs(taxonomy_dir, exist_ok=True)
            kreport_dir = os.path.join(taxonomy_dir, "kraken_reports")
            bracken_dir = os.path.join(taxonomy_dir, "bracken_output")
            os.makedirs(kreport_dir, exist_ok=True)
            os.makedirs(bracken_dir, exist_ok=True)
            
            # Run Kraken2
            if not args.skip_kraken:
                log_print("Starting Kraken2 classification...", level="info")
                
                from kraken_tools.preprocessing.kraken_run import run_kraken, run_kraken_parallel
                
                if args.use_parallel:
                    kraken_results = run_kraken_parallel(
                        input_files=input_for_classification,
                        output_dir=kreport_dir,
                        threads=args.threads_per_sample,
                        max_parallel=args.max_parallel,
                        kraken_db=args.kraken_db,
                        paired=False,  # KneadData outputs are single-end
                        logger=logger
                    )
                else:
                    kraken_results = run_kraken(
                        input_files=input_for_classification,
                        output_dir=kreport_dir,
                        threads=args.threads,
                        kraken_db=args.kraken_db,
                        paired=False,  # KneadData outputs are single-end
                        logger=logger
                    )
                
                if not kraken_results:
                    log_print("ERROR: Kraken2 classification failed", level="error")
                    sys.exit(1)
                
                kreport_files = {sample_id: results['report'] for sample_id, results in kraken_results.items()}
                log_print(f"Kraken2 completed with {len(kreport_files)} reports", level="info")
            
            # Run Bracken
            if not args.skip_bracken:
                log_print("Starting Bracken abundance estimation...", level="info")
                
                from kraken_tools.preprocessing.bracken_run import run_bracken, run_bracken_parallel
                
                if args.use_parallel:
                    bracken_results = run_bracken_parallel(
                        kreport_files=kreport_files if 'kreport_files' in locals() else None,
                        output_dir=bracken_dir,
                        threads=args.threads_per_sample,
                        max_parallel=args.max_parallel,
                        bracken_db=args.bracken_db,
                        taxonomic_level=args.taxonomic_level,
                        threshold=args.threshold,
                        logger=logger
                    )
                else:
                    bracken_results = run_bracken(
                        kreport_files=kreport_files if 'kreport_files' in locals() else None,
                        output_dir=bracken_dir,
                        threads=args.threads,
                        bracken_db=args.bracken_db,
                        taxonomic_level=args.taxonomic_level,
                        threshold=args.threshold,
                        logger=logger
                    )
                
                if not bracken_results:
                    log_print("ERROR: Bracken abundance estimation failed", level="error")
                    sys.exit(1)
                
                log_print(f"Bracken completed for {len(bracken_results)} samples", level="info")
    
    # Handle process command
    elif args.command == "process":
        if not args.kreport_dir and not args.bracken_dir:
            log_print("ERROR: At least one of --kreport-dir or --bracken-dir is required", level="error")
            sys.exit(1)
        
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Process Kraken/Bracken files
        log_print("Processing taxonomic files...", level="info")
        abundance_file = process_kraken_files_only(
            sample_key=args.sample_key,
            kreport_dir=args.kreport_dir,
            bracken_dir=args.bracken_dir,
            output_dir=args.output_dir,
            output_prefix=args.output_prefix,
            taxonomic_level=args.taxonomic_level if hasattr(args, 'taxonomic_level') else "S",
            skip_kraken=args.skip_kraken,
            skip_bracken=args.skip_bracken,
            no_interactive=args.no_interactive,
            log_file=args.log_file
        )
        
        if not abundance_file:
            log_print("ERROR: Processing taxonomic files failed", level="error")
            sys.exit(1)
        
        log_print(f"File processing completed. Abundance file: {abundance_file}", level="info")
    
    # Handle analyze command
    elif args.command == "analyze":
        if not args.abundance_file:
            log_print("ERROR: --abundance-file is required for analysis", level="error")
            sys.exit(1)
        
        if not os.path.isfile(args.abundance_file):
            log_print(f"ERROR: Abundance file not found: {args.abundance_file}", level="error")
            sys.exit(1)
        
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Run downstream analysis
        log_print("Running downstream analysis...", level="info")
        success = analyze_existing_kraken_files(
            abundance_file=args.abundance_file,
            sample_key=args.sample_key,
            output_dir=args.output_dir,
            group_col=args.group_col,
            min_abundance=args.min_abundance,
            min_prevalence=args.min_prevalence,
            log_file=args.log_file
        )
        
        if not success:
            log_print("ERROR: Downstream analysis failed", level="error")
            sys.exit(1)
        
        log_print("Downstream analysis completed successfully", level="info")
    
    # Handle diff-abundance command
    elif args.command == "diff-abundance":
        if not args.abundance_file:
            log_print("ERROR: --abundance-file is required for differential abundance analysis", level="error")
            sys.exit(1)
        
        if not os.path.isfile(args.abundance_file):
            log_print(f"ERROR: Abundance file not found: {args.abundance_file}", level="error")
            sys.exit(1)
        
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Run differential abundance analysis
        log_print("Running differential abundance analysis...", level="info")
        methods = args.methods.split(",")
        
        results = run_taxonomic_differential_abundance(
            abundance_file=args.abundance_file,
            sample_key=args.sample_key,
            output_dir=args.output_dir,
            group_col=args.group_col,
            methods=methods,
            min_abundance=args.min_abundance,
            min_prevalence=args.min_prevalence,
            log_file=args.log_file
        )
        
        if not results:
            log_print("ERROR: Differential abundance analysis failed", level="error")
            sys.exit(1)
        
        log_print(f"Differential abundance analysis completed with {len(results)} methods", level="info")
    
    # Handle glmm command
    elif args.command == "glmm":
        if not args.abundance_file:
            log_print("ERROR: --abundance-file is required for GLMM analysis", level="error")
            sys.exit(1)
        
        if not os.path.isfile(args.abundance_file):
            log_print(f"ERROR: Abundance file not found: {args.abundance_file}", level="error")
            sys.exit(1)
        
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Run GLMM analysis
        log_print("Running GLMM analysis...", level="info")
        
        try:
            # Import GLMM module
            from kraken_tools.analysis.glmm_analysis import run_glmm_analysis
            
            success = run_glmm_analysis(
                abundance_file=args.abundance_file,
                sample_key=args.sample_key,
                output_dir=args.output_dir,
                formula=args.formula,
                model=args.model,
                group_col=args.group_col,
                min_abundance=args.min_abundance,
                min_prevalence=args.min_prevalence,
                logger=logger
            )
            
            if not success:
                log_print("ERROR: GLMM analysis failed", level="error")
                sys.exit(1)
            
            log_print("GLMM analysis completed successfully", level="info")
        except ImportError:
            log_print("ERROR: GLMM analysis module not found", level="error")
            sys.exit(1)
            
    # Handle permanova command
    elif args.command == "permanova":
        if not args.abundance_file:
            log_print("ERROR: --abundance-file is required for PERMANOVA analysis", level="error")
            sys.exit(1)
        
        if not os.path.isfile(args.abundance_file):
            log_print(f"ERROR: Abundance file not found: {args.abundance_file}", level="error")
            sys.exit(1)
        
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Run PERMANOVA analysis
        log_print("Running PERMANOVA analysis...", level="info")
        results = run_permanova_analysis(
            abundance_file=args.abundance_file,
            metadata_file=args.sample_key,
            output_dir=args.output_dir,
            categorical_vars=args.categorical_vars,
            group_col=args.group_col,
            distance_metric=args.distance_metric,
            transform=args.transform,
            permutations=args.permutations,
            min_group_size=args.min_group_size,
            make_pcoa=args.make_pcoa,
            log_file=args.log_file
        )
        
        if not results:
            log_print("ERROR: PERMANOVA analysis failed", level="error")
            sys.exit(1)
        
        sig_count = sum(1 for r in results.values() if r['p_value'] < 0.05)
        log_print(f"PERMANOVA analysis completed with {sig_count} significant variables", level="info")

    # Handle feature-selection command
    elif args.command == "feature-selection":
        if not args.abundance_file:
            log_print("ERROR: --abundance-file is required for feature selection", level="error")
            sys.exit(1)
        
        if not os.path.isfile(args.abundance_file):
            log_print(f"ERROR: Abundance file not found: {args.abundance_file}", level="error")
            sys.exit(1)
        
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Run feature selection
        log_print("Running Random Forest feature selection...", level="info")
        results = run_feature_selection(
            abundance_file=args.abundance_file,
            metadata_file=args.sample_key,
            output_dir=args.output_dir,
            predictors=args.predictors,
            n_estimators=args.n_estimators,
            distance_metric=args.distance_metric,
            transform=args.transform,
            test_size=args.test_size,
            random_state=args.random_state,
            log_file=args.log_file
        )
        
        if results is None:
            log_print("ERROR: Feature selection failed", level="error")
            sys.exit(1)
        
        log_print(f"Feature selection completed with {len(results)} ranked features", level="info")

    # Handle rf-shap command
    elif args.command == "rf-shap":
        if not args.abundance_file:
            log_print("ERROR: --abundance-file is required for RF-SHAP analysis", level="error")
            sys.exit(1)
        
        if not os.path.isfile(args.abundance_file):
            log_print(f"ERROR: Abundance file not found: {args.abundance_file}", level="error")
            sys.exit(1)
        
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Run RF-SHAP analysis
        log_print("Running Random Forest with SHAP analysis...", level="info")
        results = run_rf_shap_analysis(
            abundance_file=args.abundance_file,
            metadata_file=args.sample_key,
            output_dir=args.output_dir,
            target_taxa=args.target_taxa,
            predictors=args.predictors,
            random_effects=args.random_effects,
            transform=args.transform,
            n_estimators=args.n_estimators,
            test_size=args.test_size,
            top_n=args.top_n,
            mixed_model=args.mixed_model,
            log_file=args.log_file
        )
        
        if not results:
            log_print("ERROR: RF-SHAP analysis failed", level="error")
            sys.exit(1)
        
        log_print(f"RF-SHAP analysis completed for {len(results)} taxa", level="info")

    # Calculate and report elapsed time
    elapsed = time.time() - start_time
    hh, rr = divmod(elapsed, 3600)
    mm, ss = divmod(rr, 60)
    log_print(f"Command '{args.command}' finished in {int(hh)}h {int(mm)}m {int(ss)}s", level="info")


if __name__ == "__main__":
    main()
        
        # Process files step
        abundance_file = None
        
        if not args.skip_kraken and not kreport_dir:
            kreport_dir = os.path.join(taxonomy_dir, "kraken_reports")
            
        if not args.skip_bracken and not bracken_dir:
            bracken_dir = os.path.join(taxonomy_dir, "bracken_output")
        
        log_print("Processing taxonomic files...", level="info")
        abundance_file, success = run_full_pipeline(
            sample_key=args.sample_key,
            kreport_dir=kreport_dir,
            bracken_dir=bracken_dir,
            output_dir=processed_dir,
            output_prefix=args.output_prefix,
            taxonomic_level=args.taxonomic_level,
            group_col=args.group_col,
            skip_kraken=args.skip_kraken,
            skip_bracken=args.skip_bracken,
            skip_downstream=args.skip_downstream,
            min_abundance=args.min_abundance,
            min_prevalence=args.min_prevalence,
            no_interactive=args.no_interactive,
            log_file=args.log_file
        )
        
        if not success:
            log_print("ERROR: Processing taxonomic files failed", level="error")
            sys.exit(1)
        
        log_print(f"File processing completed. Abundance file: {abundance_file}", level="info")
        
        # Run differential abundance analysis if requested
        if args.run_diff_abundance and abundance_file:
            log_print("Running differential abundance analysis...", level="info")
            diff_abundance_dir = os.path.join(args.output_dir, "DifferentialAbundance")
            os.makedirs(diff_abundance_dir, exist_ok=True)
            
            results = run_taxonomic_differential_abundance(
                abundance_file=abundance_file,
                sample_key=args.sample_key,
                output_dir=diff_abundance_dir,
                group_col=args.group_col,
                methods="aldex2,ancom,ancom-bc".split(","),
                min_abundance=args.min_abundance,
                min_prevalence=args.min_prevalence,
                log_file=args.log_file
            )
            
            if not results:
                log_print("WARNING: Differential abundance analysis did not produce results", level="warning")
    
    # Handle preprocess command
    elif args.command == "preprocess":
        if not args.input_fastq:
            log_print("ERROR: --input-fastq is required for preprocessing", level="error")
            sys.exit(1)
        
        # Check KneadData installation
        kneaddata_ok, kneaddata_msg = check_kneaddata_installation()
        if not kneaddata_ok:
            log_print(f"ERROR: KneadData not properly installed: {kneaddata_msg}", level="error")
            sys.exit(1)
        
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Run preprocessing
        log_print("Starting preprocessing (KneadData)...", level="info")
        
        if args.use_parallel:
            preprocessing_results = run_preprocessing_pipeline_parallel(
                input_files=args.input_fastq,
                output_dir=args.output_dir,
                threads_per_sample=args.threads_per_sample,
                max_parallel=args.max_parallel,
                kneaddata_dbs=args.kneaddata_dbs,
                paired=args.paired,
                logger=logger
            )
        else:
            preprocessing_results = run_preprocessing_pipeline(
                input_files=args.input_fastq,
                output_dir=args.output_dir,
                threads=args.threads,
                kneaddata_dbs=args.kneaddata_dbs,
                paired=args.paired,
                logger=logger
            )
        
        if not preprocessing_results:
            log_print("ERROR: Preprocessing failed", level="error")
            sys.exit(1)
        
        kneaddata_files = preprocessing_results.get('kneaddata_files', [])
        log_print(f"Preprocessing completed with {len(kneaddata_files)} clean read files", level="info")
    
    # Handle classify command
    elif args.command == "classify":
        if not args.input_fastq:
            log_print("ERROR: --input-fastq is required for classification", level="error")
            sys.exit(1)
        
        # Check installations
        if not args.skip_kraken:
            kraken_ok, kraken_msg = check_kraken_installation()
            if not kraken_ok:
                log_print(f"ERROR: Kraken2 not properly installed: {kraken_msg}", level="error")
                sys.exit(1)
        
        if not args.skip_bracken:
            bracken_ok, bracken_msg = check_bracken_installation()
            if not bracken_ok:
                log_print(f"ERROR: Bracken not properly installed: {bracken_msg}", level="error")
                sys.exit(1)
        
        # Create output directories
        os.makedirs(args.output_dir, exist_ok=True)
        kreport_dir = os.path.join(args.output_dir, "kraken_reports")
        bracken_dir = os.path.join(args.output_dir, "bracken_output")
        os.makedirs(kreport_dir, exist_ok=True)
        os.makedirs(bracken_dir, exist_ok=True)
        
        # Run Kraken2
        kreport_files = {}
        if not args.skip_kraken:
            log_print("Starting Kraken2 classification...", level="info")
            
            from kraken_tools.preprocessing.kraken_run import run_kraken, run_kraken_parallel
            
            if args.use_parallel:
                kraken_results = run_kraken_parallel(
                    input_files=args.input_fastq,
                    output_dir=kreport_dir,
                    threads=args.threads_per_sample,
                    max_parallel=args.max_parallel,
                    kraken_db=args.kraken_db,
                    paired=args.paired,
                    logger=logger
                )
            else:
                kraken_results = run_kraken(
                    input_files=args.input_fastq,
                    output_dir=kreport_dir,
                    threads=args.threads,
                    kraken_db=args.kraken_db,
                    paired=args.paired,
                    logger=logger
                )
            
            if not kraken_results:
                log_print("ERROR: Kraken2 classification failed", level="error")
                sys.exit(1)
            
            kreport_files = {sample_id: results['report'] for sample_id, results in kraken_results.items()}
            log_print(f"Kraken2 completed with {len(kreport_files)} reports", level="info")
        
        # Run Bracken
        if not args.skip_bracken:
            log_print("Starting Bracken abundance estimation...", level="info")
            
            from kraken_tools.preprocessing.bracken_run import run_bracken, run_bracken_parallel
            
            if args.use_parallel:
                bracken_results = run_bracken_parallel(
                    kreport_files=kreport_files if kreport_files else None,
                    output_dir=bracken_dir,
                    threads=args.threads_per_sample,
                    max_parallel=args.max_parallel,
                    bracken_db=args.bracken_db,
                    taxonomic_level=args.taxonomic_level,
                    threshold=args.threshold,
                    logger=logger
                )
            else:
                bracken_results = run_bracken(
                    kreport_files=kreport_files if kreport_files else None,
                    output_dir=bracken_dir,
                    threads=args.threads,
                    bracken_db=args.bracken_db,
                    taxonomic_level=args.taxonomic_level,
                    threshold=args.threshold,
                    logger=logger
                )
            
            if not bracken_results:
                log_print("ERROR: Bracken abundance estimation failed", level="error")
                sys.exit(1)
            
            log_print(f"Bracken completed for {len(bracken_results)} samples", level="info")