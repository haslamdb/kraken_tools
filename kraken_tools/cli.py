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


def main():
    """Main entry point for the kraken_tools CLI."""
    parser = argparse.ArgumentParser(description="Kraken Tools: Process and analyze Kraken2/Bracken taxonomic output")

    # --- 1. Global Options ---
    global_group = parser.add_argument_group("Global Options")
    global_group.add_argument("--log-file", default=None, help="Path to combined log file")
    global_group.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default=INFO)",
    )
    global_group.add_argument(
        "--max-memory", type=int, default=None, help="Maximum memory usage in MB (default: unlimited)"
    )
    global_group.add_argument(
        "--list-files", action="store_true", help="Just list input files in specified directories, then exit"
    )
    global_group.add_argument(
        "--no-interactive", action="store_true", help="Non-interactive mode for sample key column selection"
    )

    # --- 2. Input/Output Options ---
    io_group = parser.add_argument_group("Input/Output Options")
    io_group.add_argument("--sample-key", required=True, help="CSV file with columns for sample names and metadata")
    io_group.add_argument(
        "--output-dir",
        default="./KrakenOutput",
        help="Directory where processed files and downstream results will go",
    )
    io_group.add_argument(
        "--output-prefix",
        default="ProcessedFiles",
        help="Prefix for intermediate output directories/files",
    )

    # --- 3. Preprocessing Options ---
    preprocessing_group = parser.add_argument_group("Preprocessing Options")
    preprocessing_group.add_argument(
        "--run-preprocessing", action="store_true", help="Run preprocessing (KneadData and Kraken2/Bracken) on raw sequence files"
    )
    preprocessing_group.add_argument("--input-fastq", nargs="+", help="Input FASTQ file(s) for preprocessing")
    preprocessing_group.add_argument(
        "--paired", action="store_true", help="Input files are paired-end reads (default: False)"
    )
    preprocessing_group.add_argument(
        "--kneaddata-dbs", nargs="+", help="Path(s) to KneadData reference database(s). Can specify multiple databases."
    )
    preprocessing_group.add_argument(
        "--kraken-db", help="Path to Kraken2 database for taxonomic classification"
    )
    preprocessing_group.add_argument(
        "--bracken-db", help="Path to Bracken database file (e.g., database150mers.kmer_distrib)"
    )
    preprocessing_group.add_argument(
        "--threads", type=int, default=1, help="Number of threads to use for preprocessing"
    )
    preprocessing_group.add_argument(
        "--kneaddata-output-dir",
        help="Directory for KneadData output files (default: {output-dir}/PreprocessedData/kneaddata_output)",
    )
    preprocessing_group.add_argument(
        "--kraken-output-dir",
        help="Directory for Kraken2 output files (default: {output-dir}/PreprocessedData/kraken_output)",
    )
    preprocessing_group.add_argument(
        "--bracken-output-dir",
        help="Directory for Bracken output files (default: {output-dir}/PreprocessedData/bracken_output)",
    )

    # --- 4. Kraken/Bracken Processing Options ---
    kraken_group = parser.add_argument_group("Kraken/Bracken Processing Options")
    kraken_group.add_argument(
        "--kreport-dir", help="Directory containing raw Kraken2 report files (.kreport)"
    )
    kraken_group.add_argument(
        "--bracken-dir", help="Directory containing raw Bracken abundance files"
    )
    kraken_group.add_argument(
        "--taxonomic-level", 
        default="S", 
        choices=["D", "P", "C", "O", "F", "G", "S"],
        help="Taxonomic level for analysis (D=domain, P=phylum, C=class, O=order, F=family, G=genus, S=species)"
    )
    kraken_group.add_argument(
        "--threshold", 
        type=float, 
        default=10.0,
        help="Threshold for Bracken (minimum number of reads required)"
    )
    kraken_group.add_argument(
        "--skip-kraken", action="store_true", help="Skip Kraken2 report processing"
    )
    kraken_group.add_argument(
        "--skip-bracken", action="store_true", help="Skip Bracken abundance processing"
    )

    # --- 5. Downstream Analysis Options ---
    downstream_group = parser.add_argument_group("Downstream Analysis Options")
    downstream_group.add_argument(
        "--skip-downstream", action="store_true", help="Skip downstream analysis"
    )
    downstream_group.add_argument(
        "--group-col", default="Group", help="The column name to use for grouping in stats (default: 'Group')"
    )
    downstream_group.add_argument(
        "--min-abundance", type=float, default=0.01, 
        help="Minimum relative abundance threshold for inclusion in analysis (default: 0.01 = 1%%)"
    )
    downstream_group.add_argument(
        "--min-prevalence", type=float, default=0.1, 
        help="Minimum prevalence threshold (proportion of samples) for inclusion in analysis (default: 0.1 = 10%%)"
    )

    # --- 6. Differential Abundance Analysis Options ---
    diff_group = parser.add_argument_group("Differential Abundance Options")
    diff_group.add_argument(
        "--run-diff-abundance",
        action="store_true",
        help="Run differential abundance analysis using ANCOM, ALDEx2, and/or ANCOM-BC",
    )
    diff_group.add_argument(
        "--diff-methods",
        default="aldex2,ancom,ancom-bc",
        help="Comma-separated list of methods to use (default: aldex2,ancom,ancom-bc)",
    )

    # --- 7. Parallel Processing Options ---
    parallel_group = parser.add_argument_group("Parallel Processing Options")
    parallel_group.add_argument(
        "--threads-per-sample", type=int, default=1, help="Number of threads to use per sample"
    )
    parallel_group.add_argument(
        "--max-parallel", type=int, default=None, help="Maximum number of samples to process in parallel"
    )
    parallel_group.add_argument(
        "--use-parallel", action="store_true", help="Use parallel processing for preprocessing steps"
    )

    args = parser.parse_args()

    # Setup logging
    logger = setup_logger(log_file=args.log_file, log_level=getattr(logging, args.log_level))
    log_print("Starting Kraken Tools Pipeline", level="info")

    start_time = time.time()

    # If memory limit is specified, try to apply it
    if args.max_memory:
        success = limit_memory_usage(args.max_memory)
        if success:
            log_print(f"Set memory limit to {args.max_memory} MB", level="info")
        else:
            log_print("Failed to set memory limit, proceeding with unlimited memory", level="warning")

    # Make sure required directories are provided based on selected operations
    if not args.run_preprocessing and not args.list_files:
        if not args.kreport_dir and not args.bracken_dir:
            log_print("ERROR: At least one of --kreport-dir or --bracken-dir is required unless using --run-preprocessing or --list-files", level="error")
            sys.exit(1)

    # If only listing files, do that and exit
    if args.list_files:
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

    # Validate sample key first 
    samples, selected_columns = validate_sample_key(args.sample_key, no_interactive=args.no_interactive)

    # If running preprocessing, check installations and run the pipeline
    if args.run_preprocessing:
        if not args.input_fastq:
            log_print("ERROR: --input-fastq is required when using --run-preprocessing", level="error")
            sys.exit(1)

        # Check installations
        kneaddata_ok, kneaddata_msg = check_kneaddata_installation()
        if not kneaddata_ok:
            log_print(f"ERROR: KneadData not properly installed: {kneaddata_msg}", level="error")
            sys.exit(1)

        kraken_ok, kraken_msg = check_kraken_installation()
        if not kraken_ok:
            log_print(f"ERROR: Kraken2 not properly installed: {kraken_msg}", level="error")
            sys.exit(1)

        bracken_ok, bracken_msg = check_bracken_installation()
        if not bracken_ok:
            log_print(f"ERROR: Bracken not properly installed: {bracken_msg}", level="error")
            sys.exit(1)

        # Create preprocessing output directory
        preproc_dir = os.path.join(args.output_dir, "PreprocessedData")
        os.makedirs(preproc_dir, exist_ok=True)

        # Use specified output dirs or create defaults
        kneaddata_output_dir = args.kneaddata_output_dir
        if kneaddata_output_dir is None:
            kneaddata_output_dir = os.path.join(preproc_dir, "kneaddata_output")

        kraken_output_dir = args.kraken_output_dir
        if kraken_output_dir is None:
            kraken_output_dir = os.path.join(preproc_dir, "kraken_output")

        bracken_output_dir = args.bracken_output_dir
        if bracken_output_dir is None:
            bracken_output_dir = os.path.join(preproc_dir, "bracken_output")

        # Choose between regular or parallel processing
        if args.use_parallel:
            log_print("Using parallel preprocessing pipeline", level="info")
            # This function would be defined in preprocessing/pipeline.py
            preprocessing_results = run_preprocessing_pipeline_parallel(
                input_files=args.input_fastq,
                output_dir=preproc_dir,
                threads_per_sample=args.threads_per_sample,
                max_parallel=args.max_parallel,
                kneaddata_dbs=args.kneaddata_dbs,
                kraken_db=args.kraken_db,
                bracken_db=args.bracken_db,
                paired=args.paired,
                taxonomic_level=args.taxonomic_level,
                threshold=args.threshold,
                kneaddata_output_dir=kneaddata_output_dir,
                kraken_output_dir=kraken_output_dir,
                bracken_output_dir=bracken_output_dir,
                logger=logger,
            )
        else:
            log_print("Using standard preprocessing pipeline", level="info")
            preprocessing_results = run_preprocessing_pipeline(
                input_files=args.input_fastq,
                output_dir=preproc_dir,
                threads=args.threads,
                kneaddata_dbs=args.kneaddata_dbs,
                kraken_db=args.kraken_db,
                bracken_db=args.bracken_db,
                paired=args.paired,
                taxonomic_level=args.taxonomic_level,
                threshold=args.threshold,
                kneaddata_output_dir=kneaddata_output_dir,
                kraken_output_dir=kraken_output_dir,
                bracken_output_dir=bracken_output_dir,
                logger=logger,
            )

    # Here we would add the processing steps for Kraken/Bracken files
    # and downstream analysis similar to how it's done in humann3_tools
    # ...

    elapsed = time.time() - start_time
    hh, rr = divmod(elapsed, 3600)
    mm, ss = divmod(rr, 60)
    log_print(f"Pipeline finished in {int(hh)}h {int(mm)}m {int(ss)}s", level="info")


if __name__ == "__main__":
    main()
