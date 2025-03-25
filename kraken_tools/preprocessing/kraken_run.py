# kraken_tools/preprocessing/kraken_run.py
import os
import subprocess
import logging
import glob
import re
from kraken_tools.utils.cmd_utils import run_cmd_with_retry, CommandErrorType, CommandResult
from kraken_tools.utils.db_validation import validate_database, DatabaseType
from kraken_tools.logger import log_print

def check_kraken_installation():
    """Check if Kraken2 is installed and available."""
    try:
        result = subprocess.run(["kraken2", "--version"], 
                               capture_output=True, text=True, check=False)
        if result.returncode == 0:
            version_info = result.stdout.strip()
            # Extract version number
            match = re.search(r'version (\d+\.\d+\.\d+)', version_info)
            version = match.group(1) if match else "unknown"
            return True, f"Kraken2 version {version} found"
        return False, "Kraken2 command exists but returned an error"
    except FileNotFoundError:
        return False, "Kraken2 not found in PATH"

def process_single_sample_kraken(input_file, sample_id=None, output_dir=None, 
                               threads=1, kraken_db=None, paired_file=None, 
                               additional_options=None, logger=None):
    """Process a single sample with Kraken2."""
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    if sample_id is None:
        sample_id = os.path.basename(input_file).split('.')[0]
    
    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), "kraken_output", sample_id)
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Validate Kraken database before running
    db_validation = validate_database(kraken_db, DatabaseType.KRAKEN, logger)
    if not db_validation.success:
        logger.error(f"Kraken database validation failed for {kraken_db}. Cannot proceed with analysis.")
        return None
    
    # Build kraken2 command
    cmd = ["kraken2", "--db", kraken_db]
    
    # Add output files
    report_file = os.path.join(output_dir, f"{sample_id}.kreport")
    output_file = os.path.join(output_dir, f"{sample_id}.kraken")
    
    cmd.extend(["--report", report_file, "--output", output_file])
    
    # Add threads
    cmd.extend(["--threads", str(threads)])
    
    # Add input files (paired or single)
    if paired_file:
        cmd.extend(["--paired", input_file, paired_file])
    else:
        cmd.append(input_file)
    
    # Add additional options
    if additional_options:
        for key, value in additional_options.items():
            if value is True:
                cmd.append(f"--{key}")
            elif value is not None and value != "":
                cmd.extend([f"--{key}", str(value)])
    
    # Run Kraken2 with retry capability
    logger.info(f"Running Kraken2 for sample {sample_id}")
    
    # Set a reasonable timeout (20 minutes per thread as a rough baseline)
    timeout = threads * 1200  # 20 minutes per thread
    
    result = run_cmd_with_retry(cmd, timeout=timeout, exit_on_error=False, verbose=True)
    
    if not result.success:
        # Detailed error handling based on error type
        if result.error_type == CommandErrorType.DATABASE_ERROR:
            logger.error(
                f"Database error while running Kraken2 for sample {sample_id}. "
                f"Please verify your database is correctly built and accessible."
            )
        elif result.error_type == CommandErrorType.MEMORY_ERROR:
            logger.error(
                f"Out of memory while running Kraken2 for sample {sample_id}. "
                f"Try using fewer threads, a smaller database, or a machine with more RAM."
            )
        elif result.error_type == CommandErrorType.INPUT_ERROR:
            logger.error(
                f"Input file error while running Kraken2 for sample {sample_id}. "
                f"Check that input file exists and is readable: {input_file}"
            )
        else:
            logger.error(f"Kraken2 run failed for sample {sample_id}: {result.error_message}")
        
        # Check if any output was produced despite error
        if os.path.exists(report_file) and os.path.getsize(report_file) > 0:
            logger.warning(
                f"Kraken2 reported an error but produced output files. "
                f"The results may be incomplete but will be used."
            )
        else:
            return None
    
    # Verify the output files were created
    if not os.path.exists(report_file) or os.path.getsize(report_file) == 0:
        logger.error(f"Kraken2 did not produce a valid report file for sample {sample_id}")
        return None
    
    logger.info(f"Kraken2 completed for sample {sample_id}")
    
    # Return paths to output files
    return {
        'report': report_file,
        'output': output_file
    }

# Add this function at the module level (outside any other functions)
def paired_kraken_wrapper(input_file, sample_id=None, paired_file=None, **kwargs):
    """
    Wrapper to handle paired files correctly for Kraken.
    
    Args:
        input_file: Path to input file (R1 for paired-end)
        sample_id: Sample identifier
        paired_file: Path to paired file (R2)
        **kwargs: Additional arguments for process_single_sample_kraken
        
    Returns:
        Dictionary with output file paths or None if failed
    """
    # Call the appropriate function based on whether this is paired data
    if paired_file:
        return process_single_sample_kraken(
            input_file, 
            sample_id=sample_id, 
            paired_file=paired_file, 
            **kwargs
        )
    else:
        return process_single_sample_kraken(
            input_file,
            sample_id=sample_id,
            **kwargs
        )


def run_kraken_parallel(input_files, output_dir, threads=1, max_parallel=None,
                      kraken_db=None, paired=False, additional_options=None,
                      logger=None):
    """
    Run Kraken2 on multiple samples in parallel.
    
    Args:
        input_files: List of input FASTQ files
        output_dir: Base directory for outputs
        threads: Number of threads per sample
        max_parallel: Maximum number of samples to process in parallel
        kraken_db: Path to Kraken2 database
        paired: Whether input files are paired-end
        additional_options: Dict of additional Kraken2 options
        logger: Logger instance
        
    Returns:
        Dict mapping sample IDs to output files
    """
    from kraken_tools.utils.resource_utils import track_peak_memory
    from kraken_tools.preprocessing.parallel import run_parallel
    
  
    def _run_kraken_parallel(input_files, output_dir, threads_per_sample, max_parallel,
                           kraken_db, paired, additional_options, logger):
        if logger is None:
            logger = logging.getLogger('kraken_analysis')
        
        # First, validate the database once
        db_validation = validate_database(kraken_db, DatabaseType.KRAKEN, logger)
        if not db_validation.success:
            logger.error(f"Kraken database validation failed for {kraken_db}. Cannot proceed with analysis.")
            return {}
        
        # Create sample list
        sample_list = []
        
        if paired:
            # Ensure we have pairs of files
            if len(input_files) % 2 != 0:
                logger.error("Paired mode requires an even number of input files")
                return {}
            
            # Pair files in order (assuming they're in order: R1, R2, R1, R2, etc.)
            for i in range(0, len(input_files), 2):
                r1_file = input_files[i]
                r2_file = input_files[i+1]
                
                # Get sample name from the R1 file
                sample_name = extract_sample_name(r1_file)
                sample_list.append((sample_name, r1_file, r2_file))
        else:
            # Single-end reads
            for file in input_files:
                sample_name = os.path.basename(file).split('.')[0]
                sample_list.append((sample_name, file))
        
        # Prepare common arguments for all samples
        kwargs = {
            'output_dir': output_dir,
            'threads': threads_per_sample,
            'kraken_db': kraken_db,
            'additional_options': additional_options,
            'logger': logger
        }
        
        results = run_parallel(sample_list, paired_kraken_wrapper, 
                             max_workers=max_parallel, **kwargs)
        
        # Check for overall success rate
        total_samples = len(sample_list)
        successful_samples = len(results)
        success_rate = (successful_samples / total_samples) * 100 if total_samples > 0 else 0
        
        logger.info(f"Kraken2 completed with {successful_samples}/{total_samples} samples successful ({success_rate:.1f}%)")
        
        # If success rate is low, provide a warning
        if 0 < success_rate < 50:
            logger.warning(
                "Less than half of the samples were processed successfully. "
                "This may indicate issues with the database, input files, or system resources."
            )
        
        return results
    
    # Run the wrapped function
    return _run_kraken_parallel(
        input_files, output_dir, threads, max_parallel,
        kraken_db, paired, additional_options, logger
    )

def extract_sample_name(filename):
    """
    Extract sample name from filename using various heuristics.
    
    For paired files like "trimmed_read1_SampleName.fastq.gz", extracts "SampleName".
    
    Args:
        filename: Path to input file
        
    Returns:
        Extracted sample name
    """
    # Get just the filename without path
    basename = os.path.basename(filename)
    
    # Remove extensions like .fastq, .fastq.gz, .fq, .fq.gz
    for ext in ['.fastq.gz', '.fq.gz', '.fastq', '.fq']:
        if basename.lower().endswith(ext):
            basename = basename[:-len(ext)]
            break
    
    # Try to identify common paired-end patterns
    patterns = [
        # Pattern: trimmed_read1_SampleName or trimmed_read2_SampleName
        r'trimmed_read[12]_(.+)',
        # Pattern: SampleName_R1 or SampleName_R2
        r'(.+)_R[12]',
        # Pattern: SampleName_1 or SampleName_2
        r'(.+)_[12]',
        # Pattern: SampleName.R1 or SampleName.R2
        r'(.+)\.R[12]'
    ]
    
    import re
    for pattern in patterns:
        match = re.match(pattern, basename)
        if match:
            return match.group(1)
    
    # If no pattern matches, fall back to using parts
    parts = basename.split('_')
    if len(parts) > 2:
        # For "trimmed_read1_SampleName", return "SampleName"
        return parts[2]
    elif len(parts) == 2:
        # For patterns like "SampleName_R1", return "SampleName"
        if parts[1].startswith('R') or parts[1].isdigit():
            return parts[0]
        else:
            return parts[1]
    
    # Fallback to original behavior if nothing else works
    return parts[0]

def run_kraken(input_files, output_dir, threads=1, kraken_db=None,
              paired=False, additional_options=None, logger=None):
    """
    Run Kraken2 on input sequence files.
    
    Args:
        input_files: List of input FASTQ files
        output_dir: Directory for Kraken2 output
        threads: Number of threads to use
        kraken_db: Path to Kraken2 database
        paired: Whether input files are paired
        additional_options: Dict of additional Kraken2 options
        logger: Logger instance
        
    Returns:
        Dict mapping sample IDs to output files
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    os.makedirs(output_dir, exist_ok=True)
    
    # First, validate the database once
    db_validation = validate_database(kraken_db, DatabaseType.KRAKEN, logger)
    if not db_validation.success:
        logger.error(f"Kraken database validation failed for {kraken_db}. Cannot proceed with analysis.")
        return {}
    
    results = {}
    
    # Process files based on whether they're paired
    if paired:
        # Ensure we have pairs of files
        if len(input_files) % 2 != 0:
            logger.error("Paired mode requires an even number of input files")
            return {}
        
        # Process pairs
        for i in range(0, len(input_files), 2):
            r1_file = input_files[i]
            r2_file = input_files[i+1]
            
            # Get sample name from the file
            sample_name = extract_sample_name(r1_file)
            logger.info(f"Processing paired files for sample {sample_name}")
            
            sample_output_dir = os.path.join(output_dir, sample_name)
            os.makedirs(sample_output_dir, exist_ok=True)
            
            result = process_single_sample_kraken(
                r1_file,
                sample_id=sample_name,
                output_dir=sample_output_dir,
                threads=threads,
                kraken_db=kraken_db,
                paired_file=r2_file,
                additional_options=additional_options,
                logger=logger
            )
            
            if result:
                results[sample_name] = result
    else:
        # Process single-end files
        for file in input_files:
            sample_name = extract_sample_name(file)
            logger.info(f"Processing file for sample {sample_name}")
            
            sample_output_dir = os.path.join(output_dir, sample_name)
            os.makedirs(sample_output_dir, exist_ok=True)
            
            result = process_single_sample_kraken(
                file,
                sample_id=sample_name,
                output_dir=sample_output_dir,
                threads=threads,
                kraken_db=kraken_db,
                additional_options=additional_options,
                logger=logger
            )
            
            if result:
                results[sample_name] = result
    
    # Check for overall success rate
    total_samples = len(input_files) // 2 if paired else len(input_files)
    successful_samples = len(results)
    success_rate = (successful_samples / total_samples) * 100 if total_samples > 0 else 0
    
    logger.info(f"Kraken2 completed with {successful_samples}/{total_samples} samples successful ({success_rate:.1f}%)")
    
    return results
