# kraken_tools/preprocessing/kraken_run.py
import os
import subprocess
import logging
from kraken_tools.utils.cmd_utils import run_cmd
from kraken_tools.logger import log_print

def check_kraken_installation():
    """Check if Kraken2 is installed and available."""
    try:
        result = subprocess.run(["kraken2", "--version"], 
                               capture_output=True, text=True, check=False)
        if result.returncode == 0:
            return True, result.stdout.strip()
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
    
    # Build kraken2 command
    cmd = ["kraken2", "--db", kraken_db if kraken_db else "none"]
    
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
    
    # Run Kraken2
    logger.info(f"Running Kraken2 for sample {sample_id}")
    success = run_cmd(cmd, exit_on_error=False)
    
    if not success:
        logger.error(f"Kraken2 run failed for sample {sample_id}")
        return None
    
    logger.info(f"Kraken2 completed for sample {sample_id}")
    
    # Return paths to output files
    return {
        'report': report_file,
        'output': output_file
    }

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
    
    @track_peak_memory
    def _run_kraken_parallel(input_files, output_dir, threads_per_sample, max_parallel,
                           kraken_db, paired, additional_options, logger):
        if logger is None:
            logger = logging.getLogger('kraken_analysis')
        
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
                sample_name = os.path.basename(r1_file).split('_')[0]
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
        
        # Define wrapper function to handle paired/unpaired correctly
        def paired_kraken_wrapper(input_file, sample_id=None, paired_file=None, **kwargs):
            # This function is defined in preprocessing/parallel.py
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
        
        # Run in parallel
        results = run_parallel(sample_list, paired_kraken_wrapper, 
                             max_workers=max_parallel, **kwargs)
        
        return results
    
    # Run the wrapped function
    return _run_kraken_parallel(
        input_files, output_dir, threads, max_parallel,
        kraken_db, paired, additional_options, logger
    )

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
            sample_name = os.path.basename(r1_file).split('_')[0]
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
            sample_name = os.path.basename(file).split('.')[0]
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
    
    logger.info(f"Kraken2 completed for {len(results)} samples")
    return results
