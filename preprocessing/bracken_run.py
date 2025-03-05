# kraken_tools/preprocessing/bracken_run.py
import os
import subprocess
import logging
from kraken_tools.utils.cmd_utils import run_cmd
from kraken_tools.logger import log_print

def check_bracken_installation():
    """Check if Bracken is installed and available."""
    try:
        # Try to execute est_abundance.py (the main Bracken script)
        result = subprocess.run(["est_abundance.py", "--help"], 
                               capture_output=True, text=True, check=False)
        if result.returncode == 0:
            return True, "Bracken (est_abundance.py) is installed"
        return False, "Bracken command exists but returned an error"
    except FileNotFoundError:
        # Try an alternative location, as it might be installed in a different path
        try:
            result = subprocess.run(["bracken", "--help"], 
                                   capture_output=True, text=True, check=False)
            if result.returncode == 0:
                return True, "Bracken is installed"
            return False, "Bracken command exists but returned an error"
        except FileNotFoundError:
            return False, "Bracken (est_abundance.py or bracken) not found in PATH"

def process_single_sample_bracken(kreport_file, sample_id=None, output_dir=None, 
                                 bracken_db=None, taxonomic_level="S", 
                                 threshold=10, threads=1, additional_options=None, 
                                 logger=None):
    """Process a single Kraken2 report with Bracken for abundance estimation."""
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    if sample_id is None:
        sample_id = os.path.basename(kreport_file).split('.')[0]
    
    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), "bracken_output", sample_id)
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Check if input file exists
    if not os.path.exists(kreport_file):
        logger.error(f"Kraken report file not found: {kreport_file}")
        return None
    
    # Build est_abundance.py command
    cmd = ["est_abundance.py"]
    
    # Add input file
    cmd.extend(["-i", kreport_file])
    
    # Add database
    cmd.extend(["-k", bracken_db])
    
    # Add taxonomic level
    cmd.extend(["-l", taxonomic_level])
    
    # Add threshold
    cmd.extend(["-t", str(threshold)])
    
    # Add threads
    cmd.extend(["-t", str(threads)])
    
    # Add output file
    output_file = os.path.join(output_dir, f"{sample_id}_{taxonomic_level.lower()}_abundance.txt")
    cmd.extend(["-o", output_file])
    
    # Add additional options
    if additional_options:
        for key, value in additional_options.items():
            if value is True:
                cmd.append(f"--{key}")
            elif value is not None and value != "":
                cmd.extend([f"--{key}", str(value)])
    
    # Run Bracken
    logger.info(f"Running Bracken (est_abundance.py) for sample {sample_id} at {taxonomic_level} level")
    success = run_cmd(cmd, exit_on_error=False)
    
    if not success:
        logger.error(f"Bracken run failed for sample {sample_id}")
        return None
    
    logger.info(f"Bracken completed for sample {sample_id}")
    
    # Also generate MPA-style report if possible
    mpa_file = None
    try:
        mpa_output = os.path.join(output_dir, f"{sample_id}.mpa.report")
        mpa_cmd = ["kreport2mpa.py", "-r", kreport_file, "-o", mpa_output, "--no-intermediate-ranks"]
        
        logger.info(f"Generating MPA report for sample {sample_id}")
        if run_cmd(mpa_cmd, exit_on_error=False):
            mpa_file = mpa_output
            logger.info(f"MPA report generated for sample {sample_id}")
        else:
            logger.warning(f"Failed to generate MPA report for sample {sample_id}")
    except Exception as e:
        logger.warning(f"Error generating MPA report: {e}")
    
    # Return paths to output files
    return {
        'abundance': output_file,
        'mpa': mpa_file
    }

def run_bracken_parallel(kreport_files, output_dir, threads=1, max_parallel=None,
                       bracken_db=None, taxonomic_level="S", threshold=10, 
                       additional_options=None, logger=None):
    """
    Run Bracken on multiple Kraken reports in parallel.
    
    Args:
        kreport_files: List of Kraken report files or Dict mapping sample IDs to report files
        output_dir: Base directory for outputs
        threads: Number of threads per sample
        max_parallel: Maximum number of samples to process in parallel
        bracken_db: Path to Bracken database
        taxonomic_level: Taxonomic level for abundance estimation
        threshold: Minimum number of reads for abundance estimation
        additional_options: Dict of additional Bracken options
        logger: Logger instance
        
    Returns:
        Dict mapping sample IDs to output files
    """
    from kraken_tools.utils.resource_utils import track_peak_memory
    from kraken_tools.preprocessing.parallel import run_parallel
    
    @track_peak_memory
    def _run_bracken_parallel(kreport_files, output_dir, threads_per_sample, max_parallel,
                            bracken_db, taxonomic_level, threshold, additional_options, 
                            logger):
        if logger is None:
            logger = logging.getLogger('kraken_analysis')
        
        # Create sample list
        sample_list = []
        
        # Check if input is a dict or list
        if isinstance(kreport_files, dict):
            for sample_id, kreport_file in kreport_files.items():
                sample_list.append((sample_id, kreport_file))
        else:
            # Process list of report files
            for file in kreport_files:
                sample_name = os.path.basename(file).split('.')[0]
                if sample_name.endswith('.kreport'):
                    sample_name = sample_name[:-8]  # Remove .kreport
                sample_list.append((sample_name, file))
        
        # Prepare common arguments for all samples
        kwargs = {
            'output_dir': output_dir,
            'bracken_db': bracken_db,
            'taxonomic_level': taxonomic_level,
            'threshold': threshold,
            'threads': threads_per_sample,
            'additional_options': additional_options,
            'logger': logger
        }
        
        # Run in parallel
        results = run_parallel(sample_list, process_single_sample_bracken, 
                             max_workers=max_parallel, **kwargs)
        
        return results
    
    # Run the wrapped function
    return _run_bracken_parallel(
        kreport_files, output_dir, threads, max_parallel,
        bracken_db, taxonomic_level, threshold, additional_options, logger
    )

def run_bracken(kreport_files, output_dir, threads=1, bracken_db=None, 
              taxonomic_level="S", threshold=10, additional_options=None, logger=None):
    """
    Run Bracken on Kraken2 report files.
    
    Args:
        kreport_files: List of Kraken report files or Dict mapping sample IDs to report files
        output_dir: Directory for Bracken output
        threads: Number of threads to use
        bracken_db: Path to Bracken database
        taxonomic_level: Taxonomic level for abundance estimation
        threshold: Minimum number of reads for abundance estimation
        additional_options: Dict of additional Bracken options
        logger: Logger instance
        
    Returns:
        Dict mapping sample IDs to output files
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    os.makedirs(output_dir, exist_ok=True)
    
    results = {}
    
    # Check if input is a dict or list
    if isinstance(kreport_files, dict):
        input_files = [(sample_id, file_path) for sample_id, file_path in kreport_files.items()]
    else:
        # Process list of report files
        input_files = []
        for file in kreport_files:
            sample_name = os.path.basename(file).split('.')[0]
            if sample_name.endswith('.kreport'):
                sample_name = sample_name[:-8]  # Remove .kreport
            input_files.append((sample_name, file))
    
    # Process each report file
    for sample_id, kreport_file in input_files:
        logger.info(f"Processing Kraken report for sample {sample_id}")
        
        sample_output_dir = os.path.join(output_dir, sample_id)
        os.makedirs(sample_output_dir, exist_ok=True)
        
        result = process_single_sample_bracken(
            kreport_file,
            sample_id=sample_id,
            output_dir=sample_output_dir,
            bracken_db=bracken_db,
            taxonomic_level=taxonomic_level,
            threshold=threshold,
            threads=threads,
            additional_options=additional_options,
            logger=logger
        )
        
        if result:
            results[sample_id] = result
    
    logger.info(f"Bracken completed for {len(results)} samples")
    return results
