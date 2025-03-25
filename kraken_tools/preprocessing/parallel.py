# kraken_tools/preprocessing/parallel.py
import os
import time
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

def process_sample_parallel(sample_tuple, function, **kwargs):
    """
    Process a single sample with the provided function.
    
    Args:
        sample_tuple: Tuple of (sample_id, file_path) or (sample_id, r1_file, r2_file)
        function: Function to run on the sample
        **kwargs: Additional arguments to pass to the function
        
    Returns:
        Tuple of (sample_id, result)
    """
    logger = logging.getLogger('kraken_analysis')
    logger.debug(f"Processing sample tuple: {sample_tuple}")
    
    # Extract the sample_id which is always the first element
    sample_id = sample_tuple[0]
    
    # Handle different tuple formats based on length
    if len(sample_tuple) == 3:
        # This is a 3-value tuple (sample_id, r1_file, r2_file) for paired data
        _, r1_file, r2_file = sample_tuple
        
        start_time = time.time()
        logger.info(f"Started processing sample {sample_id}")
        
        # Call the function with r1_file as input_file and r2_file as paired_file
        result = function(r1_file, sample_id=sample_id, paired_file=r2_file, **kwargs)
        
        elapsed = time.time() - start_time
        logger.info(f"Finished processing sample {sample_id} in {elapsed:.2f} seconds")
        
        return sample_id, result
    
    elif len(sample_tuple) == 2:
        # This is a 2-value tuple (sample_id, file_path)
        _, file_path = sample_tuple
        
        start_time = time.time()
        logger.info(f"Started processing sample {sample_id}")
        
        result = function(file_path, sample_id=sample_id, **kwargs)
        
        elapsed = time.time() - start_time
        logger.info(f"Finished processing sample {sample_id} in {elapsed:.2f} seconds")
        
        return sample_id, result
    
    else:
        # Invalid tuple format
        logger.error(f"Invalid sample tuple format: {sample_tuple}")
        return sample_id, None

def run_parallel(sample_list, function, max_workers=None, **kwargs):
    """
    Run a function on multiple samples in parallel with progress bar.
    
    Args:
        sample_list: List of tuples with sample information (could be 2 or 3 elements)
        function: Function to run on each sample
        max_workers: Maximum number of parallel processes (None = CPU count)
        **kwargs: Additional arguments to pass to the function
        
    Returns:
        Dictionary mapping sample_ids to results
    """
    logger = logging.getLogger('kraken_analysis')
    logger.info(f"Starting parallel processing of {len(sample_list)} samples with {max_workers} workers")
    
    results = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit tasks with proper handling for different tuple formats
        future_to_sample = {}
        for sample_tuple in sample_list:
            # Extract sample_id which should always be the first element
            if len(sample_tuple) >= 1:
                sample_id = sample_tuple[0]
                future = executor.submit(process_sample_parallel, sample_tuple, function, **kwargs)
                future_to_sample[future] = sample_id
            else:
                logger.error(f"Invalid sample tuple format: {sample_tuple}")
        
        # Process results with progress bar
        for future in tqdm(as_completed(future_to_sample), total=len(future_to_sample), 
                          desc="Processing samples", unit="sample"):
            sample_id = future_to_sample[future]
            try:
                sample_id, result = future.result()
                if result is not None:
                    results[sample_id] = result
                    logger.info(f"Successfully processed sample {sample_id}")
                else:
                    logger.error(f"Failed to process sample {sample_id}")
            except Exception as e:
                logger.error(f"Error processing sample {sample_id}: {str(e)}")
                logger.debug(f"Error details:", exc_info=True)
    
    logger.info(f"Completed parallel processing. Successfully processed {len(results)} of {len(sample_list)} samples")
    return results
