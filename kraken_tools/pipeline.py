# kraken_tools/preprocessing/pipeline.py
import os
import logging
from kraken_tools.preprocessing.kneaddata import run_kneaddata, run_kneaddata_parallel
from kraken_tools.preprocessing.kraken_run import run_kraken, run_kraken_parallel
from kraken_tools.preprocessing.bracken_run import run_bracken, run_bracken_parallel
from kraken_tools.logger import log_print
from kraken_tools.utils.resource_utils import track_peak_memory

def run_preprocessing_pipeline(
    input_files,
    output_dir,
    threads=1,
    kneaddata_dbs=None,
    kraken_db=None,
    bracken_db=None,
    taxonomic_level="S",
    threshold=10,
    paired=False,
    kneaddata_options=None,
    kraken_options=None,
    bracken_options=None,
    kneaddata_output_dir=None,
    kraken_output_dir=None,
    bracken_output_dir=None,
    logger=None
):
    """
    Run the full preprocessing pipeline: KneadData → Kraken2 → Bracken.

    Args:
        input_files: List of input FASTQ files
        output_dir: Base directory for outputs
        threads: Number of threads to use
        kneaddata_dbs: Path(s) to KneadData reference database(s)
        kraken_db: Path to Kraken2 database
        bracken_db: Path to Bracken database
        taxonomic_level: Taxonomic level for Bracken (D, P, C, O, F, G, S)
        threshold: Threshold for Bracken (minimum number of reads)
        paired: Whether input files are paired-end
        kneaddata_options: Dict of additional KneadData options
        kraken_options: Dict of additional Kraken2 options
        bracken_options: Dict of additional Bracken options
        kneaddata_output_dir: Custom directory for KneadData outputs
        kraken_output_dir: Custom directory for Kraken2 outputs
        bracken_output_dir: Custom directory for Bracken outputs
        logger: Logger instance
        
    Returns:
        Dict with results from each step
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    # 1. Create output directories
    os.makedirs(output_dir, exist_ok=True)
    
    if kneaddata_output_dir is None:
        kneaddata_output_dir = os.path.join(output_dir, "kneaddata_output")
        
    if kraken_output_dir is None:
        kraken_output_dir = os.path.join(output_dir, "kraken_output")
        
    if bracken_output_dir is None:
        bracken_output_dir = os.path.join(output_dir, "bracken_output")
    
    os.makedirs(kneaddata_output_dir, exist_ok=True)
    os.makedirs(kraken_output_dir, exist_ok=True)
    os.makedirs(bracken_output_dir, exist_ok=True)
    
    logger.info(f"Starting preprocessing pipeline with {len(input_files)} input files")
    
    # 2. Run KneadData
    logger.info("Starting KneadData step...")
    kneaddata_files = run_kneaddata(
        input_files=input_files,
        output_dir=kneaddata_output_dir,
        threads=threads,
        reference_dbs=kneaddata_dbs,
        paired=paired,
        additional_options=kneaddata_options,
        logger=logger
    )
    
    if not kneaddata_files:
        logger.error("KneadData step failed; stopping pipeline")
        return None
    
    logger.info(f"KneadData produced {len(kneaddata_files)} output files")
    
    # 3. Run Kraken2
    logger.info("Starting Kraken2 step...")
    kraken_results = run_kraken(
        input_files=kneaddata_files,
        output_dir=kraken_output_dir,
        threads=threads,
        kraken_db=kraken_db,
        paired=False,  # KneadData outputs are always single-end
        additional_options=kraken_options,
        logger=logger
    )
    
    if not kraken_results:
        logger.error("Kraken2 step failed; stopping pipeline")
        return {
            'kneaddata_files': kneaddata_files
        }
    
    logger.info(f"Kraken2 completed for {len(kraken_results)} samples")
    
    # 4. Run Bracken
    logger.info("Starting Bracken step...")
    kreport_files = {sample_id: results['report'] for sample_id, results in kraken_results.items()}
    
    bracken_results = run_bracken(
        kreport_files=kreport_files,
        output_dir=bracken_output_dir,
        threads=threads,
        bracken_db=bracken_db,
        taxonomic_level=taxonomic_level,
        threshold=threshold,
        additional_options=bracken_options,
        logger=logger
    )
    
    if not bracken_results:
        logger.error("Bracken step failed")
        return {
            'kneaddata_files': kneaddata_files,
            'kraken_results': kraken_results
        }
    
    logger.info(f"Bracken completed for {len(bracken_results)} samples")
    
    # 5. Return results
    return {
        'kneaddata_files': kneaddata_files,
        'kraken_results': kraken_results,
        'bracken_results': bracken_results
    }

@track_peak_memory
def run_preprocessing_pipeline_parallel(
    input_files,
    output_dir,
    threads_per_sample=1,
    max_parallel=None,
    kneaddata_dbs=None,
    kraken_db=None,
    bracken_db=None,
    taxonomic_level="S",
    threshold=10,
    paired=False,
    kneaddata_options=None,
    kraken_options=None,
    bracken_options=None,
    kneaddata_output_dir=None,
    kraken_output_dir=None,
    bracken_output_dir=None,
    logger=None
):
    """
    Run the full preprocessing pipeline in parallel:
    KneadData → Kraken2 → Bracken.
    
    Args:
        Same as run_preprocessing_pipeline, with additional parallelization parameters
        
    Returns:
        Dict with results from each step
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    # 1. Create output directories
    os.makedirs(output_dir, exist_ok=True)
    
    if kneaddata_output_dir is None:
        kneaddata_output_dir = os.path.join(output_dir, "kneaddata_output")
        
    if kraken_output_dir is None:
        kraken_output_dir = os.path.join(output_dir, "kraken_output")
        
    if bracken_output_dir is None:
        bracken_output_dir = os.path.join(output_dir, "bracken_output")
    
    os.makedirs(kneaddata_output_dir, exist_ok=True)
    os.makedirs(kraken_output_dir, exist_ok=True)
    os.makedirs(bracken_output_dir, exist_ok=True)
    
    logger.info(f"Starting parallel preprocessing pipeline with {len(input_files)} input files")
    
    # 2. Run KneadData in parallel
    logger.info("Starting KneadData step in parallel...")
    kneaddata_results = run_kneaddata_parallel(
        input_files=input_files,
        output_dir=kneaddata_output_dir,
        threads=threads_per_sample,
        max_parallel=max_parallel,
        reference_dbs=kneaddata_dbs,
        paired=paired,
        additional_options=kneaddata_options,
        logger=logger
    )
    
    if not kneaddata_results:
        logger.error("KneadData step failed; stopping pipeline")
        return None
    
    # Collect KneadData outputs
    kneaddata_files = []
    for sample_files in kneaddata_results.values():
        if isinstance(sample_files, list):
            kneaddata_files.extend(sample_files)
    
    if not kneaddata_files:
        logger.error("No KneadData output files found")
        return None
        
    logger.info(f"KneadData produced {len(kneaddata_files)} output files")
    
    # 3. Run Kraken2 in parallel
    logger.info("Starting Kraken2 step in parallel...")
    kraken_results = run_kraken_parallel(
        input_files=kneaddata_files,
        output_dir=kraken_output_dir,
        threads=threads_per_sample,
        max_parallel=max_parallel,
        kraken_db=kraken_db,
        paired=False,  # KneadData outputs are always single-end
        additional_options=kraken_options,
        logger=logger
    )
    
    if not kraken_results:
        logger.error("Kraken2 step failed; stopping pipeline")
        return {
            'kneaddata_files': kneaddata_files
        }
    
    logger.info(f"Kraken2 completed for {len(kraken_results)} samples")
    
    # 4. Run Bracken in parallel
    logger.info("Starting Bracken step in parallel...")
    kreport_files = {}
    for sample_id, results in kraken_results.items():
        if 'report' in results:
            kreport_files[sample_id] = results['report']
    
    bracken_results = run_bracken_parallel(
        kreport_files=kreport_files,
        output_dir=bracken_output_dir,
        threads=threads_per_sample,
        max_parallel=max_parallel,
        bracken_db=bracken_db,
        taxonomic_level=taxonomic_level,
        threshold=threshold,
        additional_options=bracken_options,
        logger=logger
    )
    
    if not bracken_results:
        logger.error("Bracken step failed")
        return {
            'kneaddata_files': kneaddata_files,
            'kraken_results': kraken_results
        }
    
    logger.info(f"Bracken completed for {len(bracken_results)} samples")
    
    # 5. Return results
    return {
        'kneaddata_files': kneaddata_files,
        'kraken_results': kraken_results,
        'bracken_results': bracken_results
    }
