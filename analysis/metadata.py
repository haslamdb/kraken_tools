# kraken_tools/analysis/metadata.py
import pandas as pd
import logging
import traceback

def read_and_process_metadata(sample_key, logger):
    """
    Read and process sample metadata file.
    
    Args:
        sample_key: Path to the sample key CSV file
        logger: Logger instance for logging
        
    Returns:
        DataFrame containing the sample metadata
    """
    try:
        df = pd.read_csv(sample_key)
        logger.info(f"Loaded sample key with {len(df)} rows, columns: {list(df.columns)}")
        
        # Try to automatically identify the sample ID column if not explicitly named "SampleName"
        if "SampleName" not in df.columns:
            potential_id_cols = ["SampleID", "Sample_ID", "Sample", "ID", "sample_id", "sample_name"]
            for col in potential_id_cols:
                if col in df.columns:
                    logger.info(f"Renaming column '{col}' to 'SampleName'")
                    df = df.rename(columns={col: "SampleName"})
                    break
        
        # Check if we have a SampleName column now
        if "SampleName" not in df.columns:
            logger.warning(f"Could not identify sample ID column. Using first column as SampleName: {df.columns[0]}")
            df = df.rename(columns={df.columns[0]: "SampleName"})
        
        return df
    except Exception as e:
        logger.error(f"Error reading sample key: {str(e)}")
        logger.error(traceback.format_exc())
        raise RuntimeError(f"Failed to process metadata: {e}")
