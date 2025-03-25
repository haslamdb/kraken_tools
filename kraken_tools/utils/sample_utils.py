# kraken_tools/utils/sample_utils.py
import os
import sys
import csv
import logging

from kraken_tools.utils.file_utils import check_file_exists
from kraken_tools.logger import log_print

def validate_sample_key_noninteractive(sample_key_file):
    """Simpler version of reading sample key in non-interactive mode."""
    if not check_file_exists(sample_key_file, "Sample key"):
        sys.exit(1)
    try:
        with open(sample_key_file, "r", encoding="utf-8-sig") as csvfile:
            reader = csv.DictReader(csvfile)
            csv_columns = reader.fieldnames if reader.fieldnames else []
            if not csv_columns:
                log_print("ERROR: Sample key has no columns", level='error')
                sys.exit(1)
            
            # Try to find common sample ID
            sample_id_col = None
            common_id_names = ["SampleName", "Sample", "SampleID", "Sample_ID", "sample_name", "sample_id"]
            for col in common_id_names:
                if col in csv_columns:
                    sample_id_col = col
                    break
            
            # If still none found, use the first column
            if not sample_id_col:
                sample_id_col = csv_columns[0]
            
            # Gather sample IDs
            csvfile.seek(0)
            next(reader)
            samples = []
            for row in reader:
                if row.get(sample_id_col):
                    samples.append(row[sample_id_col])
            return samples
    except Exception as e:
        log_print(f"ERROR reading sample key (non-interactive): {str(e)}", level='error')
        sys.exit(1)


def validate_sample_key(sample_key_file, no_interactive=False):
    """
    Validate the sample key CSV file. 
    If no_interactive=True, do a simpler check; else attempt user interaction.
    Returns (samples, selected_columns) or (samples, None).
    """
    if not check_file_exists(sample_key_file, "Sample key"):
        sys.exit(1)
    
    # Non-interactive version
    if no_interactive:
        return validate_sample_key_noninteractive(sample_key_file), None

    # -- Interactive version (no outer try here) --
    with open(sample_key_file, "r", encoding="utf-8-sig") as csvfile:
        reader = csv.DictReader(csvfile)
        csv_columns = reader.fieldnames if reader.fieldnames else []
        if not csv_columns:
            log_print(f"ERROR: Sample key has no columns: {sample_key_file}", level='error')
            sys.exit(1)
        
        log_print(f"\nAvailable columns in sample key ({len(csv_columns)}):")
        for i, col in enumerate(csv_columns, 1):
            log_print(f"  {i}. {col}")
        
        # Try to find a sample ID column
        sample_id_col = None
        common_id_names = ["SampleName", "Sample", "SampleID", "Sample_ID", "sample_name", "sample_id"]
        for col in common_id_names:
            if col in csv_columns:
                sample_id_col = col
                log_print(f"\nAutomatic selection: '{col}' as sample identifier", level='info')
                break
        
        # If no common name found, prompt user
        if not sample_id_col:
            while True:
                try:
                    log_print("\nWhich column has the sample identifiers?", level='info')
                    selection = input("Enter column number or name: ")
                    # Try as number
                    try:
                        idx = int(selection) - 1
                        if 0 <= idx < len(csv_columns):
                            sample_id_col = csv_columns[idx]
                            break
                        else:
                            log_print(f"Error: number must be between 1 and {len(csv_columns)}", level='error')
                    except ValueError:
                        # Try as column name
                        if selection in csv_columns:
                            sample_id_col = selection
                            break
                        else:
                            log_print(f"Error: '{selection}' is not valid", level='error')
                except KeyboardInterrupt:
                    log_print("Operation aborted by user", level='warning')
                    sys.exit(1)
        
        # Optional grouping columns
        group_cols = {}
        remaining_cols = [c for c in csv_columns if c != sample_id_col]
        if remaining_cols:
            log_print("\nSelect columns for grouping (optional). Press Enter to skip:", level='info')
            for i, col in enumerate(remaining_cols, 1):
                log_print(f"  {i}. {col}")
            
            log_print("Example: '1,3,5' or 'Group,Treatment,Site'", level='info')
            try:
                selection = input("Selection: ")
                if selection.strip():
                    items = [it.strip() for it in selection.split(',')]
                    for it in items:
                        try:
                            idx = int(it) - 1
                            if 0 <= idx < len(remaining_cols):
                                col_name = remaining_cols[idx]
                                group_cols[col_name] = col_name
                        except ValueError:
                            if it in remaining_cols:
                                group_cols[it] = it
            except Exception as e:
                log_print(f"ERROR reading grouping columns: {str(e)}", level='error')
                sys.exit(1)

        # After optional grouping, read sample IDs
        csvfile.seek(0)
        next(reader)
        samples = []
        for row in reader:
            if row.get(sample_id_col):
                samples.append(row[sample_id_col])
        
        if not samples:
            log_print(f"ERROR: No samples found in column '{sample_id_col}'", level='error')
            sys.exit(1)
        
        selected_columns = {
            "sample_id": sample_id_col,
            "grouping": group_cols
        }
        return samples, selected_columns


def check_input_files_exist(samples, kreport_dir, bracken_dir):
    """
    Check if all required input files exist for each sample.
    Returns:
        valid_kreport_samples: list of (sample, kreport_file)
        valid_bracken_samples: list of (sample, bracken_file)
    """
    missing_kreport_files = []
    missing_bracken_files = []
    valid_kreport_samples = []
    valid_bracken_samples = []
    
    kreport_patterns = [
        "{sample}.kreport",
        "{sample}.kreport.txt",
        "{sample}_kreport.txt",
        "kreport_{sample}.txt"
    ]
    bracken_patterns = [
        "{sample}.bracken",
        "{sample}_abundance.txt",
        "{sample}_bracken_abundance.txt",
        "bracken_{sample}.txt"
    ]
    
    logger = logging.getLogger('kraken_analysis')
    logger.info(f"Checking kreport file patterns: {kreport_patterns}")
    logger.info(f"Checking bracken file patterns: {bracken_patterns}")
    
    for sample in samples:
        # Check kreport
        found_kreport = False
        for pattern in kreport_patterns:
            kreport_file = os.path.join(kreport_dir, pattern.format(sample=sample))
            if os.path.isfile(kreport_file):
                valid_kreport_samples.append((sample, kreport_file))
                found_kreport = True
                break
        if not found_kreport:
            missing_kreport_files.append(sample)
        
        # Check bracken
        found_bracken = False
        for pattern in bracken_patterns:
            bracken_file = os.path.join(bracken_dir, pattern.format(sample=sample))
            if os.path.isfile(bracken_file):
                valid_bracken_samples.append((sample, bracken_file))
                found_bracken = True
                break
        if not found_bracken:
            missing_bracken_files.append(sample)
    
    if missing_kreport_files:
        log_print(f"WARNING: {len(missing_kreport_files)} samples missing kreport files", level='warning')
        if len(missing_kreport_files) <= 10:
            log_print(f"Missing: {missing_kreport_files}", level='warning')
    
    if missing_bracken_files:
        log_print(f"WARNING: {len(missing_bracken_files)} samples missing bracken files", level='warning')
        if len(missing_bracken_files) <= 10:
            log_print(f"Missing: {missing_bracken_files}", level='warning')
    
    if not valid_kreport_samples and not valid_bracken_samples:
        log_print("ERROR: No valid kreport or bracken files found for any samples. Exiting.", level='error')
        sys.exit(1)
    
    log_print(f"Found valid kreport files for {len(valid_kreport_samples)} / {len(samples)} samples", level='info')
    log_print(f"Found valid bracken files for {len(valid_bracken_samples)} / {len(samples)} samples", level='info')
    return valid_kreport_samples, valid_bracken_samples