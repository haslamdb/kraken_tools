# kraken_tools/utils/db_validation.py
import os
import logging
import hashlib
import glob
import json
from enum import Enum

class DatabaseType(Enum):
    """Enum for different database types used in kraken-tools."""
    KRAKEN = "kraken"
    BRACKEN = "bracken"
    KNEADDATA = "kneaddata"

class DatabaseValidationResult:
    """Object to store the result of database validation."""
    def __init__(self, success=False, db_type=None, 
                 error_message=None, recommendations=None, 
                 files_checked=None, db_info=None):
        self.success = success
        self.db_type = db_type
        self.error_message = error_message
        self.recommendations = recommendations or []
        self.files_checked = files_checked or []
        self.db_info = db_info or {}
        
    def __bool__(self):
        return self.success

def get_directory_size(path):
    """Calculate the total size of a directory in bytes."""
    total_size = 0
    for dirpath, _, filenames in os.walk(path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            if os.path.exists(fp) and not os.path.islink(fp):
                total_size += os.path.getsize(fp)
    return total_size

def format_size(size_bytes):
    """Format size in bytes to human-readable format."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.2f} PB"

def validate_kraken_db(db_path, logger=None):
    """
    Validate a Kraken2 database directory.
    
    Args:
        db_path: Path to Kraken2 database directory
        logger: Logger instance
        
    Returns:
        DatabaseValidationResult object
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    # Initialize recommendations list
    recommendations = []
    
    # Check if directory exists
    if not os.path.exists(db_path):
        return DatabaseValidationResult(
            success=False,
            db_type=DatabaseType.KRAKEN,
            error_message=f"Kraken database directory does not exist: {db_path}",
            recommendations=[
                "Verify the database path is correct",
                "Download or build the Kraken database",
                "Check for typos in the path"
            ]
        )
    
    if not os.path.isdir(db_path):
        return DatabaseValidationResult(
            success=False,
            db_type=DatabaseType.KRAKEN,
            error_message=f"Kraken database path is not a directory: {db_path}",
            recommendations=[
                "Provide the path to the directory containing the Kraken database",
                "The path should be to the database directory, not a specific file"
            ]
        )
    
    # Critical files that should be present in a Kraken2 database
    required_files = [
        "hash.k2d",
        "opts.k2d",
        "taxo.k2d"
    ]
    
    # Check for each required file
    missing_files = []
    found_files = []
    
    for req_file in required_files:
        file_path = os.path.join(db_path, req_file)
        if not os.path.exists(file_path):
            missing_files.append(req_file)
        else:
            found_files.append(file_path)
    
    # Get database info
    db_info = {
        "path": db_path,
        "size": format_size(get_directory_size(db_path))
    }
    
    # Check for taxonomy files
    taxonomy_files = glob.glob(os.path.join(db_path, "taxonomy", "*.dmp"))
    if not taxonomy_files:
        if os.path.exists(os.path.join(db_path, "taxo.k2d")):
            logger.info("taxo.k2d exists but taxonomy directory not found (might be built-in)")
        else:
            missing_files.append("taxonomy/*.dmp")
            recommendations.append("The taxonomy directory may be missing or incomplete")
    else:
        found_files.extend(taxonomy_files[:3])  # Add only few taxonomy files to the found list
    
    # Check for library files
    library_files = glob.glob(os.path.join(db_path, "library", "*.fna"))
    if not library_files:
        # This is not necessarily an error, as the library may have been built into the database
        db_info["library_files"] = "Not found (possibly built-in)"
    else:
        db_info["library_files"] = f"{len(library_files)} files found"
        found_files.extend(library_files[:3])  # Add only a few library files to the found list
    
    # Check database consistency
    if missing_files:
        return DatabaseValidationResult(
            success=False,
            db_type=DatabaseType.KRAKEN,
            error_message=f"Kraken database is incomplete. Missing files: {', '.join(missing_files)}",
            recommendations=[
                "The database seems to be incomplete or corrupted",
                "Re-download or rebuild the Kraken database",
                "Check that the database was built correctly"
            ],
            files_checked=found_files + [os.path.join(db_path, f) for f in missing_files],
            db_info=db_info
        )
    
    # Calculate a hash of opts.k2d for basic integrity check
    try:
        opts_path = os.path.join(db_path, "opts.k2d")
        with open(opts_path, 'rb') as f:
            file_hash = hashlib.md5(f.read()).hexdigest()
        db_info["opts_hash"] = file_hash
    except Exception as e:
        logger.warning(f"Could not calculate hash for opts.k2d: {e}")
        recommendations.append("Warning: Could not verify database integrity")
    
    # Check for Bracken-compatible files
    bracken_files = glob.glob(os.path.join(db_path, "database*mers.kmer_distrib"))
    if not bracken_files:
        recommendations.append(
            "No Bracken kmer distribution files found. " 
            "If you plan to use Bracken, you need to build a Bracken database."
        )
    else:
        db_info["bracken_files"] = [os.path.basename(f) for f in bracken_files]
        found_files.extend(bracken_files)
    
    # Check for seqid2taxid.map if library files exist
    seqid_file = os.path.join(db_path, "seqid2taxid.map")
    if os.path.exists(seqid_file):
        found_files.append(seqid_file)
    
    # If we get this far, the database seems valid
    return DatabaseValidationResult(
        success=True,
        db_type=DatabaseType.KRAKEN,
        error_message=None,
        recommendations=recommendations,
        files_checked=found_files,
        db_info=db_info
    )

def validate_bracken_db(db_path, logger=None):
    """
    Validate a Bracken database file.
    
    Args:
        db_path: Path to Bracken kmer distribution file
        logger: Logger instance
        
    Returns:
        DatabaseValidationResult object
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    # Initialize recommendations list
    recommendations = []
    
    # Check if file exists
    if not os.path.exists(db_path):
        return DatabaseValidationResult(
            success=False,
            db_type=DatabaseType.BRACKEN,
            error_message=f"Bracken database file does not exist: {db_path}",
            recommendations=[
                "Verify the database path is correct",
                "Run bracken-build to create the kmer distribution file",
                "The file should be named something like 'database150mers.kmer_distrib'"
            ]
        )
    
    # Check if it's a file (not a directory)
    if os.path.isdir(db_path):
        # User might have provided the Kraken DB path instead of the specific file
        kmer_distrib_files = glob.glob(os.path.join(db_path, "database*mers.kmer_distrib"))
        
        if kmer_distrib_files:
            # Found potential Bracken files, suggest the first one
            return DatabaseValidationResult(
                success=False,
                db_type=DatabaseType.BRACKEN,
                error_message=(
                    f"Bracken database path is a directory, but should be a file. "
                    f"Found {len(kmer_distrib_files)} potential Bracken database files."
                ),
                recommendations=[
                    f"Use this specific file instead: {kmer_distrib_files[0]}",
                    "Bracken requires the path to the kmer distribution file, not a directory"
                ],
                files_checked=kmer_distrib_files[:5]  # Show up to 5 found files
            )
        else:
            # No Bracken files found in the directory
            return DatabaseValidationResult(
                success=False,
                db_type=DatabaseType.BRACKEN,
                error_message=(
                    f"Bracken database path is a directory, but should be a file. "
                    f"No Bracken database files found in this directory."
                ),
                recommendations=[
                    "Bracken requires the path to the kmer distribution file, not a directory",
                    "Run bracken-build to create the kmer distribution file",
                    "Check that you're using the correct database path"
                ]
            )
    
    # Check file naming pattern
    filename = os.path.basename(db_path)
    if not filename.startswith("database") or not filename.endswith(".kmer_distrib"):
        recommendations.append(
            "The file name doesn't follow the typical Bracken database naming pattern "
            "(e.g., database150mers.kmer_distrib). This might still work, but is unusual."
        )
    
    # Check for k-mer length in filename
    kmer_match = re.search(r'database(\d+)mers', filename)
    kmer_length = int(kmer_match.group(1)) if kmer_match else None
    
    # Check file size (a Bracken database should be at least a few MB)
    file_size = os.path.getsize(db_path)
    if file_size < 1024 * 1024:  # Less than 1 MB
        recommendations.append(
            f"The file size ({format_size(file_size)}) seems unusually small for a Bracken database. "
            f"This might indicate an incomplete or corrupted database."
        )
    
    # Look for the corresponding Kraken database
    kraken_db_dir = os.path.dirname(db_path)
    required_kraken_files = ["hash.k2d", "opts.k2d", "taxo.k2d"]
    missing_kraken_files = []
    
    for req_file in required_kraken_files:
        if not os.path.exists(os.path.join(kraken_db_dir, req_file)):
            missing_kraken_files.append(req_file)
    
    if missing_kraken_files:
        recommendations.append(
            "The Bracken database doesn't appear to be in a valid Kraken database directory. "
            "Bracken requires the kmer distribution file to be in the same directory as the Kraken database."
        )
    
    # Database info
    db_info = {
        "path": db_path,
        "size": format_size(file_size),
        "kraken_db_dir": kraken_db_dir
    }
    
    if kmer_length:
        db_info["kmer_length"] = kmer_length
    
    # If we get this far, the database seems valid
    return DatabaseValidationResult(
        success=len(missing_kraken_files) == 0,  # Valid only if all Kraken files present
        db_type=DatabaseType.BRACKEN,
        error_message=None if len(missing_kraken_files) == 0 else 
                     f"Bracken database found, but Kraken database is incomplete: {', '.join(missing_kraken_files)}",
        recommendations=recommendations,
        files_checked=[db_path] + [os.path.join(kraken_db_dir, f) for f in required_kraken_files],
        db_info=db_info
    )

def validate_kneaddata_db(db_path, logger=None):
    """
    Validate a KneadData reference database.
    
    Args:
        db_path: Path to KneadData database directory
        logger: Logger instance
        
    Returns:
        DatabaseValidationResult object
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    # Initialize recommendations list
    recommendations = []
    
    # Check if directory exists
    if not os.path.exists(db_path):
        return DatabaseValidationResult(
            success=False,
            db_type=DatabaseType.KNEADDATA,
            error_message=f"KneadData database does not exist: {db_path}",
            recommendations=[
                "Verify the database path is correct",
                "Download the KneadData database using kneaddata_database command",
                "Check for typos in the path"
            ]
        )
    
    # KneadData databases could be either a directory or a Bowtie2 index prefix
    if os.path.isdir(db_path):
        # Check for Bowtie2 index files in the directory
        bt2_files = glob.glob(os.path.join(db_path, "*.bt2"))
        bt2l_files = glob.glob(os.path.join(db_path, "*.bt2l"))  # Large index files
        
        if not bt2_files and not bt2l_files:
            # Look for specific files that might suggest this is a valid database
            fasta_files = glob.glob(os.path.join(db_path, "*.fasta")) + glob.glob(os.path.join(db_path, "*.fa"))
            
            if fasta_files:
                return DatabaseValidationResult(
                    success=False,
                    db_type=DatabaseType.KNEADDATA,
                    error_message=f"Found FASTA files but no Bowtie2 index files in: {db_path}",
                    recommendations=[
                        "This directory contains FASTA files but not Bowtie2 indexes",
                        "You need to build Bowtie2 indexes using 'bowtie2-build'",
                        "Or use kneaddata_database to download pre-built indexes"
                    ],
                    files_checked=fasta_files[:5],
                    db_info={"path": db_path, "fasta_files": len(fasta_files)}
                )
            else:
                return DatabaseValidationResult(
                    success=False,
                    db_type=DatabaseType.KNEADDATA,
                    error_message=f"No Bowtie2 index or FASTA files found in: {db_path}",
                    recommendations=[
                        "KneadData requires Bowtie2 index files for the reference database",
                        "Download pre-built indexes using 'kneaddata_database --download'",
                        "Check that you're providing the correct database path"
                    ],
                    files_checked=[],
                    db_info={"path": db_path}
                )
        
        # Get base name from first BT2 file
        all_bt2_files = bt2_files + bt2l_files
        if all_bt2_files:
            # Example: /path/to/db/human.1.bt2 -> human
            first_file = os.path.basename(all_bt2_files[0])
            base_name = first_file.split('.')[0]
            
            # Verify if we have a complete set of BT2 files
            expected_suffixes = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
            alternative_suffixes = ['.1.bt2l', '.2.bt2l', '.3.bt2l', '.4.bt2l', '.rev.1.bt2l', '.rev.2.bt2l']
            
            missing_files = []
            for suffix in expected_suffixes:
                if not os.path.exists(os.path.join(db_path, base_name + suffix)):
                    # Check for large index alternative
                    alt_suffix = suffix.replace('.bt2', '.bt2l')
                    if not os.path.exists(os.path.join(db_path, base_name + alt_suffix)):
                        missing_files.append(base_name + suffix)
            
            if missing_files:
                return DatabaseValidationResult(
                    success=False,
                    db_type=DatabaseType.KNEADDATA,
                    error_message=f"Incomplete Bowtie2 index files for {base_name}. Missing: {', '.join(missing_files)}",
                    recommendations=[
                        "The Bowtie2 index appears to be incomplete",
                        "Rebuild the index using 'bowtie2-build'",
                        "Or download a complete pre-built database"
                    ],
                    files_checked=all_bt2_files,
                    db_info={"path": db_path, "base_name": base_name}
                )
            
            # If we got here, the database seems valid
            db_info = {
                "path": db_path,
                "base_name": base_name,
                "type": "bt2l" if bt2l_files else "bt2",
                "size": format_size(get_directory_size(db_path))
            }
            
            return DatabaseValidationResult(
                success=True,
                db_type=DatabaseType.KNEADDATA,
                error_message=None,
                recommendations=[],
                files_checked=all_bt2_files[:6],  # Show up to 6 files (complete set)
                db_info=db_info
            )
    else:
        # The provided path is not a directory, it might be a Bowtie2 index prefix
        db_dir = os.path.dirname(db_path)
        base_name = os.path.basename(db_path)
        
        # Check if the corresponding Bowtie2 index files exist
        expected_suffixes = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
        alternative_suffixes = ['.1.bt2l', '.2.bt2l', '.3.bt2l', '.4.bt2l', '.rev.1.bt2l', '.rev.2.bt2l']
        
        found_files = []
        missing_files = []
        
        # Check for regular or large index files
        is_large_index = False
        for suffix in expected_suffixes:
            file_path = db_path + suffix
            alt_file_path = db_path + suffix.replace('.bt2', '.bt2l')
            
            if os.path.exists(file_path):
                found_files.append(file_path)
            elif os.path.exists(alt_file_path):
                found_files.append(alt_file_path)
                is_large_index = True
            else:
                missing_files.append(base_name + suffix)
        
        if missing_files:
            return DatabaseValidationResult(
                success=False,
                db_type=DatabaseType.KNEADDATA,
                error_message=f"Incomplete Bowtie2 index files for {base_name}. Missing: {', '.join(missing_files)}",
                recommendations=[
                    "The Bowtie2 index appears to be incomplete",
                    "Rebuild the index using 'bowtie2-build'",
                    "Or download a complete pre-built database"
                ],
                files_checked=found_files,
                db_info={"path": db_path, "base_name": base_name}
            )
        
        # If we got here, the database seems valid
        db_info = {
            "path": db_path,
            "base_name": base_name,
            "type": "bt2l" if is_large_index else "bt2",
            "size": format_size(sum(os.path.getsize(f) for f in found_files))
        }
        
        return DatabaseValidationResult(
            success=True,
            db_type=DatabaseType.KNEADDATA,
            error_message=None,
            recommendations=[],
            files_checked=found_files,
            db_info=db_info
        )

def print_database_validation_result(result, logger=None):
    """
    Print a detailed database validation result to logs.
    
    Args:
        result: DatabaseValidationResult object
        logger: Logger instance
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    db_type_name = result.db_type.value.capitalize() if result.db_type else "Unknown"
    
    if result.success:
        logger.info(f"✅ {db_type_name} database validation successful")
        
        # Print database info
        if result.db_info:
            logger.info("Database information:")
            for key, value in result.db_info.items():
                logger.info(f"  - {key}: {value}")
    else:
        logger.error(f"❌ {db_type_name} database validation failed: {result.error_message}")
    
    # Print recommendations
    if result.recommendations:
        logger.info("Recommendations:")
        for i, rec in enumerate(result.recommendations, 1):
            logger.info(f"  {i}. {rec}")
    
    # Print files checked (debug level)
    if result.files_checked:
        logger.debug("Files checked:")
        for f in result.files_checked:
            logger.debug(f"  - {f}")

def validate_database(db_path, db_type, logger=None):
    """
    Validate a database of specified type.
    
    Args:
        db_path: Path to database
        db_type: DatabaseType enum
        logger: Logger instance
        
    Returns:
        DatabaseValidationResult object
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    logger.info(f"Validating {db_type.value} database at: {db_path}")
    
    if db_type == DatabaseType.KRAKEN:
        result = validate_kraken_db(db_path, logger)
    elif db_type == DatabaseType.BRACKEN:
        result = validate_bracken_db(db_path, logger)
    elif db_type == DatabaseType.KNEADDATA:
        result = validate_kneaddata_db(db_path, logger)
    else:
        result = DatabaseValidationResult(
            success=False,
            db_type=None,
            error_message=f"Unknown database type: {db_type}"
        )
    
    print_database_validation_result(result, logger)
    return result