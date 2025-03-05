# kraken_tools/__init__.py
"""
Kraken Tools - A package for processing and analyzing Kraken2/Bracken taxonomic output data.

This package provides functions for:
1. Running Kraken2 for taxonomic classification
2. Running Bracken for abundance estimation
3. Processing Kraken2/Bracken output files
4. Downstream analysis (statistical tests, visualization)
5. Differential abundance analysis
"""

__version__ = "0.1.0"

# Import main functions for easy access
from kraken_tools.main import (
    run_full_pipeline,
    process_kraken_files_only,
    analyze_existing_kraken_files,
    run_taxonomic_differential_abundance,
)

# Make logger functions available at the package level
from kraken_tools.logger import setup_logger, log_print

# Import key utility functions
from kraken_tools.utils.file_utils import check_file_exists, sanitize_filename

# Import differential abundance functions
from kraken_tools.analysis.differential import (
    aldex2_like,
    ancom,
    ancom_bc,
    run_differential_abundance_analysis
)

# Import preprocessing functions for easier access
from kraken_tools.preprocessing.kneaddata import (
    run_kneaddata,
    run_kneaddata_parallel
)

from kraken_tools.preprocessing.kraken_run import (
    run_kraken,
    run_kraken_parallel
)

from kraken_tools.preprocessing.bracken_run import (
    run_bracken,
    run_bracken_parallel
)
