# Kraken Tools - Installation Guide

This guide provides detailed instructions for installing the Kraken Tools package using various methods. Choose the installation method that works best for your system and requirements.

## Prerequisites

Before installing Kraken Tools, ensure you have the following prerequisites:

### For All Methods:
- Git (to clone the repository)
- Internet connection (to download dependencies)

### For Conda Installation:
- Conda (Anaconda or Miniconda)
- Python 3.9 or newer (3.12 recommended)

### For Pip Installation:
- Python 3.9 or newer (3.12 recommended)
- pip (Python package installer)

## Installation Methods

Kraken Tools offers three primary installation methods:

1. **Conda Setup Script (Python)**: Recommended for Windows users or those who prefer Python scripts
2. **Conda Setup Shell Script**: Recommended for Unix/Linux/macOS users
3. **Pip Installation**: For users who prefer to manage dependencies manually or don't use conda

## 1. Conda Setup Script (Python)

This method uses a Python script to create a conda environment with all required dependencies.

```bash
# Clone the repository
git clone https://github.com/haslamdb/kraken_tools.git

# Navigate to the package directory
cd kraken_tools

# Run the setup script
python conda_setup.py
```

### Options

The conda setup script supports several command-line options:

```bash
python conda_setup.py [--name ENV_NAME] [--python PYTHON_VERSION] [--bioconda-channel CHANNEL] [--force]
```

- `--name`: Name of the conda environment (default: kraken-tools)
- `--python`: Python version to use (default: 3.12)
- `--bioconda-channel`: Conda channel for bioconda tools (default: bioconda)
- `--force`: Force recreation of an existing environment

### Examples

```bash
# Create a custom environment with Python 3.10
python conda_setup.py --name my-kraken-env --python 3.10

# Recreate an existing environment
python conda_setup.py --force
```

## 2. Conda Setup Shell Script

For Unix/Linux/macOS users, a shell script is provided for convenience.

```bash
# Clone the repository
git clone https://github.com/haslamdb/kraken_tools.git

# Navigate to the package directory
cd kraken_tools

# Make the script executable
chmod +x conda_setup.sh

# Run the setup script
./conda_setup.sh
```

### Options

The shell script supports the same options as the Python script:

```bash
./conda_setup.sh [--name ENV_NAME] [--python PYTHON_VERSION] [--bioconda-channel CHANNEL] [--force]
```

### Examples

```bash
# Create a custom environment with Python 3.10
./conda_setup.sh --name my-kraken-env --python 3.10

# Recreate an existing environment
./conda_setup.sh --force
```

## 3. Manual Conda Installation

If you prefer to set up the conda environment step by step:

```bash
# Clone the repository
git clone https://github.com/haslamdb/kraken_tools.git
cd kraken_tools

# Create a conda environment
conda create -n kraken-tools python=3.12 -y
conda activate kraken-tools

# Install core dependencies from conda-forge
conda install -c conda-forge pandas numpy scipy scikit-bio scikit-learn scikit-posthocs statsmodels matplotlib seaborn matplotlib-venn tqdm psutil -y

# Install bioconda tools
conda install -c bioconda kraken2 bracken kneaddata -y

# Install the package
pip install -e .
```

## 4. Pip Installation

For users who prefer pip or want to manage dependencies manually:

```bash
# Clone the repository
git clone https://github.com/haslamdb/kraken_tools.git

# Navigate to the package directory
cd kraken_tools

# Install the package and dependencies
pip install -e .

# Or install with development dependencies
pip install -e ".[dev]"
```

### Note on External Tools

When using pip installation, you'll need to manually install Kraken2, Bracken, and KneadData. These tools are typically installed via conda or from their respective repositories:

```bash
# Install bioconda tools with conda
conda install -c bioconda kraken2 bracken kneaddata
```

Or follow the installation instructions for each tool:
- [Kraken2 Installation](https://github.com/DerrickWood/kraken2/wiki/Manual)
- [Bracken Installation](https://github.com/jenniferlu717/Bracken)
- [KneadData Installation](https://github.com/biobakery/kneaddata)

## Verifying the Installation

After installation, verify that Kraken Tools is correctly installed:

```bash
# Activate the conda environment (if using conda)
conda activate kraken-tools  # or your custom environment name

# Test the command-line tool
kraken-tools --help
```

You should see the help message for the `kraken-tools` command.

## Troubleshooting

### Common Installation Issues

#### Conda Environment Creation Fails

If you encounter errors during conda environment creation:

1. Check that conda is properly installed and in your PATH
2. Try specifying a different Python version: `--python 3.9`
3. On Windows, run the script from Anaconda Prompt
4. Check that you have administrative privileges if needed

#### Dependency Installation Failures

If some dependencies fail to install:

1. Make sure you have internet connectivity
2. Try updating conda: `conda update -n base conda`
3. Try installing the problematic package manually:
   ```bash
   conda install -n kraken-tools -c conda-forge package_name
   ```

#### Bioconda Tools Installation Issues

If Kraken2, Bracken, or KneadData fail to install:

1. Check the bioconda channel: `conda config --show channels`
2. Add the bioconda channel: `conda config --add channels bioconda`
3. Try installing each tool individually:
   ```bash
   conda install -n kraken-tools -c bioconda kraken2
   conda install -n kraken-tools -c bioconda bracken
   conda install -n kraken-tools -c bioconda kneaddata
   ```

### Getting Help

If you encounter issues not covered in this guide:

1. Check the [GitHub repository issues](https://github.com/haslamdb/kraken_tools/issues)
2. Run commands with verbose logging: `--log-level DEBUG`
3. Open a new issue on GitHub with details about your environment and the error messages

## Post-Installation Setup

After installing Kraken Tools, you'll need reference databases for various steps:

### 1. Kraken2 Database

Kraken2 requires a database for taxonomic classification. You can:
- Download a pre-built database from the [Kraken2 website](https://github.com/DerrickWood/kraken2/wiki/Manual#standard-kraken-2-database)
- Build your own custom database using the `kraken2-build` command

```bash
# Example of downloading a pre-built database
mkdir -p /path/to/kraken_db
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230605.tar.gz
tar -xzf k2_standard_20230605.tar.gz -C /path/to/kraken_db
```

### 2. Bracken Database

Bracken uses the Kraken2 database with additional files for abundance estimation:

```bash
# Generate Bracken database with read length of 150bp
bracken-build -d /path/to/kraken_db -k 35 -l 150 -t 8
```

### 3. KneadData Reference Database

KneadData requires reference databases for host sequence removal:

```bash
# Download human genome database
kneaddata_database --download human_genome bowtie2 /path/to/kneaddata_db
```

## Environment Variables (Optional)

You can set environment variables to avoid typing database paths repeatedly:

```bash
# Add to your .bashrc or .bash_profile
export KRAKEN_DB_PATH="/path/to/kraken_db"
export BRACKEN_DB_PATH="/path/to/kraken_db/database150mers.kmer_distrib"
export KNEADDATA_DB_PATH="/path/to/kneaddata_db"
```

Then you can use them in your commands:

```bash
kraken-tools full-pipeline \
    --input-fastq reads_1.fastq.gz reads_2.fastq.gz \
    --paired \
    --kneaddata-dbs $KNEADDATA_DB_PATH \
    --kraken-db $KRAKEN_DB_PATH \
    --bracken-db $BRACKEN_DB_PATH \
    --sample-key metadata.csv \
    --output-dir results/
```
