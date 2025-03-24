#!/bin/bash
#
# conda_setup.sh - Set up Conda environment for Kraken Tools
#
# This script creates a conda environment with all necessary dependencies
# for running the Kraken Tools package on Unix/Linux/macOS systems.
#
# Usage:
#   ./conda_setup.sh [--name ENV_NAME] [--python PYTHON_VERSION] 
#                    [--bioconda-channel CHANNEL] [--force]
#
# Examples:
#   # Create default environment (kraken-tools)
#   ./conda_setup.sh
#   
#   # Create environment with specific name and Python version
#   ./conda_setup.sh --name my-kraken-env --python 3.12
#   
#   # Force recreation of an existing environment
#   ./conda_setup.sh --force

set -e

# Default values
ENV_NAME="kraken-tools"
PYTHON_VERSION="3.12"
BIOCONDA_CHANNEL="bioconda"
FORCE=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --name)
      ENV_NAME="$2"
      shift 2
      ;;
    --python)
      PYTHON_VERSION="$2"
      shift 2
      ;;
    --bioconda-channel)
      BIOCONDA_CHANNEL="$2"
      shift 2
      ;;
    --force)
      FORCE=true
      shift
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: ./conda_setup.sh [--name ENV_NAME] [--python PYTHON_VERSION] [--bioconda-channel CHANNEL] [--force]"
      exit 1
      ;;
  esac
done

# Print banner
echo "================================================================================"
echo "Kraken Tools Conda Environment Setup"
echo "================================================================================"
echo "Environment name: $ENV_NAME"
echo "Python version: $PYTHON_VERSION"
echo "Bioconda channel: $BIOCONDA_CHANNEL"
echo "Force recreation: $FORCE"
echo "================================================================================"

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found in PATH. Please install conda first."
    exit 1
fi

echo "Using conda installation at: $(which conda)"

# Create or use existing environment
if conda env list | grep -q "^$ENV_NAME "; then
    if [ "$FORCE" = true ]; then
        echo "Environment '$ENV_NAME' already exists. Removing it as requested..."
        conda env remove --name $ENV_NAME
        conda create --name $ENV_NAME python=$PYTHON_VERSION -y
    else
        echo "Environment '$ENV_NAME' already exists. Use --force to recreate it."
        echo "Proceeding with dependency installation..."
    fi
else
    echo "Creating conda environment: $ENV_NAME with Python $PYTHON_VERSION"
    conda create --name $ENV_NAME python=$PYTHON_VERSION -y
fi

# Install core dependencies from conda-forge
echo "Installing core dependencies..."
conda install --name $ENV_NAME -y -c conda-forge \
    pandas \
    numpy \
    matplotlib \
    seaborn \
    scipy \
    scikit-learn \
    scikit-posthocs \
    statsmodels \
    tqdm \
    psutil \
    matplotlib-venn || {
        echo "Warning: Some dependencies failed to install together. Trying individually..."
        for pkg in pandas numpy matplotlib seaborn scipy scikit-learn scikit-posthocs statsmodels tqdm psutil matplotlib-venn; do
            conda install --name $ENV_NAME -y -c conda-forge $pkg || echo "Warning: Failed to install $pkg"
        done
    }

# Install bioconda tools
echo "Installing bioconda tools from $BIOCONDA_CHANNEL channel..."
conda install --name $ENV_NAME -y -c $BIOCONDA_CHANNEL \
    kraken2 \
    bracken \
    kneaddata || {
        echo "Warning: Some bioconda tools failed to install together. Trying individually..."
        for pkg in kraken2 bracken kneaddata; do
            conda install --name $ENV_NAME -y -c $BIOCONDA_CHANNEL $pkg || echo "Warning: Failed to install $pkg"
        done
    }

# Install pip dependencies
echo "Installing additional dependencies with pip..."
conda run --name $ENV_NAME pip install \
    scikit-bio || {
        echo "Warning: Some pip dependencies failed to install together. Trying individually..."
        conda run --name $ENV_NAME pip install scikit-bio || echo "Warning: Failed to install scikit-bio"
    }

# Install current package in development mode
echo "Installing kraken-tools package..."
conda run --name $ENV_NAME pip install -e . || {
    echo "Warning: Could not install kraken-tools in development mode."
    echo "Please install it manually after activation:"
    echo "    conda activate $ENV_NAME"
    echo "    pip install -e ."
}

# Create activation script
cat > activate_env.sh << EOF
#!/bin/bash
# Activate the $ENV_NAME conda environment
conda activate $ENV_NAME
EOF

chmod +x activate_env.sh
echo "Created activation script: activate_env.sh"

# Print success message and next steps
echo
echo "================================================================================"
echo "Kraken Tools conda environment '$ENV_NAME' has been set up successfully!"
echo "================================================================================"
echo
echo "To activate the environment, run:"
echo "    conda activate $ENV_NAME"
echo
echo "To verify the installation, run:"
echo "    kraken-tools --help"
echo
echo "To get started with a full pipeline workflow, run:"
echo "    kraken-tools full-pipeline \\"
echo "        --input-fastq reads_1.fastq.gz reads_2.fastq.gz \\"
echo "        --paired \\"
echo "        --kneaddata-dbs /path/to/kneaddata_db \\"
echo "        --kraken-db /path/to/kraken_db \\"
echo "        --bracken-db /path/to/kraken_db/database150mers.kmer_distrib \\"
echo "        --sample-key metadata.csv \\"
echo "        --output-dir results/ \\"
echo "        --group-col 'Group' \\"
echo "        --threads 8"
echo
echo "Or for differential abundance analysis:"
echo "    kraken-tools diff-abundance \\"
echo "        --abundance-file abundance.tsv \\"
echo "        --sample-key metadata.csv \\"
echo "        --output-dir results/diff_abundance/ \\"
echo "        --group-col 'Group' \\"
echo "        --methods aldex2,ancom,ancom-bc"
echo
echo "Refer to the README.md for more examples and options."
echo "================================================================================"
