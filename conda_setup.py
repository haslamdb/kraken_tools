#!/usr/bin/env python3
"""
conda_setup.py - Set up Conda environment for Kraken Tools

This script creates a conda environment with all necessary dependencies
for running the Kraken Tools package. It handles both the core
taxonomic analysis tools (Kraken2, Bracken, KneadData) and the additional
analysis dependencies.

Usage:
    python conda_setup.py [--name ENV_NAME] [--python PYTHON_VERSION]
                         [--bioconda-channel CHANNEL] [--force]

Examples:
    # Create default environment (kraken-tools)
    python conda_setup.py
    
    # Create environment with specific name and Python version
    python conda_setup.py --name my-kraken-env --python 3.12
    
    # Force recreation of an existing environment
    python conda_setup.py --force
"""

import os
import sys
import subprocess
import argparse
import platform
import shutil
from pathlib import Path


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Set up Conda environment for Kraken Tools"
    )
    
    parser.add_argument(
        "--name", default="kraken-tools",
        help="Name of the conda environment to create (default: kraken-tools)"
    )
    parser.add_argument(
        "--python", default="3.12",
        help="Python version to use (default: 3.12)"
    )
    parser.add_argument(
        "--bioconda-channel", default="bioconda",
        help="Conda channel for bioconda tools (default: bioconda)"
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Force recreation of the environment if it already exists"
    )
    
    return parser.parse_args()


def check_conda():
    """Check if conda is installed and available."""
    # On Windows, try both conda and conda.exe
    if platform.system() == "Windows":
        conda_path = shutil.which("conda.exe") or shutil.which("conda") or shutil.which("conda.bat")
    else:
        conda_path = shutil.which("conda")
    
    if not conda_path:
        print("ERROR: conda not found in PATH. Please install conda first.")
        sys.exit(1)
    
    print(f"Found conda at: {conda_path}")
    return conda_path


def get_os_info():
    """Get information about the operating system."""
    system = platform.system()
    if system == "Linux":
        return "linux"
    elif system == "Darwin":
        return "macos"
    elif system == "Windows":
        return "windows"
    else:
        print(f"WARNING: Unsupported OS: {system}. Some features may not work.")
        return "unknown"


def create_environment(env_name, python_version, force=False):
    """Create a new conda environment with the specified Python version."""
    print(f"Creating conda environment: {env_name} with Python {python_version}")
    
    # Get conda executable path
    conda_path = check_conda()
    
    # Check if environment already exists
    try:
        # First try using conda env list
        result = subprocess.run(
            [conda_path, "env", "list", "--json"],
            capture_output=True, text=True, check=True
        )
        
        import json
        envs = json.loads(result.stdout)["envs"]
        env_exists = any(env_name in env for env in envs)
    except (subprocess.SubprocessError, json.JSONDecodeError):
        # Fallback if json method fails
        print("Warning: Could not get environment list in JSON format. Using text output.")
        result = subprocess.run(
            [conda_path, "env", "list"],
            capture_output=True, text=True
        )
        env_exists = env_name in result.stdout
    
    if env_exists:
        if force:
            print(f"Environment '{env_name}' already exists. Removing it as requested...")
            subprocess.run([conda_path, "env", "remove", "--name", env_name, "-y"], check=True)
        else:
            print(f"Environment '{env_name}' already exists. Use --force to recreate it.")
            return False
    
    # Create the environment
    print(f"Running: {conda_path} create --name {env_name} python={python_version} -y")
    subprocess.run(
        [conda_path, "create", "--name", env_name, f"python={python_version}", "-y"],
        check=True
    )
    
    return True


def install_dependencies(env_name, bioconda_channel):
    """Install all dependencies in the conda environment."""
    print("Installing dependencies...")
    
    # Get conda executable path
    conda_path = check_conda()
    
    # Define the commands to run with conda and pip
    conda_cmd = [conda_path, "install", "--name", env_name, "-y"]
    
    # Core dependencies available from conda main channels
    core_deps = [
        "pandas",
        "numpy",
        "matplotlib",
        "seaborn",
        "scipy",
        "scikit-learn",
        "statsmodels",
        "tqdm",
        "psutil"
    ]
    
    # Bioconda tools
    bioconda_deps = [
        "kraken2",
        "bracken",
        "kneaddata"
    ]
    
    # Dependencies that are better installed with pip
    pip_deps = [
        "scikit-posthocs",   # For Dunn's test
        "scikit-bio",        
        "matplotlib-venn",
        "shap"   
    ]
    
    # Install core dependencies from conda
    print("Installing core dependencies...")
    try:
        subprocess.run(conda_cmd + ["-c", "conda-forge"] + core_deps, check=True)
    except subprocess.CalledProcessError:
        print("Warning: Some core dependencies could not be installed. Trying individually...")
        for dep in core_deps:
            try:
                subprocess.run(conda_cmd + ["-c", "conda-forge", dep], check=False)
            except subprocess.CalledProcessError:
                print(f"  Could not install {dep} from conda. Will try pip later.")
    
    # Install bioconda tools
    print(f"Installing bioconda tools from {bioconda_channel} channel...")
    try:
        subprocess.run(conda_cmd + ["-c", bioconda_channel] + bioconda_deps, check=True)
    except subprocess.CalledProcessError:
        print("Warning: Some bioconda tools could not be installed. Trying individually...")
        for dep in bioconda_deps:
            try:
                subprocess.run(conda_cmd + ["-c", bioconda_channel, dep], check=False)
            except subprocess.CalledProcessError:
                print(f"  Could not install {dep}. Please install it manually after setup.")
    
    # Install pip dependencies
    print("Installing additional dependencies with pip...")
    pip_cmd = [
        conda_path, "run", "--name", env_name,
        "python", "-m", "pip", "install"
    ] + pip_deps
    
    # Try installing all pip dependencies at once, fall back to individual installation if needed
    try:
        subprocess.run(pip_cmd, check=True)
    except subprocess.CalledProcessError:
        print("Warning: Could not install all pip dependencies together. Trying individually...")
        for dep in pip_deps:
            try:
                subprocess.run([conda_path, "run", "--name", env_name, 
                              "python", "-m", "pip", "install", dep], check=False)
            except subprocess.CalledProcessError:
                print(f"  Could not install {dep} with pip.")
    
    # Install current package in development mode
    print("Installing kraken-tools package...")
    try:
        subprocess.run([
            conda_path, "run", "--name", env_name,
            "python", "-m", "pip", "install", "-e", "."
        ], check=True)
    except subprocess.CalledProcessError:
        print("Warning: Could not install kraken-tools in development mode. Please install it manually.")
        print("    cd kraken-tools")
        print("    conda activate", env_name)
        print("    pip install -e .")


def setup_post_installation(env_name):
    """Perform post-installation setup tasks."""
    print("Setting up post-installation configurations...")
    
    # Create a convenience script for activating the environment
    script_name = "activate_env.sh" if get_os_info() != "windows" else "activate_env.bat"
    
    try:
        with open(script_name, "w") as f:
            if get_os_info() != "windows":
                f.write(f"#!/bin/bash\n")
                f.write(f"# Activate the {env_name} conda environment\n")
                f.write(f"conda activate {env_name}\n")
            else:
                f.write(f"@echo off\n")
                f.write(f":: Activate the {env_name} conda environment\n")
                # On Windows, we need to handle Anaconda vs Miniconda differences
                f.write(f"call conda activate {env_name}\n")
                f.write(f"if errorlevel 1 (\n")
                # Fallback for older Anaconda on Windows or when activate script has different path
                f.write(f"    echo Trying alternative activation method...\n")
                f.write(f"    call activate {env_name}\n")
                f.write(f")\n")
        
        if get_os_info() != "windows":
            os.chmod(script_name, 0o755)
        
        print(f"Created activation script: {script_name}")
    except Exception as e:
        print(f"Warning: Could not create activation script: {e}")
        print(f"You can still activate the environment with: conda activate {env_name}")


def print_next_steps(env_name):
    """Print information about next steps after installation."""
    print("\n" + "="*80)
    print(f"Kraken Tools conda environment '{env_name}' has been set up successfully!")
    print("="*80)
    print("\nTo activate the environment, run:")
    print(f"    conda activate {env_name}")
    print("\nTo verify the installation, run:")
    print("    kraken-tools --help")
    print("\nTo get started with a full pipeline workflow, run:")
    print("    kraken-tools full-pipeline \\")
    print("        --input-fastq reads_1.fastq.gz reads_2.fastq.gz \\")
    print("        --paired \\")
    print("        --kneaddata-dbs /path/to/kneaddata_db \\")
    print("        --kraken-db /path/to/kraken_db \\")
    print("        --bracken-db /path/to/kraken_db/database150mers.kmer_distrib \\")
    print("        --sample-key metadata.csv \\")
    print("        --output-dir results/ \\")
    print("        --group-col 'Group' \\")
    print("        --threads 8")
    print("\nRefer to the README.md for more examples and options.")
    print("="*80)


def main():
    """Main function to set up the conda environment."""
    args = parse_args()
    
    # We moved the check_conda() to the individual functions where conda is used
    
    try:
        # Create the environment
        if not create_environment(args.name, args.python, args.force):
            print("\nUsing existing environment. Proceeding with dependency installation...")
        
        # Install dependencies
        install_dependencies(args.name, args.bioconda_channel)
        
        # Perform post-installation setup
        setup_post_installation(args.name)
        
        # Print next steps
        print_next_steps(args.name)
    except Exception as e:
        print(f"\nERROR: Setup failed: {e}")
        print("\nTroubleshooting tips:")
        print("1. Make sure conda is properly installed and in your PATH")
        print("2. Try running with --force to recreate the environment")
        print("3. Try specifying an older Python version: --python 3.8")
        print("4. If on Windows, try running the script from Anaconda Prompt")
        print("5. Check that you have administrative privileges if needed")
        sys.exit(1)


if __name__ == "__main__":
    main()
