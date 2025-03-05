# kraken_tools/utils/cmd_utils.py
import os
import sys
import subprocess
import logging

def run_cmd(cmd, exit_on_error=True, verbose=True):
    """
    Utility function to run a shell command with subprocess.
    
    Args:
        cmd (list): Command to run as a list of strings
        exit_on_error (bool): Whether to exit the program if the command fails
        verbose (bool): Whether to print/log the command being run
        
    Returns:
        bool: True if command executed successfully, False otherwise
    """
    logger = logging.getLogger('kraken_analysis')
    
    if verbose:
        logger.info(f"Running: {' '.join(cmd)}")
    
    # Additional check for 'cp' commands
    if cmd[0] == "cp" and len(cmd) >= 3:
        src = cmd[1]
        if not os.path.exists(src):
            logger.error(f"ERROR: Source file does not exist: {src}")
            if exit_on_error:
                sys.exit(1)
            return False
        
        dst_dir = os.path.dirname(cmd[2])
        if not os.path.exists(dst_dir):
            logger.info(f"Creating directory: {dst_dir}")
            os.makedirs(dst_dir, exist_ok=True)
    
    try:
        process = subprocess.run(cmd, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stdout = process.stdout.decode('utf-8', errors='replace')
        stderr = process.stderr.decode('utf-8', errors='replace')
        if stdout.strip():
            logger.debug(f"Command stdout: {stdout}")
        if stderr.strip():
            logger.debug(f"Command stderr: {stderr}")
        return True
    except subprocess.CalledProcessError as e:
        error_msg = e.stderr.decode('utf-8', errors='replace')
        logger.error(f"ERROR: Command failed with exit code {e.returncode}")
        logger.error(f"Error message: {error_msg}")
        logger.error(f"Failed command: {' '.join(cmd)}")
        if exit_on_error:
            sys.exit(1)
        return False