# kraken_tools/utils/cmd_utils.py
import os
import sys
import subprocess
import logging
import re
import time
from functools import wraps
from enum import Enum

class CommandErrorType(Enum):
    """Categorize different types of command errors for better handling."""
    SUCCESS = 0
    NOT_FOUND = 1           # Command not found (installation issue)
    PERMISSION_ERROR = 2    # Permission denied
    DATABASE_ERROR = 3      # Database issues
    MEMORY_ERROR = 4        # Out of memory
    DISK_ERROR = 5          # Out of disk space
    NETWORK_ERROR = 6       # Network connectivity issues
    TIMEOUT_ERROR = 7       # Command timed out
    INPUT_ERROR = 8         # Issues with input files
    OUTPUT_ERROR = 9        # Issues with output files/directories
    CONFIG_ERROR = 10       # Configuration problems
    UNKNOWN_ERROR = 99      # Unclassified error

class CommandResult:
    """Object to store detailed information about command execution results."""
    def __init__(self, success=False, error_type=None, returncode=None, 
                 stdout=None, stderr=None, error_message=None, cmd=None):
        self.success = success
        self.error_type = error_type or CommandErrorType.UNKNOWN_ERROR
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        self.error_message = error_message
        self.cmd = cmd
        
    def __bool__(self):
        return self.success

def classify_error(stderr, returncode):
    """
    Classify the error based on stderr content and return code.
    
    Args:
        stderr: Standard error output as string
        returncode: Process return code
        
    Returns:
        Tuple of (CommandErrorType, error_message)
    """
    # Convert to lowercase for case-insensitive matching
    stderr_lower = stderr.lower() if stderr else ""
    
    # Check for specific error patterns
    if returncode == 127 or "command not found" in stderr_lower:
        return CommandErrorType.NOT_FOUND, "Command not found (check installation)"
    
    if "permission denied" in stderr_lower:
        return CommandErrorType.PERMISSION_ERROR, "Permission denied (check file permissions)"
    
    if any(x in stderr_lower for x in ["database", "db", "reference"]) and any(
            x in stderr_lower for x in ["not found", "cannot open", "does not exist", "failed to load"]):
        return CommandErrorType.DATABASE_ERROR, "Database error (check database path and structure)"
    
    if any(x in stderr_lower for x in ["memory", "out of memory", "cannot allocate"]):
        return CommandErrorType.MEMORY_ERROR, "Out of memory (try using less threads or a smaller dataset)"
    
    if any(x in stderr_lower for x in ["disk", "no space", "cannot write"]):
        return CommandErrorType.DISK_ERROR, "Disk error (check available disk space)"
    
    if any(x in stderr_lower for x in ["network", "connection", "timeout", "unreachable"]):
        return CommandErrorType.NETWORK_ERROR, "Network error (check connectivity)"
    
    if "timeout" in stderr_lower:
        return CommandErrorType.TIMEOUT_ERROR, "Command timed out (try increasing timeout)"
    
    if any(x in stderr_lower for x in ["input file", "no such file", "invalid", "corrupt"]):
        return CommandErrorType.INPUT_ERROR, "Input file error (check input files)"
    
    if any(x in stderr_lower for x in ["output", "cannot create", "directory"]):
        return CommandErrorType.OUTPUT_ERROR, "Output error (check output directory permissions)"
    
    if any(x in stderr_lower for x in ["config", "configuration", "parameter", "option", "invalid"]):
        return CommandErrorType.CONFIG_ERROR, "Configuration error (check command options)"
    
    # If no specific match, return unknown error
    return CommandErrorType.UNKNOWN_ERROR, f"Unknown error (return code: {returncode})"

def with_retries(max_attempts=3, initial_delay=1, backoff_factor=2, retryable_errors=None):
    """
    Decorator for functions that should be retried on failure.
    
    Args:
        max_attempts: Maximum number of retry attempts
        initial_delay: Initial delay between retries in seconds
        backoff_factor: Factor by which the delay increases for each retry
        retryable_errors: List of CommandErrorType that should trigger retries
        
    Returns:
        Decorated function
    """
    if retryable_errors is None:
        # By default, retry network, timeout and memory errors
        retryable_errors = [
            CommandErrorType.NETWORK_ERROR,
            CommandErrorType.TIMEOUT_ERROR,
            CommandErrorType.MEMORY_ERROR
        ]
    
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            logger = logging.getLogger('kraken_analysis')
            attempts = 0
            delay = initial_delay
            
            while attempts < max_attempts:
                result = func(*args, **kwargs)
                
                # If successful or not a retryable error, return immediately
                if result.success or result.error_type not in retryable_errors:
                    return result
                
                attempts += 1
                if attempts < max_attempts:
                    logger.warning(
                        f"Retryable error: {result.error_message}. "
                        f"Retrying in {delay} seconds (attempt {attempts}/{max_attempts-1})..."
                    )
                    time.sleep(delay)
                    delay *= backoff_factor
                else:
                    logger.error(f"Failed after {max_attempts} attempts: {result.error_message}")
            
            return result
        
        return wrapper
    
    return decorator

@with_retries(max_attempts=3, initial_delay=1, backoff_factor=2)
def run_cmd_with_retry(cmd, timeout=None, **kwargs):
    """
    Run a command with automatic retry mechanism.
    
    Args:
        cmd: Command to run as a list of strings
        timeout: Command timeout in seconds
        **kwargs: Additional args passed to run_cmd
        
    Returns:
        CommandResult object
    """
    return run_cmd(cmd, timeout=timeout, **kwargs)

def run_cmd(cmd, exit_on_error=True, verbose=True, timeout=None):
    """
    Run a shell command with detailed error handling.
    
    Args:
        cmd: Command to run as a list of strings
        exit_on_error: Whether to exit the program if the command fails
        verbose: Whether to log the command being run
        timeout: Command timeout in seconds
        
    Returns:
        CommandResult object containing detailed execution information
    """
    logger = logging.getLogger('kraken_analysis')
    
    if verbose:
        logger.info(f"Running: {' '.join(cmd)}")
    
    # Check for commands that need special handling
    if cmd[0] == "cp" and len(cmd) >= 3:
        src = cmd[1]
        if not os.path.exists(src):
            error_msg = f"Source file does not exist: {src}"
            logger.error(f"ERROR: {error_msg}")
            
            result = CommandResult(
                success=False,
                error_type=CommandErrorType.INPUT_ERROR,
                returncode=1,
                stderr=error_msg,
                error_message=error_msg,
                cmd=cmd
            )
            
            if exit_on_error:
                sys.exit(1)
            return result
        
        dst_dir = os.path.dirname(cmd[2])
        if not os.path.exists(dst_dir):
            logger.info(f"Creating directory: {dst_dir}")
            os.makedirs(dst_dir, exist_ok=True)
    
    try:
        process = subprocess.run(
            cmd, 
            check=False,  # Don't raise exception, we'll handle errors ourselves
            stderr=subprocess.PIPE, 
            stdout=subprocess.PIPE,
            timeout=timeout
        )
        
        stdout = process.stdout.decode('utf-8', errors='replace')
        stderr = process.stderr.decode('utf-8', errors='replace')
        
        if stdout.strip() and verbose:
            logger.debug(f"Command stdout: {stdout}")
        
        # Check if command was successful
        if process.returncode == 0:
            if stderr.strip() and verbose:
                logger.debug(f"Command stderr (non-fatal): {stderr}")
            
            return CommandResult(
                success=True,
                error_type=CommandErrorType.SUCCESS,
                returncode=0,
                stdout=stdout,
                stderr=stderr,
                cmd=cmd
            )
        else:
            # Command failed, classify the error
            error_type, error_msg = classify_error(stderr, process.returncode)
            
            # Create a more detailed error message
            detailed_error = f"Command failed: {error_msg}"
            if stderr.strip():
                # Try to extract the most relevant part of stderr
                # Often, the last few lines contain the most useful info
                stderr_lines = stderr.strip().split('\n')
                if len(stderr_lines) > 5:
                    stderr_extract = '\n'.join(stderr_lines[-5:])
                    detailed_error += f"\nError details (last 5 lines):\n{stderr_extract}"
                else:
                    detailed_error += f"\nError details:\n{stderr}"
            
            logger.error(detailed_error)
            logger.error(f"Failed command: {' '.join(cmd)}")
            
            result = CommandResult(
                success=False,
                error_type=error_type,
                returncode=process.returncode,
                stdout=stdout,
                stderr=stderr,
                error_message=error_msg,
                cmd=cmd
            )
            
            if exit_on_error:
                sys.exit(1)
            
            return result
            
    except subprocess.TimeoutExpired:
        error_msg = f"Command timed out after {timeout} seconds"
        logger.error(f"ERROR: {error_msg}")
        logger.error(f"Timed out command: {' '.join(cmd)}")
        
        result = CommandResult(
            success=False,
            error_type=CommandErrorType.TIMEOUT_ERROR,
            returncode=None,
            stderr=error_msg,
            error_message=error_msg,
            cmd=cmd
        )
        
        if exit_on_error:
            sys.exit(1)
        
        return result
        
    except Exception as e:
        error_msg = f"Error executing command: {str(e)}"
        logger.error(f"ERROR: {error_msg}")
        logger.error(f"Failed command: {' '.join(cmd)}")
        
        result = CommandResult(
            success=False,
            error_type=CommandErrorType.UNKNOWN_ERROR,
            returncode=None,
            stderr=str(e),
            error_message=error_msg,
            cmd=cmd
        )
        
        if exit_on_error:
            sys.exit(1)
        
        return result