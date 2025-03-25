# kraken_tools/utils/resource_utils.py

import logging
import os
import multiprocessing
import psutil
import resource
import threading
import time
from functools import wraps

def get_memory_usage():
    """
    Get current memory usage in MB for the current process.
    """
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB

def monitor_memory_usage(logger, threshold_mb=1000, interval=60):
    """
    Start a background thread to monitor memory usage and log warnings 
    if it exceeds the threshold.
    
    Args:
        logger: Logger instance.
        threshold_mb: Memory threshold in MB to trigger warnings.
        interval: Check interval in seconds.
        
    Returns:
        A tuple (monitor_thread, stop_event) that can be used to stop monitoring.
    """
    stop_event = threading.Event()

    def monitoring_worker():
        while not stop_event.is_set():
            mem_usage = get_memory_usage()
            if mem_usage > threshold_mb:
                logger.warning(
                    f"High memory usage detected: {mem_usage:.2f} MB "
                    f"(threshold: {threshold_mb} MB)"
                )
            time.sleep(interval)

    monitor_thread = threading.Thread(target=monitoring_worker, daemon=True)
    monitor_thread.start()
    
    return monitor_thread, stop_event

def stop_memory_monitoring(thread_data):
    """
    Stop the memory monitoring thread.
    
    Args:
        thread_data: A tuple (thread, stop_event) as returned by monitor_memory_usage.
    """
    if thread_data:
        thread, stop_event = thread_data
        stop_event.set()
        thread.join(timeout=1)

def track_peak_memory(func):
    """
    Decorator to track peak memory usage during function execution.
    
    Usage:
        @track_peak_memory
        def my_function(logger, ...):
            # function code
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        # Find logger in args or kwargs
        logger = None
        for arg in args:
            if isinstance(arg, logging.Logger):
                logger = arg
                break
        if not logger and 'logger' in kwargs:
            logger = kwargs['logger']
        if not logger:
            logger = logging.getLogger('kraken_analysis')

        # Record starting memory
        start_mem = get_memory_usage()
        logger.info(f"Starting memory: {start_mem:.2f} MB")

        # Track peak memory
        peak_mem = start_mem
        
        def memory_checker():
            nonlocal peak_mem
            current = get_memory_usage()
            if current > peak_mem:
                peak_mem = current

        stop_event = threading.Event()

        def check_loop():
            while not stop_event.is_set():
                memory_checker()
                time.sleep(1)

        monitor_thread = threading.Thread(target=check_loop, daemon=True)
        monitor_thread.start()

        try:
            result = func(*args, **kwargs)  # Run the decorated function
            return result
        finally:
            # Stop monitoring
            stop_event.set()
            monitor_thread.join(timeout=1)
            
            # Log results
            end_mem = get_memory_usage()
            logger.info(
                f"Memory usage: start={start_mem:.2f} MB, "
                f"peak={peak_mem:.2f} MB, end={end_mem:.2f} MB"
            )
    return wrapper

def limit_memory_usage(max_memory_mb=None):
    """
    Attempt to limit memory usage of the current process (Unix/Linux only).
    
    Args:
        max_memory_mb: Maximum memory in MB (None = no limit).
    
    Returns:
        True if the limit was successfully set, False otherwise.
    """
    if max_memory_mb is None:
        return False

    try:
        max_memory_bytes = max_memory_mb * 1024 * 1024
        # Set both soft and hard limits
        resource.setrlimit(resource.RLIMIT_AS, (max_memory_bytes, max_memory_bytes))
        return True
    except (ImportError, ValueError, OSError, resource.error):
        return False

def calculate_optimal_resources(available_threads, num_samples, min_threads_per_sample=1):
    """
    Calculate optimal thread allocation between samples and processes.
    
    Args:
        available_threads: Total available CPU threads.
        num_samples: Number of samples to process.
        min_threads_per_sample: Minimum threads to allocate per sample.
        
    Returns:
        A tuple (threads_per_sample, max_parallel_samples).
    """
    if available_threads is None:
        available_threads = multiprocessing.cpu_count()
    
    # Start with maximum parallelism
    max_parallel = min(num_samples, available_threads)
    threads_per_sample = max(min_threads_per_sample, available_threads // max_parallel)
    
    # Recalculate max_parallel based on threads_per_sample
    max_parallel = min(max_parallel, available_threads // threads_per_sample)
    
    return threads_per_sample, max_parallel

def log_resource_usage(logger, sample_id=None):
    """
    Log current resource usage (memory and CPU).
    
    Args:
        logger: Logger instance to write the usage info.
        sample_id: Optional sample identifier to include in the log.
    """
    mem_usage = get_memory_usage()
    cpu_percent = psutil.cpu_percent()
    
    message = f"Resource usage: Memory: {mem_usage:.2f} MB, CPU: {cpu_percent}%"
    if sample_id:
        message = f"Sample {sample_id}: {message}"
    
    logger.info(message)

def estimate_memory_requirements(num_samples, tool="kraken2"):
    """
    Estimate memory requirements for a tool based on the number of samples.
    
    Args:
        num_samples: Number of samples to process.
        tool: Tool name ("kneaddata", "kraken2", or "bracken").
        
    Returns:
        Estimated memory in MB for the given number of samples.
    """
    # These are rough estimates and may need adjustment based on real-world usage
    if tool.lower() == "kneaddata":
        # KneadData is somewhat less memory-intensive
        return 2000 + (500 * num_samples)  # Base 2GB + 500MB per sample
    elif tool.lower() == "kraken2":
        # Kraken2 can be quite memory-intensive depending on the database
        return 8000 + (1000 * num_samples)  # Base 8GB + 1GB per sample
    elif tool.lower() == "bracken":
        # Bracken typically uses less memory than Kraken2
        return 4000 + (500 * num_samples)  # Base 4GB + 500MB per sample
    else:
        # Generic estimate if tool is unrecognized
        return 4000 + (1000 * num_samples)

def check_resource_availability(required_memory, required_threads):
    """
    Check if the system has enough memory and CPU threads available.
    
    Args:
        required_memory: Required memory in MB.
        required_threads: Required CPU threads.
        
    Returns:
        A tuple (memory_ok, cpu_ok, available_memory, available_threads) indicating:
          - memory_ok: bool, whether enough memory is available
          - cpu_ok: bool, whether enough CPU threads are available
          - available_memory: total available system memory in MB
          - available_threads: total available CPU threads
    """
    available_memory = psutil.virtual_memory().available / (1024 * 1024)  # MB
    available_threads = multiprocessing.cpu_count()
    
    memory_ok = available_memory >= required_memory
    cpu_ok = available_threads >= required_threads
    
    return memory_ok, cpu_ok, available_memory, available_threads