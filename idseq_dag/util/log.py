import logging
import multiprocessing
import os
import sys
import json
import datetime
import functools

from contextlib import contextmanager

print_lock = multiprocessing.RLock()


def configure_logger(log_file=None):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    if log_file:
        handler = logging.FileHandler(log_file)
        formatter = logging.Formatter("%(asctime)s: [%(threadName)12s] %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    # Echo to stdout so they get to CloudWatch
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s: [%(threadName)12s] %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def write(message, warning=False, flush=True):
    logger = logging.getLogger()
    with print_lock:
        if warning:
            logger.warning(message)
        else:
            logger.info(message)
        if flush:
            sys.stdout.flush()

def log_event(event_name, values=None, start_time=None, warning=False, flush=True):
    '''
    Write a new log event.

    Parameters:
    event_name (str): name of the event
    values(dict): Optional. Values associated to that event. Will be logged in json format.
    start_time(datetime): Optional. If given, it will calc the elapsed time from this datetime to now and write it to the log line.

    Returns:
    datetime: Now. It can be used to pass to parameter start_time in a future log_event call.

    Example:
    start = log_event("downloaded_started", {"file":"abc"})
    # ... download your file here ...
    log_event("downloaded_completed", {"file":"abc"}, start_time=start)
    '''
    fmt_values = "{}" if values is None else json.dumps(values)
    if start_time is None:
        fmt_message = "Event '%s' %s" % (event_name, fmt_values)
    else:
        duration = (datetime.datetime.now()-start_time).total_seconds()
        fmt_message = "Event '%s' %s (%.1f seconds)" % (event_name, fmt_values, duration)
    write(fmt_message, warning, flush)
    return datetime.datetime.now()

@contextmanager
def log_context(context_name, values=None, log_caller_info=True):
    '''
    Log context manager to track started and completed events for a block of code

    Parameters:
    context_name (str): Name for this context.
    values(dict): Optional. Values associated to that event. Will be logged in json format.
    log_caller_info(bool): Log the caller file_name and method name.

    Example:
    # some_file.py
    from log import log_context

    def some_method
        # ... some code here ...
        with log_context("sub_step", {"abc": "123"}):
            # ... your code goes here ...

    The code above will generate two log entries:
    Event ctx_start {"n": "sub_step", "v": {abc": 123}, "caller": {"f": ["some_file.py", "some_method"]}
    Event ctx_end {"n": "sub_step", "v": {abc": 123}, "caller": {"f": ["some_file.py", "some_method"]} (0.3 sec)
    '''
    val = {"n": context_name}
    if values is not None:
        val["v"] = values
    if log_caller_info:
        f_code = sys._getframe(2).f_code
        val["caller"] = {"f": [os.path.basename(f_code.co_filename), f_code.co_name]}
    start = log_event("ctx_start", val)
    try:
        yield
        log_event("ctx_end", val, start_time=start)
    except Exception as e:
        val["error_type"] = type(e).__name__
        val["error_args"] = e.args
        log_event("ctx_error", val, start_time=start)
        raise e

def set_up_stdout():
    # Unbuffer stdout and redirect stderr into stdout. This helps observe logged
    # events in realtime.
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    os.dup2(sys.stdout.fileno(), sys.stderr.fileno())
