import logging
import multiprocessing
import os
import sys
import threading
import idseq_dag.util.command as command

print_lock = multiprocessing.RLock()


def configure_logger(log_file):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    handler = logging.FileHandler(log_file)
    formatter = logging.Formatter("%(asctime)s: %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # Echo to stdout so they get to CloudWatch
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def upload_log_file(sample_s3_output_path, lock=threading.RLock()):
    with lock:
        logh = logging.getLogger().handlers[0]
        logh.flush()
        command.execute("aws s3 cp --quiet %s %s/;" % (logh.baseFilename,
                                                       sample_s3_output_path))


def write_to_log(message, warning=False, flush=True):
    logger = logging.getLogger()
    with print_lock:
        if warning:
            logger.warning(message)
        else:
            logger.info(message)
        if flush:
            sys.stdout.flush()


def set_up_stdout():
    # Unbuffer stdout and redirect stderr into stdout. This helps observe logged
    # events in realtime.
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    os.dup2(sys.stdout.fileno(), sys.stderr.fileno())
