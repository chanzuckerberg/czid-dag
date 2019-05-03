import multiprocessing
import threading
import idseq_dag.util.log as log

class TraceLock():
    def __init__(self, lock_name, lock=multiprocessing.RLock()):
        self._lock = lock
        self._lock_name = lock_name

    def acquire(self):
        v = {"lock_name": self._lock_name, "thread_name": threading.current_thread().name}
        log.log_event("trace_lock", values={**v, "state": "acquiring"})
        if self._lock.acquire(False):
            log.log_event("trace_lock", values={**v, "state": "acquired"})
        else:
            log.log_event("trace_lock", values={**v, "state": "waiting"})
            self._lock.acquire(True)
            log.log_event("trace_lock", values={**v, "state": "acquired_after_wait"})

    def __enter__(self):
        self.acquire()

    def release(self):
        log.log_event("trace_lock", values={"lock_name": self._lock_name,
                                            "thread_name": threading.current_thread().name,
                                            "state": "release"})
        self._lock.release()

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.release()