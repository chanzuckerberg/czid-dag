import threading
import time

class PeriodicThread(threading.Thread):
    '''
    Thread that keeps executing its target func periodically until it is told to stop.
    Usage example:
        t = PeriodicThread(target=poll_something, sleep_seconds=60, args=(some_argument_1, some_argument_2))
        t.start()
        ...
        t.stop()
        t.join()
    '''

    def __init__(self, target, sleep_seconds, args=(), kwargs=None):
        super(PeriodicThread, self).__init__()
        self.target = target
        self.sleep_seconds = sleep_seconds
        self.args = args
        self.kwargs = kwargs
        self._stop_event = threading.Event()

    def stop(self):
        self._stop_event.set()

    def stopped(self):
        return self._stop_event.is_set()

    def run(self) -> None:
        while not self.stopped():
            self.target(*self.args, **self.kwargs)
            time.sleep(self.sleep_seconds)
