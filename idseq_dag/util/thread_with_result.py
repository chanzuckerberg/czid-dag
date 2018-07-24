#!/bin/env python3
import threading
import traceback
import random
import time

class ThreadWithResult(threading.Thread):
    """A simple replacement for Python's built-in threading.Thread class,
    adding support for target funcs that return results or raise exceptions.
    Usage example:

        t = ThreadWithResult(target=fetch_something)
        t.start()
        t.join()
        assert t.completed and not t.exception
        print(t.result)
    """
    def __init__(self, target, args=(), kwargs=None):
        super(ThreadWithResult, self).__init__()
        self.args = args
        if kwargs is None:
            # This is a bit arcane, but we want to create a new empty dict
            # here, rather than having kwargs={} above, because in python
            # {} default values are singleton globals.  This matters
            # only in the most rare and exceptional circumstances.
            self.kwargs = {}
        else:
            self.kwargs = kwargs
        self.target = target
        self.exception = None
        self.completed = False
        self.result = None
        self.print_traceback = True

    def run(self):
        try:
            self.result = self.target(*self.args, **self.kwargs)
            self.exception = False
        except: #pylint: disable=bare-except
            if self.print_traceback:
                traceback.print_exc()
            self.exception = True
        finally:
            self.completed = True

def execute_all(threads):
    """Run the provided threads.  If all complete without raising an exception,
    return a list of their results."""
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    for i, t in enumerate(threads):
        assert t.completed and not t.exception, f"Problem in thread {i}."
    return [t.result for t in threads]

# ONLY TESTS BELOW THIS LINE
# Run this module as a command to trigger the unit test.

def unit_test(r=random.Random(time.time())):
    "Test thread_with_result.py"
    try:
        threads_to_fail = set([2, 3, 5, 7])
        def target(i):
            time.sleep(r.random())
            if i in threads_to_fail:
                raise RuntimeError()
            return i
        threads = [
            ThreadWithResult(target, (i,))
            for i in range(16)
        ]
        for i, t in enumerate(threads):
            if i in threads_to_fail:
                t.print_traceback = False
        try:
            execute_all(threads)
        except AssertionError:
            pass
        except:  #pylint: disable=bare-except
            print("Unexpected exception.")
            raise
        else:
            assert False, "execute_all should have raised an AssertionError"
        for i, t in enumerate(threads):
            assert t.completed
            if i in threads_to_fail:
                assert t.exception, f"Thread {i} should have raised an exception."
                assert t.result is None, f"Thread {i} should not have returned a result."
            else:
                assert not t.exception, f"Thread {i} should not have raised an exception."
                assert t.result == i, f"Thread {i} should have returned result {i}."
        print("Tests passed.")
    except:
        print("Tests failed.")
        raise

if __name__ == "__main__":
    unit_test()
