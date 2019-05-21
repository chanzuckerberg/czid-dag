import unittest
from time import sleep
from unittest.mock import patch, call, ANY
from multiprocessing import Value


# Module under test
import idseq_dag.util.command as command

class TestCommand(unittest.TestCase):
    '''Tests for `util/command.py`'''

    @patch('idseq_dag.util.log.log_event')
    def test_run_in_subprocess_invoke_db_hack_only_once(self, mock_log_event):
        '''Test if run in subprocess invokes db_hack only once'''
        multiprocessing_counter = Value('i', 0)
        def increment_counter(counter):
            with counter.get_lock():
                counter.value += 1
        command._db_hack_invocation_flag.value = False

        command.run_in_subprocess(increment_counter)(multiprocessing_counter)
        command.run_in_subprocess(increment_counter)(multiprocessing_counter)

        self.assertEqual(multiprocessing_counter.value, 2)
        mock_log_event.assert_has_calls([
            call('ctx_start', None, extra_fields={'context_name': 'db_hack', 'caller': ANY}),
            call('ctx_end', None, extra_fields={'context_name': 'db_hack', 'caller': ANY}, start_time=ANY)
        ])
        self.assertEqual(mock_log_event.call_count, 2)
