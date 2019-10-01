import os
import time
import unittest

from tests.unit.unittest_helpers import relative_file_path

from idseq_dag.steps.generate_lz4 import PipelineStepGenerateLZ4

INPUT_FILE = relative_file_path(__file__, 'dummy testfile.txt')
OUTPUT_FILE = INPUT_FILE + '.lz4'
INPUT_FILE_GZIP = INPUT_FILE + '.gz'
OUTPUT_FILE_GZIP = INPUT_FILE_GZIP + '.lz4'

class TestPipelineStepGenerateLZ4(unittest.TestCase):

    def setUp(self):
        self.step = PipelineStepGenerateLZ4(
            name='test_generate_lz4',
            input_files=[[INPUT_FILE, INPUT_FILE, INPUT_FILE_GZIP]],
            output_files=[OUTPUT_FILE, OUTPUT_FILE, OUTPUT_FILE_GZIP],
            output_dir_local='',
            ref_dir_local='',
            output_dir_s3='',
            additional_files={},
            additional_attributes={},
            step_status_lock=None,
            step_status_local=''
        )

    def test_get_command(self):
        self.assertEqual(
            self.step.get_command(INPUT_FILE).as_test_str(),
            'lz4 -9 -f {}'.format(INPUT_FILE),
        )

    def test_run(self):
        self.step.run()

        self.assertTrue(os.path.exists(OUTPUT_FILE))
        os.remove(OUTPUT_FILE)

        # Skip gzip files
        self.assertFalse(os.path.exists(OUTPUT_FILE_GZIP))