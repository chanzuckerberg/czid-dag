import os
import sys
import tempfile
import unittest

from idseq_dag.steps.run_validate_input import PipelineStepRunValidateInput
from idseq_dag.exceptions import InsufficientReadsError, InvalidFileFormatError

fasta_filename = os.path.join(os.path.dirname(__file__), "fixtures", "reads.fasta")
fastq_filename = os.path.join(os.path.dirname(__file__), "fixtures", "reads.fastq")

class TestValidateInput(unittest.TestCase):
    def test_validate_input(self):
        PipelineStepRunValidateInput.quick_check_file(None, fasta_filename, is_fastq=False, max_fragments_to_check=100)
        PipelineStepRunValidateInput.quick_check_file(None, fastq_filename, is_fastq=True, max_fragments_to_check=100)
        with self.assertRaises(FileNotFoundError):
            PipelineStepRunValidateInput.quick_check_file(None, "nonexistent", is_fastq=True, max_fragments_to_check=100)
        with tempfile.NamedTemporaryFile() as tf:
            with self.assertRaises(InsufficientReadsError):
                PipelineStepRunValidateInput.quick_check_file(None, tf.name, is_fastq=True, max_fragments_to_check=100)
            tf.write(b">@" * 9000)
            tf.flush()
            with self.assertRaises(InvalidFileFormatError):
                PipelineStepRunValidateInput.quick_check_file(None, tf.name, is_fastq=True, max_fragments_to_check=100)
            tf.seek(0)
            tf.write(">read\n\n".encode())
            tf.flush()
            with self.assertRaises(InvalidFileFormatError):
                PipelineStepRunValidateInput.quick_check_file(None, tf.name, is_fastq=True, max_fragments_to_check=100)
