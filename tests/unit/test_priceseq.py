import os
import sys
import tempfile
import unittest

from idseq_dag.steps.run_priceseq import PipelineStepRunPriceSeq
from idseq_dag.exceptions import InvalidInputFileError

fasta_filename = os.path.join(os.path.dirname(__file__), "fixtures", "reads.fasta")
fastq_filename = os.path.join(os.path.dirname(__file__), "fixtures", "reads.fastq")

@unittest.skipIf(os.uname().sysname != "Linux", "Skipping test on incompatible platform")
class TestCountReads(unittest.TestCase):
    def test_run_priceseqfilter(self):
        with tempfile.NamedTemporaryFile() as tf:
            PipelineStepRunPriceSeq.run_priceseqfilter(None, [fasta_filename], [tf.name], is_paired=False, file_type="fastq")
