import unittest

from idseq_dag.engine.pipeline_flow import PipelineFlow
import pdb

from idseq_dag.steps.run_star import PipelineStepRunStar
from idseq_dag.util import command
import idseq_dag.util.test as test
from tests.idseq_step_setup import IdseqStepSetup


class CustomTest1(unittest.TestCase):
    def test_custom1(self):
        # TODO: Check results
        # TODO: Clean up the folder

        runstep = IdseqStepSetup.get_step_object(PipelineStepRunStar, "star_out", paired=True, dag_file="examples/custom_test_1.json")
        runstep.start()
        runstep.wait_until_finished()
        actual = runstep.output_files_local()
        A_actual = actual[0]
        truth_base = "s3://idseq-samples-staging/samples-test/truth/"
        A_expected = truth_base + runstep.output_files[0]
        test.should_match_exactly(A_expected, A_actual)
