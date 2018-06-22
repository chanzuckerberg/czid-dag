import unittest

from idseq_dag.engine.pipeline_flow import PipelineFlow
import pdb

from idseq_dag.steps.run_star import PipelineStepRunStar
from idseq_dag.util import command
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
        self.compare(A_expected, A_actual)

        

    def compare(self, expected, actual):
        expected_content = command.execute_with_output(f"aws s3 cp {expected} -")
        actual_content = command.execute_with_output(f"cat {actual}")
        if expected_content == actual_content:
            print(f"File {expected} is the same as {actual}")
        else:
            raise RuntimeError(f"{actual} not equal to {expected}")
