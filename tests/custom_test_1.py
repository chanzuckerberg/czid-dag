import unittest
import os

from idseq_dag.steps.run_star import PipelineStepRunStar
from idseq_dag.steps.run_priceseq import PipelineStepRunPriceSeq
from idseq_dag.steps.run_cdhitdup import PipelineStepRunCDHitDup
from idseq_dag.steps.run_lzw import PipelineStepRunLZW
from idseq_dag.steps.run_bowtie2 import PipelineStepRunBowtie2
from idseq_dag.steps.run_subsample import PipelineStepRunSubsample
from idseq_dag.steps.run_gsnap_filter import PipelineStepRunGsnapFilter
from idseq_dag.steps.generate_taxid_fasta import PipelineStepGenerateTaxidFasta
from idseq_dag.steps.generate_taxid_locator import PipelineStepGenerateTaxidLocator
from idseq_dag.steps.generate_alignment_viz import PipelineStepGenerateAlignmentViz
from idseq_dag.steps.combine_taxon_counts import PipelineStepCombineTaxonCounts
from idseq_dag.steps.generate_annotated_fasta import PipelineStepGenerateAnnotatedFasta
import idseq_dag.util.test as test
from tests.idseq_step_setup import IdseqStepSetup


class CustomTest1(unittest.TestCase):
    def test_custom1(self):
        dag_file = "examples/generic_test_dag.json"
        test_set = "s3://idseq-samples-test/test-sets/RR004_water_2_S23"
        output_dir_s3 = "s3://idseq-samples-test/scratch"

        step_pairs = [
            # (PipelineStepRunStar, "star_out"),
            # (PipelineStepRunPriceSeq, "priceseq_out"),
            # (PipelineStepRunCDHitDup, "cdhitdup_out"),
            # (PipelineStepRunLZW, "lzw_out"),
            # (PipelineStepRunBowtie2, "bowtie2_out"),
            # (PipelineStepRunSubsample, "subsampled_out"),
            # (PipelineStepRunGsnapFilter, "gsnap_filter_out"),
        ]

        for p in step_pairs:
            runstep = IdseqStepSetup.get_test_step_object(
                p[0], p[1], dag_file, test_set, output_dir_s3)
            runstep.start()
            runstep.wait_until_finished()
            runstep.output_files.append(f"{p[1]}.count")
            expected_files = [
                os.path.join(test_set, f) for f in runstep.output_files
            ]
            actual_files = runstep.output_files_local()

            for expected, actual in zip(expected_files, actual_files):
                test.should_match_exactly(expected, actual)

        step_pairs = [
            # (PipelineStepCombineTaxonCounts, "taxon_count_out"),
            # (PipelineStepGenerateAnnotatedFasta, "annotated_out"),
            (PipelineStepGenerateTaxidFasta, "taxid_fasta_out")
            # (PipelineStepGenerateTaxidLocator, "taxid_locator_out")
        ]

        for p in step_pairs:
            runstep = IdseqStepSetup.get_test_step_object(
                p[0], p[1], dag_file, test_set, output_dir_s3)
            runstep.start()
            runstep.wait_until_finished()
            expected_files = [
                os.path.join(test_set, f) for f in runstep.output_files
            ]
            actual_files = runstep.output_files_local()

            for expected, actual in zip(expected_files, actual_files):
                test.should_match_exactly(expected, actual)

        # Verify align viz output
        # runstep = IdseqStepSetup.get_test_step_object(
        #     PipelineStepGenerateAlignmentViz, "alignment_viz_out", dag_file, test_set, output_dir_s3)
        # runstep.start()
        # runstep.wait_until_finished()
        # actual_files = runstep.additional_files_to_upload
        # expected_files = [
        #     test_set + "/align_viz/" + os.path.basename(f) for f in runstep.additional_files_to_upload
        # ]
        #
        # for expected, actual in zip(expected_files, actual_files):
        #     test.should_match_exactly(expected, actual)
