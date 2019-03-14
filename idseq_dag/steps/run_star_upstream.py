from idseq_dag.engine.pipeline_step import PipelineStep

class PipelineStepRunStarUpstream(PipelineStep):
    """ Runs STAR as a step with a single input target:
    namely the output from the validation step, which includes
    both sanitized versions of the original sequence files and validation counts.
    """

    def __init__(self, *args):
        super().__init__(*args)
        self.sequence_input_files = self.input_files_local[0][1:3]
        self.validated_input_counts_file = self.input_files_local[0][0]
