from idseq_dag.engine.pipeline_step import PipelineStep

class PipelineStepRunStarDownstream(PipelineStep):
    """ Runs STAR as a step with two input targets:
    (1) sequence files (output from any step)
    (2) validation counts (output from validation step)
    """ 

    def __init__(self, *args):
        super().__init__(*args)
        self.sequence_input_files = self.input_files_local[0]
        self.validated_input_counts_file = self.input_files_local[1][0]
