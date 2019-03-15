from idseq_dag.steps.run_star import PipelineStepRunStar

class PipelineStepRunStarDownstream(PipelineStepRunStar):
    """ Runs STAR as a step with two input targets:
    (1) sequence files (output from any step)
    (2) validation counts (output from validation step)
    """ 

    def run(self):
        self.sequence_input_files = self.input_files_local[0]
        self.validated_input_counts_file = self.input_files_local[1][0]
        super().run()
