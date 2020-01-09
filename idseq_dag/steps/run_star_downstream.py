from idseq_dag.steps.run_star import PipelineStepRunStar

class PipelineStepRunStarDownstream(PipelineStepRunStar):
    """ Runs STAR as a step with two input targets:
    (1) sequence files (output from any step)
        --> can be either (a) a single file, or (b) paired files plus a merged file which will be ignored
    (2) validation counts (output from validation step)
        --> first file is a count json; remaining files are sanitized sequence files which we don't need
            because we use the files from (1)
    """

    def run(self):
        self.sequence_input_files = self.input_files_local[0][:2]
        self.validated_input_counts_file = self.input_files_local[1][0]
        # Don't collect insert size metrics at this stage
        super(insert_size_metrics_enabled=False).run()
