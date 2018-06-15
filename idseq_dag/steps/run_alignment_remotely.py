from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command

class PipelineStepRunAlignmentRemotely(PipelineStep):
    '''
    Run gsnap/rapsearch2 remotely
    '''
    def run(self):
        ''' Run alignmment remotely '''
        input_fas = self.input_files_local[0]
        output_fas = self.output_files_local()
