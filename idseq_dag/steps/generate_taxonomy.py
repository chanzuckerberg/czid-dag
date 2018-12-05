''' Generate Taxonomy  '''
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.count as count

class PipelineStepGenerateTaxonomy(PipelineStep):
    ''' Generate Taxonomy '''
    def run(self):
        """

        """
        input_files = self.input_files_local[0][0:2]
        output_files = self.output_files_local()


    def count_reads(self):
        ''' Count reads '''
        pass

