''' Generate GSNAP index given NT '''

import os
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.log as log
import idseq_dag.util.command as command
import idseq_dag.util.count as count

class PipelineStepGenerateGsnapIndex(PipelineStep):
    ''' Generate  gsnap index '''
    def run(self):
        """
          Generate GSNAP index. To be called from idseq-infra

        """
        nt_db = self.input_files_local[0][0]
        output_nt_index_tar = self.output_files_local()[0]
        k = self.additional_attributes.get("k", 16) # kmer k
        PipelineStepGenerateGsnapIndex.generate_gsnap_index(nt_db, k, output_nt_index_tar)

    @staticmethod
    def generate_gsnap_index(input_db, k, output_tar):
        output_index_dir = os.path.dirname(output_tar)
        output_base = os.path.basename(output_tar)
        log.write(f"input: {input_db} output: {output_tar}")
        command.execute(f"gmap_build -D {output_index_dir} -d {output_base[:-4]} -k {k} {input_db} ")
        command.execute(f"cd {output_index_dir}; tar cvf {output_base} {output_base[:-4]}")


    def count_reads(self):
        ''' Count reads '''
        pass
