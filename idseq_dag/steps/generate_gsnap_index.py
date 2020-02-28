''' Generate GSNAP index given NT '''

import os
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.log as log
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns

class PipelineStepGenerateGsnapIndex(PipelineStep):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.upload_results_with_checksum = True

    ''' Generate  gsnap index '''
    def run(self):
        """
          Generate GSNAP index. To be called from idseq-infra

        """
        nt_db = self.input_files_local[0][0]
        output_nt_index_parent_dir = self.output_dir_local
        output_nt_index_dir = self.additional_attributes["output_nt_index"]
        k = self.additional_attributes.get("k", 16)  # kmer k
        log.write(f"input: {nt_db} output: {output_nt_index_dir}")
        # DUMMY OPTS
        command.execute(
            command_patterns.SingleCommand(
                cmd="mkdir",
                args=[output_nt_index_dir],
                cd=output_nt_index_parent_dir
            )
        )
        command.execute(
            command_patterns.SingleCommand(
                cmd="touch",
                args=[os.path.join(output_nt_index_dir, "foo")],
                cd=output_nt_index_parent_dir
            )
        )
        self.additional_output_folders_hidden.append(output_nt_index_dir)

    def count_reads(self):
        ''' Count reads '''
        pass
