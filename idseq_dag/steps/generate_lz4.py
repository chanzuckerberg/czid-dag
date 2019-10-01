"""Generate lz4 file given input file"""

import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.log as log

from idseq_dag.engine.pipeline_step import PipelineStep

class PipelineStepGenerateLZ4(PipelineStep):

    def run(self):
        command.execute(self.get_command())

    def get_command(self):
        input_file = self.input_files[0][0]
        output_file = self.output_files_local()[0]
        log.write(f"input: {input_file} output: {output_file}")
        return command_patterns.SingleCommand(
            cmd="lz4",
            args=[
                "-9", # max compression
                "-f", # force overwrite output file
                input_file,
                output_file
            ]
        )

    def count_reads(self):
        pass
