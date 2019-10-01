"""Generate lz4 file given input file"""

import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.log as log

from idseq_dag.engine.pipeline_step import PipelineStep

class PipelineStepGenerateLZ4(PipelineStep):

    def run(self):
        # lz4 is "scalable with multi-cores CPU" so we let it parallelize
        # itself. See https://github.com/lz4/lz4 .
        for file_list in self.input_files:
            for input_file in file_list:
                if input_file.endswith(('.gz', '.zip', '.lz4')):
                    log.write(f'Skipping already-compressed file {input_file}')
                else:
                    command.execute(self.get_command(input_file))

    def get_command(self, input_file):
        log.write(f"input: {input_file} output: {input_file}.lz4")
        return command_patterns.SingleCommand(
            cmd="lz4",
            args=[
                "-9", # max compression
                "-f", # force overwrite output file
                input_file
            ]
        )

    def count_reads(self):
        pass
