import re
import subprocess
import idseq_dag.util.counting
from idseq_dag.engine.pipeline_step import PipelineStep

class PipelineStepTruncate(PipelineStep):
    def run(self):
        total_reads = 0
        read_counts_observed = set()
        for i, input in enumerate(self.input_files_local[0]):
            output = self.output_files_local()[i]

            # Unzip input
            if input[-3:] == '.gz':
              subprocess.check_call("gunzip {}".format(input), shell=True)
              unzipped_input = input[:-3]
            else:
              unzipped_input = input

            # Count reads and truncate
            file_format = unzipped_input.split(".")[-1]
            assert file_format in ["fq", "fastq"]
            max_lines = idseq_dag.util.counting.reads2lines(self.additional_attributes["truncate_reads_to"],
                                                            file_format)
            input_lines = subprocess.check_output("wc -l {}".format(unzipped_input), shell=True)
            if input_lines > max_lines:
                subprocess.check_call("head -{} {} > {}".format(max_lines, unzipped_input, output), shell=True)
                read_count = idseq_dag.util.counting.lines2reads(max_lines)
            else:
                subprocess.check_call("mv {} {}".format(unzipped_input, output), shell=True)
                read_count = idseq_dag.util.counting.lines2reads(input_lines, file_format)
            total_reads += read_count
            read_counts_observed.add(read_count)

            # Validity check
            assert len(read_counts_observed) == 1, "Differing read counts in inputs: {}".format(read_counts_observed)
