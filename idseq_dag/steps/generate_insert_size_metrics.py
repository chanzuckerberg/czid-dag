import multiprocessing
import os
import json
import re

from idseq_dag.engine.pipeline_step import PipelineStep, StepStatus, InputFileErrors
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns


class GenerateInsertSizeMetrics(PipelineStep):
    """ Implements the step for generating insert size metrics.

    Picard computes insert size size metrics.
    It runs out the output bam from the STAR step.

    ```
    java -jar picard.jar CollectInsertSizeMetrics
        I=Aligned.out.sam
        O=picard_insert_metrics.txt
        H=insert_size_histogram.pdf
    ```

    picard documentation can be found [here](https://broadinstitute.github.io/picard/)
    """


    def run(self):
        """Run picard to generate insert metrics."""

        cd = self.output_dir_local
        cmd = 'picard'

        input_file = self.input_files_local[0][2]
        metrics_file = self.output_files[0]
        histogram_file = self.output_files[1]

        params = [
            "CollectInsertSizeMetrics",
            f"I={input_file}",
            f"O={metrics_file}",
            f"H={histogram_file}"
        ]

        command.execute(
            command_patterns.SingleCommand(
                cd=cd,
                cmd=cmd,
                args=params
            )
        )

    def count_reads(self):
        pass

