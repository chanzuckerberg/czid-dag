import os

from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.s3 as s3

class PipelineStepRunPriceSeq(PipelineStep):
    def run(self):
        """PriceSeqFilter is used to filter input data based on quality. Two FASTQ
        inputs means paired reads.

        See: http://derisilab.ucsf.edu/software/price/
        """
        for f in self.output_files_local():
            command.execute("date > %s" % f)

        input_fqs = self.input_files[0]
        FILE_TYPE = "" # additional attributes

        # PriceSeqFilter determines input type based on extension.
        # It will throw an exception if output extension doesn't match input.
        correct_file_extension = os.path.splitext(FILE_TYPE)[0]
        input_files = ["%s.%s" % (fq, correct_file_extension) for fq in input_fqs]

        for fq, f in zip(input_fqs, input_files):
            command.execute("ln %s %s" % (fq, f))
        priceseq_params = ["PriceSeqFilter", '-a', '12', '-rnf', '90', '-log', 'c']
        if len(input_fqs) == 2:
            priceseq_params.extend([
                '-fp', input_files[0], input_files[1], '-op', self.output_files_local()[0],
                self.output_files_local()[1]
            ])
        else:
            priceseq_params.extend(['-f', input_files[0], '-o', self.output_files_local()[0]])
        if "fastq" in FILE_TYPE:
            priceseq_params.extend(['-rqf', '85', '0.98'])
        command.execute(" ".join(priceseq_params))
        log.write("Finished running PriceSeqFilter.")
