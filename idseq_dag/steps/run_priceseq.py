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

        input_fqs = self.input_files
        FILE_TYPE = ""
        PRICESEQFILTER_OUT1 = 'priceseqfilter.unmapped.star.1.fq'
        PRICESEQFILTER_OUT2 = 'priceseqfilter.unmapped.star.2.fq'

        def result_path(basename):
            return os.path.join(self.output_dir_local, basename)

        # PriceSeqFilter determines input type based on extension.
        # It will throw an exception if output extension doesn't match input.
        correct_file_extension = os.path.splitext(FILE_TYPE)[0]
        input_files = ["%s.%s" % (fq, correct_file_extension) for fq in input_fqs]
        output_files = [
            "%s_priceseqfilter_output.%s" % (f, correct_file_extension)
            for f in input_files
        ]

        for fq, f in zip(input_fqs, input_files):
            command.execute("ln %s %s" % (fq, f))
        priceseq_params = ["PriceSeqFilter", '-a', '12', '-rnf', '90', '-log', 'c']
        if len(input_fqs) == 2:
            priceseq_params.extend([
                '-fp', input_files[0], input_files[1], '-op', output_files[0],
                output_files[1]
            ])
        else:
            priceseq_params.extend(['-f', input_files[0], '-o', output_files[0]])
        if "fastq" in FILE_TYPE:
            priceseq_params.extend(['-rqf', '85', '0.98'])
        command.execute(" ".join(priceseq_params))
        log.write("Finished running PriceSeqFilter.")

        out_path = result_path(PRICESEQFILTER_OUT1)
        command.execute("mv %s %s" % (output_files[0], out_path))
        s3_dst = self.output_dir_s3 + "/"
        s3.upload_with_retries(out_path, s3_dst)

        if len(input_fqs) == 2:
            out_path = result_path(PRICESEQFILTER_OUT2)
            command.execute("mv %s %s" % (output_files[1], out_path))
            s3.upload_with_retries(out_path, s3_dst)
