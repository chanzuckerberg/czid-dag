import os

from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.log as log


class PipelineStepRunPriceSeq(PipelineStep):
    def run(self):
        for f in self.output_files_local():
            command.execute("date > %s" % f)

        def result_path(basename):
            return os.path.join(self.output_dir_local, basename)

        """PriceSeqFilter is used to filter input data based on quality. Two FASTQ
        inputs means paired reads.
    
        See: http://derisilab.ucsf.edu/software/price/
        """
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
        priceseq_params = [PRICESEQ_FILTER, '-a', '12', '-rnf', '90', '-log', 'c']
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
        log.write("mv %s %s" % (output_files[0], out_path))
        s3_dst = SAMPLE_S3_OUTPUT_PATH + "/"
        uploader_start(out_path, s3_dst)

        if len(input_fqs) == 2:
            out_path = result_path(PRICESEQFILTER_OUT2)
            command.execute("mv %s %s" % (output_files[1], out_path))
            uploader_start(out_path, s3_dst)
