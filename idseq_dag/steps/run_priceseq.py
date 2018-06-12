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
        print("start of PipelineStepRunPriceSeq")
        print("Input files: " + str(self.input_files))
        print("Output files: " + str(self.output_files))
        print("Output dir local: " + self.output_dir_local)
        print("Output dir s3: " + self.output_dir_s3)
        print("Ref dir local: " + self.ref_dir_local)
        print("Additional files: " + str(self.additional_files))
        print("Additional attributes: " + str(self.additional_attributes))
        print("Output files local: " + str(self.output_files_local()))
        print("Input files local: " + str(self.input_files_local))

        input_fqs = self.input_files_local[0][0:2]

        # PriceSeqFilter determines input type based on extension.
        # It will throw an exception if output extension doesn't match input.
        file_type = self.additional_attributes["file_type"]
        correct_file_extension = os.path.splitext(file_type)[0]
        input_files = ["%s.%s" % (fq, correct_file_extension) for fq in input_fqs]
        output_files = [
            "%s_priceseqfilter_output.%s" % (f, correct_file_extension)
            for f in input_files
        ]

        for fq, f in zip(input_fqs, input_files):
            command.execute("ln %s %s" % (fq, f))
        priceseq_params = ["PriceSeqFilter", '-a', '12', '-rnf', '90', '-log', 'c']
        if len(input_files) == 2:
            priceseq_params.extend([
                '-fp', input_files[0], input_files[1], '-op', self.output_files[0],
                self.output_files[1]
            ])
        else:
            priceseq_params.extend(['-f', input_files[0], '-o', self.output_files[0]])
        if "fastq" in input_files:
            priceseq_params.extend(['-rqf', '85', '0.98'])
        command.execute(" ".join(priceseq_params))
        log.write("Finished running PriceSeqFilter.")
    #
    #     # Run FASTQ to FASTA
    #     if "fastq" in file_type:
    #         log_params = return_merged_dict(DEFAULT_LOG_PARAMS,
    #                                         {"title": "FASTQ to FASTA"})
    #         input_files = [result_path(PRICESEQFILTER_OUT1)]
    #         next_inputs = [result_path(FQ2FA_OUT1)]
    #         if number_of_input_files == 2:
    #             input_files.append(result_path(PRICESEQFILTER_OUT2))
    #             next_inputs.append(result_path(FQ2FA_OUT2))
    #         run_and_log_s3(log_params, target_outputs["run_fq2fa"], lazy_run,
    #                        SAMPLE_S3_OUTPUT_PATH, run_fq2fa, input_files,
    #                        uploader_start)
    #     else:
    #         next_inputs = [result_path(PRICESEQFILTER_OUT1)]
    #         if number_of_input_files == 2:
    #             next_inputs.append(result_path(PRICESEQFILTER_OUT2))
    #
    #
    #     # Okay this would basically be ready to go once we resolve the file type issue and do code clean up.
    #
    # def self.run_fq2fa(input_fqs, uploader_start):
    #     """FASTQ to FASTA conversion."""
    #     fq2fa(input_fqs[0], result_path(FQ2FA_OUT1))
    #     if len(input_fqs) == 2:
    #         fq2fa(input_fqs[1], result_path(FQ2FA_OUT2))
    #     write_to_log("Finished FASTQ to FASTA conversion.")
    #
    #     dst = SAMPLE_S3_OUTPUT_PATH + "/"
    #     uploader_start(result_path(FQ2FA_OUT1), dst)
    #     if len(input_fqs) == 2:
    #         uploader_start(result_path(FQ2FA_OUT2), dst)
