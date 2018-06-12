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
        # It will throw an exception if output extension doesn't match input extension.
        file_type = self.additional_attributes["file_type"]
        correct_file_extension = os.path.splitext(file_type)[0]
        input_files = ["%s.%s" % (fq, correct_file_extension) for fq in input_fqs]
        priceseq_out = [
            "%s_priceseqfilter_output.%s" % (f, correct_file_extension)
            for f in input_files
        ]

        for fq, f in zip(input_fqs, input_files):
            command.execute("ln %s %s" % (fq, f))
        priceseq_params = ["PriceSeqFilter", '-a', '12', '-rnf', '90', '-log', 'c']
        if len(input_files) == 2:
            priceseq_params.extend([
                '-fp', input_files[0], input_files[1], '-op', priceseq_out[0],
                priceseq_out[1]
            ])
        else:
            priceseq_params.extend(['-f', input_files[0], '-o', priceseq_out[0]])
        if "fastq" in file_type:
            priceseq_params.extend(['-rqf', '85', '0.98'])
        command.execute(" ".join(priceseq_params))

        # Run FASTQ to FASTA if needed
        if "fastq" in file_type:
            self.fq2fa(priceseq_out[0], self.output_files_local()[0])
            if len(input_fqs) == 2:
                self.fq2fa(priceseq_out[1], self.output_files_local()[1])
            log.write("Finished FASTQ to FASTA conversion.")

        log.write("Finished running PriceSeqFilter.")

    def fq2fa(self, input_fastq, output_fasta):
        cmd = "sed -n '1~4s/^@/>/p;2~4p' <%s >%s" % (input_fastq, output_fasta)
        command.execute(cmd)
