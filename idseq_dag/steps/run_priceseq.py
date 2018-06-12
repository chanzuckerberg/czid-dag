from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.log as log


class PipelineStepRunPriceSeq(PipelineStep):
    def run(self):
        """PriceSeqFilter is used to filter input data based on quality. Two FASTQ
        inputs means paired reads.

        See: http://derisilab.ucsf.edu/software/price/
        """
        original_inputs = self.input_files_local[0][0:2]
        output_files = self.output_files_local()
        is_paired = (len(original_inputs) == 2)

        # PriceSeqFilter determines input type based on extension. It will
        # throw an exception if output extension doesn't match input
        # extension. Set the extensions to the same thing and use ln instead
        # of renaming files.
        file_type = self.additional_attributes.get("file_type", "fastq")
        input_files = [f"{orig}.{file_type}" for orig in original_inputs]
        price_out = [
            f"{file}_priceseqfilter_output.{file_type}" for file in input_files
        ]
        for orig, file in zip(original_inputs, input_files):
            command.execute(f"ln {orig} {file}")

        params = ["PriceSeqFilter", '-a', '12', '-rnf', '90', '-log', 'c']
        if is_paired:
            params.extend([
                '-fp', input_files[0], input_files[1], '-op', price_out[0],
                price_out[1]
            ])
        else:
            params.extend(['-f', input_files[0], '-o', price_out[0]])
        if "fasta" not in file_type:  # Default fastq. Explicitly specify fasta.
            params.extend(['-rqf', '85', '0.98'])
        cmd = " ".join(params)
        command.execute(cmd)

        # Run FASTQ to FASTA if needed
        if "fasta" not in file_type:
            self.fq2fa(price_out[0], output_files[0])
            if is_paired:
                self.fq2fa(price_out[1], output_files[1])
            log.write("Finished FASTQ to FASTA conversion.")
        log.write("Finished running PriceSeqFilter.")

    def fq2fa(self, input_fastq, output_fasta):
        # FASTQ to FASTA conversion
        cmd = f"sed -n '1~4s/^@/>/p;2~4p' <{input_fastq} >{output_fasta}"
        command.execute(cmd)
