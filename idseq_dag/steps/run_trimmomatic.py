''' Run Trimmomatic '''
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.count as count

class PipelineStepRunPriceSeq(PipelineStep):
    ''' Trimmomatic PipelineStep implementation '''
    def run(self):
        """
        Trim low-quality ends from the reads.
        Illumina read quality tends to deteriorate towards the 3' end in both read-1 and read-2, particularly read-2.
        So cut off the end after the quality starts to drop below a threshold.
        Discard any reads that become too short.
        Also trim any residual Illumina adapters.

        java -jar trimmomatic-0.38.jar PE -phred33
          input_R1.fq.gz input_R2.fq.gz
          output_R1_paired.fq.gz output_R1_unpaired.fq.gz
          output_R2_paired.fq.gz output_R2_unpaired.fq.gz
          ILLUMINACLIP:TruSeq3-PE.fa:2:30:10
          LEADING:25 TRAILING:25 SLIDINGWINDOW:4:25 MINLEN:35

        See: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
        """
        input_files = self.input_files_local[0][0:2]
        output_files = self.output_files_local()
        is_paired = (len(input_files) == 2)

        # PriceSeqFilter determines input type based on extension. It will
        # throw an exception if output extension doesn't match input
        # extension.
        file_type = self.additional_attributes.get("file_type", "fastq")
        price_out = [
            f"{f}_priceseqfilter_output.{file_type}" for f in input_files
        ]

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
        if file_type != 'fasta' and file_type != 'fa':
            # Fastq
            self.fq2fa(price_out[0], output_files[0])
            if is_paired:
                self.fq2fa(price_out[1], output_files[1])
        else:
            command.execute(f"mv {price_out[0]} {output_files[0]}")
            if is_paired:
                command.execute(f"mv {price_out[1]} {output_files[1]}")

    def count_reads(self):
        ''' Count reads '''
        self.should_count_reads = True
        self.counts_dict[self.name] = count.reads_in_group(self.output_files_local()[0:2])

    @staticmethod
    def fq2fa(input_fastq, output_fasta):
        ''' FASTQ to FASTA conversion '''
        step = "FASTQ to FASTA conversion"
        log.write(f"Starting {step}...")
        cmd = f"sed -n '1~4s/^@/>/p;2~4p' <{input_fastq} >{output_fasta}"
        command.execute(cmd)
        log.write(f"Finished {step}.")
