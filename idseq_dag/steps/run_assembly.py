''' Run PriceSeq Filter '''
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.count as count

class PipelineStepRunAssembly(PipelineStep):
    ''' PriceSeq PipelineStep implementation '''
    def run(self):
        """
           Run Assembly
        """
        input_fasta = self.input_files_local[0][-1]
        assembled_contig, assembled_scaffold, bowtie_sam, contig_stats = self.output_files_local()
        read2contig = {}
        assembled_dir =
        try:

        except:
            # Common Assembly Error
            traceback.print_exc()

    @staticmethod
    def assemble(input_fasta,
                 assembled_contig,
                 assembled_scaffold,
                 bowtie_sam,
                 contig_stats,
                 read2contig,
                 assemble_dir,
                 check_presence = True):
            assembled_dir = os.path.join(genus_dir, 'spades')
            output = [None, None, None, None]
            command.execute(f"mkdir -p {assembled_dir}")
            assembled_contig_tmp = os.path.join(assembled_dir, 'contigs.fasta')
            assembled_scaffold_tmp = os.path.join(assembled_dir, 'scaffolds.fasta')
            assembled_contig = os.path.join(genus_dir, 'contigs.fasta')
            assembled_scaffold = os.path.join(genus_dir, 'scaffolds.fasta')

            try:
                if s3.check_s3_presence(self.s3_path(assembled_contig)) and \
                    s3.check_s3_presence(self.s3_path(assembled_scaffold)):
                    # check if file already assembled before, if so, reuse.
                    assembled_contig = s3.fetch_from_s3(self.s3_path(assembled_contig),
                                                        genus_dir)
                    assembled_scaffold = s3.fetch_from_s3(self.s3_path(assembled_scaffold),
                                                          genus_dir)
                else:
                    command.execute(f"spades.py -s {fasta_file} -o {assembled_dir} -m 60 -t 32 --only-assembler")
                    command.execute(f"mv {assembled_contig_tmp} {assembled_contig}")
                    command.execute(f"mv {assembled_scaffold_tmp} {assembled_scaffold}")

                # build the bowtie index based on the contigs
                (read2contig, bowtie_sam) = self.generate_read_to_contig_mapping(assembled_contig, fasta_file)
                output = [assembled_contig, assembled_scaffold, read2contig, bowtie_sam]
    @staticmethod
    def generate_read_to_contig_mapping(assembled_contig, fasta_file):
        ''' read -> contig mapping through bowtie2 alignment '''
        genus_dir = os.path.dirname(fasta_file)
        # build bowtie index based on assembled_contig
        bowtie_index_path = os.path.join(genus_dir, 'bowtie-contig')
        bowtie_sam = os.path.join(genus_dir, "read-contig.sam")
        command.execute(f"mkdir -p {bowtie_index_path}; bowtie2-build {assembled_contig} {bowtie_index_path}")
        command.execute(f"bowtie2 -x {bowtie_index_path} -f -U {fasta_file} --very-sensitive -p 32 > {bowtie_sam}")
        read2contig = {}
        with open(bowtie_sam, "r", encoding='utf-8') as samf:
            for line in samf:
                if line[0] == '@':
                    continue
                fields = line.split("\t")
                read = fields[0]
                contig = fields[2]
                if contig != '*':
                    read2contig[read] = contig
        return (read2contig, bowtie_sam)

    def count_reads(self):
        ''' count reads '''
        pass
