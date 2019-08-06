''' Run Assembly (spades) '''
import json
import os
import traceback
from collections import defaultdict
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.count as count

class PipelineStepRunAssembly(PipelineStep):
    """ To obtain longer contigs for improved sensitivity in mapping, short reads must be 
    de novo assembled using SPADES. 
    The SPADES output loses the information about which contig each individual read belongs to. 
    Therefore, we use  bowtie2 to map the original reads onto their assembled contigs.

    First, the short reads are assembled into contigs using SPADES.

    ```
    spades.py 
    -1 {input_fasta}
    -2 {input_fasta2}
    -o {assembled_dir}
    -m {memory} 
    -t 32 
    —only-assembler
    ```

    SPADES documentation can be found [here](http://cab.spbu.ru/software/spades/)

    The single-read identity of reads merged into each contig are lost by SPADES. 
    To recover this information and identify which contig each read belongs to, 
    the contigs are then used to build a Bowtie2 database:

    ```
    bowtie2-build {assembled_contig} {bowtie_index_path}
    ```

    Finally, the original reads are mapped back to their assembled contigs:

    ```
    bowtie2 
    -x {bowtie_index_path} 
    -f 
    -U {fasta_file} 
    --very-sensitive 
    -p 32 > {output_bowtie_sam}
    ```
    """

    def run(self):
        """
           Run Assembly
        """
        input_fasta = self.input_files_local[0][-1]
        bowtie_fasta = self.input_files_local[0][-1]
        input_fasta2 = None
        if len(self.input_files_local[0]) >= 2:
            input_fasta = self.input_files_local[0][0]
            input_fasta2 = self.input_files_local[0][1]

        assembled_contig, assembled_scaffold, bowtie_sam, contig_stats = self.output_files_local()
        read2contig = {}
        memory = self.additional_attributes.get('memory', 100)
        self.assemble(input_fasta, input_fasta2, bowtie_fasta, assembled_contig, assembled_scaffold,
                      bowtie_sam, contig_stats, read2contig, int(memory))
    @staticmethod
    def assemble(input_fasta,
                 input_fasta2,
                 bowtie_fasta, # fasta file for running bowtie against contigs
                 assembled_contig,
                 assembled_scaffold,
                 bowtie_sam,
                 contig_stats,
                 read2contig,
                 memory=100):
        basedir = os.path.dirname(assembled_contig)
        assembled_dir = os.path.join(basedir, 'spades')
        command.execute(f"mkdir -p {assembled_dir}")
        assembled_contig_tmp = os.path.join(assembled_dir, 'contigs.fasta')
        assembled_scaffold_tmp = os.path.join(assembled_dir, 'scaffolds.fasta')

        try:
            if input_fasta2:
                command.execute(f"spades.py -1 {input_fasta} -2 {input_fasta2} -o {assembled_dir} -m {memory} -t 32 --only-assembler")
            else:
                command.execute(f"spades.py -s {input_fasta} -o {assembled_dir} -m {memory} -t 32 --only-assembler")
            command.execute(f"mv {assembled_contig_tmp} {assembled_contig}")
            command.execute(f"mv {assembled_scaffold_tmp} {assembled_scaffold}")

            # build the bowtie index based on the contigs
            PipelineStepRunAssembly.generate_read_to_contig_mapping(assembled_contig, bowtie_fasta,
                                                                    read2contig, bowtie_sam, contig_stats)
        except:
            # Assembly failed. create dummy output files
            command.execute(f"echo ';ASSEMBLY FAILED' > {assembled_contig}")
            command.execute(f"echo ';ASSEMBLY FAILED' > {assembled_scaffold}")
            command.execute(f"echo '@NO INFO' > {bowtie_sam}")
            command.execute("echo '{}' > " +  contig_stats)
            traceback.print_exc()

        command.execute(f"rm -rf {assembled_dir}")

    @staticmethod
    def generate_read_to_contig_mapping(assembled_contig,
                                        fasta_file,
                                        read2contig,
                                        output_bowtie_sam,
                                        output_contig_stats):
        ''' read -> contig mapping through bowtie2 alignment '''
        base_output_dir = os.path.dirname(fasta_file)
        # build bowtie index based on assembled_contig
        bowtie_index_path = os.path.join(base_output_dir, 'bowtie-contig')
        command.execute(f"mkdir -p {bowtie_index_path}; bowtie2-build {assembled_contig} {bowtie_index_path}")
        command.execute(f"bowtie2 -x {bowtie_index_path} -f -U {fasta_file} --very-sensitive -p 32 > {output_bowtie_sam}")
        contig_stats = defaultdict(int)
        PipelineStepRunAssembly.generate_info_from_sam(output_bowtie_sam, read2contig, contig_stats)
        with open(output_contig_stats, 'w') as ocf:
            json.dump(contig_stats, ocf)

    @staticmethod
    def generate_info_from_sam(bowtie_sam_file, read2contig, contig_stats):
        with open(bowtie_sam_file, "r", encoding='utf-8') as samf:
            for line in samf:
                if line[0] == '@':
                    continue
                fields = line.split("\t")
                read = fields[0]
                contig = fields[2]
                contig_stats[contig] += 1
                if contig != '*':
                    read2contig[read] = contig



    def count_reads(self):
        ''' count reads '''
        pass
