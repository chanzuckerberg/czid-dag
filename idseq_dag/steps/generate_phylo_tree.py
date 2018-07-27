''' Generate Phylogenetic tree '''
import os

from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.log as log

class PipelineStepGeneratePhyloTree(PipelineStep):
    ''' 
    Generate a phylogenetic tree from the input fasta files using kSNP3:
    http://gensoft.pasteur.fr/docs/kSNP3/01/kSNP3.01%20User%20Guide%20.pdf
    Augment the inputs with a NCBI reference genomes from the same taxid.
    BELOW TO BE IMPLEMENTED
    '''
    def run(self):
        input_files = self.input_files_local[0]
        output_files = self.output_files_local()

        taxid = self.additional_attributes["taxid"]
        local_ncbi_fastas = self.get_ncbi_genomes(taxid)

        command.execute(f"echo {local_ncbi_fastas} > {output_files[0]}") # temporary

    def count_reads(self):
        pass

    def get_ncbi_genomes(self, taxid, n=10):
        '''
        Retrieve up to n GenBank reference genomes under taxid.
        Assumes taxid is species-level.
        '''
        categories = ["bacteria", "viral", "fungi", "protozoa"]
        # additional options in genbank that probably don't need right now:
        # ["archaea", "plant", 
        # "vertebrate_mammalian", "vertebrate_other", "invertebrate",
        # "other", "metagenomes"]
        for cat in categories:
            genome_list_path = f"ftp://ftp.ncbi.nih.gov/genomes/genbank/{cat}/assembly_summary.txt"
            genome_list_local = f"{self.output_dir_local}/{os.path.basename(genome_list_path)}"
            cmd = f"wget -O {genome_list_local} {genome_list_path}; "
            cmd += f"cut -f6,7,20 {genome_list_local}" # columns: 6 = taxid; 7 = species_taxid, 20 = ftp_path
            cmd += f" | grep '\\t{taxid}\\t'" # try to find taxid in the species_taxids
            cmd += f" | head -n {n} | cut -f1,3" # take only top n results
            genomes = list(filter(None, command.execute_with_output(cmd).split("\n")))
            if genomes:
                local_ncbi_fastas = []
                for line in genomes:
                    taxid, ftp_path = genomes.split("\t")
                    ftp_fasta_gz = f"{ftp_path}/{os.path.basename(ftp_path)}_genomic.fna.gz"
                    local_fasta = f"{self.output_dir_local}/refgenome_taxid_{taxid}.fasta"
                    command.execute(f"wget -O {local_fasta}.gz {ftp_fasta_gz}")
                    command.execute(f"gunzip {local_fasta}.gz")
                    local_ncbi_fastas.append(local_fasta)
                return local_ncbi_fastas
        return []
