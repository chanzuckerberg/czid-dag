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

        # knsp3 has a command for making a ksnp3-compatible input file from a directory of fasta files.
        # So copy/symlink all fasta files to dedicated directory, then run that command (MakeKSNP3infile).
        # The command makes certain unreasonable assumptions:
        # - current directory is parent directory of the fasta file directory
        # - file names do not have dots except before extension (also no spaces), so we need to enforce that in the get_ncbi_genomes method.
        input_dir_for_ksnp3 = f"{self.output_dir_local}/inputs_for_ksnp3"
        command.execute(f"mkdir {input_dir_for_ksnp3}")
        for local_file in input_files:
            command.execute(f"ln -s {local_file} {input_dir_for_ksnp3}/{os.path.basename(local_file)}")
        local_ncbi_fastas = self.get_ncbi_genomes(taxid, input_dir_for_ksnp3)
        command.execute(f"cd {input_dir_for_ksnp3}/..; MakeKSNP3infile {os.path.basename(input_dir_for_ksnp3)} {self.output_dir_local}/inputs.txt A")

        # Now run ksnp3.
        command.execute(f"cd {self.output_dir_local}; mkdir ksnp3_outputs; kSNP3 -in inputs.txt -outdir ksnp3_outputs -k 13")
        command.execute(f"mv {self.output_dir_local}/ksnp3_outputs/tree.parsimony.tre {output_files[0]}")

    def count_reads(self):
        pass

    def get_ncbi_genomes(self, taxid, destination_dir, n=10):
        '''
        Retrieve up to n GenBank reference genomes under taxid.
        Assumes taxid is species-level.
        Saves the references under file names compatible with MakeKSNP3infile.
        '''
        categories = ["bacteria", "viral", "fungi", "protozoa"]
        # additional options in genbank that probably don't need right now:
        # ["archaea", "plant", 
        # "vertebrate_mammalian", "vertebrate_other", "invertebrate",
        # "other", "metagenomes"]
        for cat in categories:
            genome_list_path = f"ftp://ftp.ncbi.nih.gov/genomes/genbank/{cat}/assembly_summary.txt"
            genome_list_local = f"{destination_dir}/{os.path.basename(genome_list_path)}"
            cmd = f"wget -O {genome_list_local} {genome_list_path}; "
            cmd += f"cut -f6,7,8,20 {genome_list_local}" # columns: 6 = taxid; 7 = species_taxid, 8 = organism name, 20 = ftp_path
            cmd += f" | grep -P '\\t{taxid}\\t'" # try to find taxid in the species_taxids
            cmd += f" | head -n {n} | cut -f3,4" # take only top n results, keep name and ftp_path
            genomes = list(filter(None, command.execute_with_output(cmd).split("\n")))
            command.execute_with_output(f"rm {genome_list_local}")
            if genomes:
                local_ncbi_fastas = []
                for line in genomes:
                    organism_name, ftp_path = line.split("\t")
                    ftp_fasta_gz = f"{ftp_path}/{os.path.basename(ftp_path)}_genomic.fna.gz"
                    local_fasta = f"{destination_dir}/genbank__{organism_name.replace(' ', '-').replace('.', ''}.fasta"
                    command.execute(f"wget -O {local_fasta}.gz {ftp_fasta_gz}")
                    command.execute(f"gunzip {local_fasta}.gz")
                    local_ncbi_fastas.append(local_fasta)
                return local_ncbi_fastas
        return []
