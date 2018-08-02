''' Generate Phylogenetic tree '''
import os
import json
import shelve

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.steps.generate_alignment_viz import PipelineStepGenerateAlignmentViz
import idseq_dag.util.command as command
import idseq_dag.util.s3 as s3

class PipelineStepGeneratePhyloTree(PipelineStep):
    ''' 
    Generate a phylogenetic tree from the input fasta files using kSNP3:
    http://gensoft.pasteur.fr/docs/kSNP3/01/kSNP3.01%20User%20Guide%20.pdf
    Augment the inputs with
      (a) NCBI sequences for the accession IDs specified in align_viz_json
    and/or
      (b) Genbank full reference genomes from the same taxid.
    '''
    def run(self):
        input_files = self.input_files_local[0]
        output_files = self.output_files_local()
        taxid = self.additional_attributes["taxid"]

        # knsp3 has a command (MakeKSNP3infile) for making a ksnp3-compatible input file from a directory of fasta files.
        # Before we can use the command, we copy/symlink all fasta files to a dedicated directory.
        # The command makes certain unreasonable assumptions we'll need to enfore:
        # - current directory is parent directory of the fasta file directory
        # - file names do not have dots except before extension (also no spaces).
        input_dir_for_ksnp3 = f"{self.output_dir_local}/inputs_for_ksnp3"
        command.execute(f"mkdir {input_dir_for_ksnp3}")
        for local_file in input_files:
            command.execute(f"ln -s {local_file} {input_dir_for_ksnp3}/{os.path.basename(local_file)}")

        # Retrieve Genbank references (full assembled genomes).
        # For now, we skip this using the option n=0 because
        # (a) sequences for the accession IDs actually matched by the sample are likely to be more relevant initially
        # (b) the downloads are slow.
        _local_genbank_fastas = self.get_genbank_genomes(taxid, input_dir_for_ksnp3, 0)

        # Retrieve NCBI NT references for the accessions in the alignment viz files.
        # These are the accessions (not necessarily full genomes) that were actually matched
        # by the sample's reads during GSNAP alignment.
        _local_accession_fastas = self.get_accession_sequences(input_dir_for_ksnp3)

        # Run MakeKSNP3infile.
        command.execute(f"cd {input_dir_for_ksnp3}/..; MakeKSNP3infile {os.path.basename(input_dir_for_ksnp3)} {self.output_dir_local}/inputs.txt A")

        # Now run ksnp3.
        command.execute(f"cd {self.output_dir_local}; mkdir ksnp3_outputs; kSNP3 -in inputs.txt -outdir ksnp3_outputs -k 13")
        command.execute(f"mv {self.output_dir_local}/ksnp3_outputs/tree.parsimony.tre {output_files[0]}")

    def count_reads(self):
        pass

    def get_genbank_genomes(self, taxid, destination_dir, n=10):
        '''
        Retrieve up to n GenBank reference genomes under taxid.
        Assumes taxid is species-level.
        Saves the references under file names compatible with MakeKSNP3infile.
        '''
        if n == 0:
            return []
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
            cmd += f" | head -n {n} | cut -f1,3,4" # take only top n results, keep name and ftp_path
            genomes = list(filter(None, command.execute_with_output(cmd).split("\n")))
            command.execute_with_output(f"rm {genome_list_local}")
            if genomes:
                local_genbank_fastas = []
                for line in genomes:
                    taxid, organism_name, ftp_path = line.split("\t")
                    clean_organism_name = PipelineStepGeneratePhyloTree.clean_name_for_ksnp3(organism_name)
                    ftp_fasta_gz = f"{ftp_path}/{os.path.basename(ftp_path)}_genomic.fna.gz"
                    local_fasta = f"{destination_dir}/genbank__{clean_organism_name}__taxid-{taxid}.fasta"
                    if os.path.isfile(local_fasta):
                        local_fasta = f"{local_fasta.split('.')[0]}__I.fasta"
                    command.execute(f"wget -O {local_fasta}.gz {ftp_fasta_gz}")
                    command.execute(f"gunzip {local_fasta}.gz")
                    local_genbank_fastas.append(local_fasta)
                return local_genbank_fastas
        return []

    def get_accession_sequences(self, dest_dir):
        '''
        Retrieve NCBI NT references for the accessions in the alignment viz files
        and write them to separate fasta files.
        '''
        # Retrieve files
        nt_db = self.additional_attributes["nt_db"]
        nt_loc_db = s3.fetch_from_s3(
            self.additional_files["nt_loc_db"],
            self.ref_dir_local,
            allow_s3mi=True)
        align_viz_s3_files = [key for key in list(self.additional_files.keys()) if key != "nt_loc_db"]
        local_align_viz_files = []
        for s3_file in align_viz_s3_files:
            local_file = s3.fetch_from_s3(
                self.additional_files[s3_file],
                self.ref_dir_local)
            local_align_viz_files.append(local_file)

        # Make map of accession to sequence
        accessions = set()
        for local_file in local_align_viz_files:
            with open(local_file, 'rb') as f:
                align_viz_dict = json.load(f)
            accessions |= set(align_viz_dict.keys())
        accession2info = dict((acc, {}) for acc in accessions)
        nt_loc_dict = shelve.open(nt_loc_db.replace(".db", ""))
        PipelineStepGenerateAlignmentViz.get_sequences_by_accession_list_from_s3(
            accession2info, nt_loc_dict, nt_db)

        # Write 1 fasta file per accession
        local_accession_fastas = []
        for acc, info in accession2info.items():
            clean_accession = PipelineStepGeneratePhyloTree.clean_name_for_ksnp3(acc)
            local_fasta = f"{dest_dir}/NCBI_NT_accession_{clean_accession}"
            with open(local_fasta, 'wb'):
                local_fasta.write(f">{acc}\n")
                local_fasta.write(info['ref_seq'])
            local_accession_fastas += local_fasta

        # Return paths of the new fasta files
        return local_accession_fastas

    @staticmethod
    def clean_name_for_ksnp3(name):
        return name.replace(' ', '-').replace('.', '')
