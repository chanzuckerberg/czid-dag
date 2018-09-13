import json
from collections import defaultdict
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.count as count


MIN_READS_PER_GENUS = 100
MIN_PCT_PER_GENUS = 1.

class PipelineStepReclassifyReads(PipelineStep):
    '''
        Reclassify reads after alignment
    '''
    def run(self):
        '''
            1. organize the reads into different genus and get the relevant accessions from first pass
                1.1 decide based on # reads in a genus if we want to reclassify. min(100, 1% nonhost)  reads
                1.2 decide based on # species under a genus. > 1
            2a. download the actual reference sequences based on the accession data.
                build blastindex per genus
            2b. assemble  per genus at the same time from 1
            3. balst the contigs from 2b to 2a. get the one with the best bitscore
            4. build the bowtie index from 2b
            5. align the short reads from 1 to 4.
            6. for reads that align to the contigs, we assign the m8 entry
               based on where the best m8 entry the contig blast to and the species accordingly
            7. for reads that do not align to the contigs, we blast it and try to re-assign it

        '''
        (align_m8, deduped_m8, hit_summary, _orig_counts) = self.input_files_local[0]
        input_fasta = self.input_files_local[1][-1]
        (blast_m8, refined_m8, refined_hit_summary, refined_counts, genus_to_accession_list) =
        self.output_files_local()
        (read_dict, accession_dict, selected_genera) = self.summarize_hits(hit_summary)
        genus_fasta = self.group_reads_by_genus(input_fasta, read_dict, selected_genera)
        genus_references = self.download_ref_sequences(selected_genera)

        genus_assembled_contig = self.assemble_all(genus_fasta)
        genus_blast_m8 = self.blast_all(genus_assembled_contig, genus_references)
        genus_bowtie_index = self.bowtie_build_all(genus_assembled_contig)
        genus_bowtie_align_sam = self.bowtie_align_all(genus_bowtie_index, genus_fasta)
        consolidated_dict = self.consolidate_all_results(genus_bowtie_align_sam, genus_blast_m8, read_dict)



    def group_reads_by_genus(input_fasta, read_dict, selected_genera):
        genus_fasta = {}

        return genus_fasta


    def download_ref_sequences(self, accession_dict):


    @staticmethod
    def summarize_hits(hit_summary_file):

    @staticmethod
    def summarize_hits(hit_summary_file):
        read_dict = {}
        accession_dict = {}
        genus_read_counts = defaultdict(0)
        genus_species = defaultdict(set)
        genus_accessions = defaultdict(set)
        total_reads = 0
        with open(hit_summary_file, 'r') as hsf:
            for line in hsf:
                read = line.rstrip().split("\t")
                accession_id = read[3]
                species_taxid = read[4]
                genus_taxid = read[5]
                read_dict[read[0]] = read
                total_reads += 1
                if accession_id == 'None' or genus_taxid < 0:
                    continue
                accession_dict[accession_id] = (species_taxid, genus_taxid)
                genus_read_counts[genus_taxid] += 1
                genus_species[genus_taxid].add(species_taxid)
                genus_accessions[genus_taxid].add(accession_id)
        min_reads_per_genus = min(MIN_READS_PER_GENUS, total_reads/100*MIN_PCT_PER_GENUS)
        selected_genera = {}
        for genus_taxid, reads in genus_read_counts.items():
            if reads >= min_reads_per_genus and len(genus_species[genus_taxid]) > 1:
                selected_genera[genus_taxid] = list(genus_accessions[genus_taxid])

        return (read_dict, accession_dict, selected_genera)


    def count_reads(self):
        pass

