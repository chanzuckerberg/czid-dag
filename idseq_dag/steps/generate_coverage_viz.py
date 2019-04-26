from collections import defaultdict
import json
import math
import os

import idseq_dag.util.m8 as m8
import idseq_dag.util.log as log
import idseq_dag.util.s3 as s3
import idseq_dag.util.command as command
from idseq_dag.engine.pipeline_step import PipelineStep
# coverage_utils contains the bulk of the business logic for this pipeline step
import idseq_dag.util.coverage as coverage_utils
from idseq_dag.util.dict import IdSeqDictValue, open_file_db_by_extension

MIN_M8_FILE_SIZE = 25 # M8 file should have at least a single line.

# These can be overridden with additional_attributes
MAX_NUM_BINS_COVERAGE = 500
NUM_ACCESSIONS_PER_TAXON = 10
MIN_CONTIG_SIZE = 4

class PipelineStepGenerateCoverageViz(PipelineStep):
    """Pipeline step to generate JSON file for read alignment visualizations to
    be consumed by the web app.
    """

    def run(self):
        (reassigned_m8, hit_summary, blast_top_m8) = self.input_files_local[0]
        (contig_coverage_json, contig_stats_json, contigs_fasta) = self.input_files_local[1]
        (gsnap_deduped_m8, ) = self.input_files_local[2]

        info_db = s3.fetch_from_s3(
            self.additional_files["info_db"],
            self.ref_dir_local,
            allow_s3mi=True)
        info_dict = open_file_db_by_extension(info_db, IdSeqDictValue.VALUE_TYPE_ARRAY)

        coverage_viz_summary = self.output_files_local()[0]

        max_num_bins_coverage = self.additional_attributes.get("max_num_bins_coverage", MAX_NUM_BINS_COVERAGE)
        num_accessions_per_taxon = self.additional_attributes.get("num_accessions_per_taxon", NUM_ACCESSIONS_PER_TAXON)
        min_contig_size = self.additional_attributes.get("min_contig_size", MIN_CONTIG_SIZE)

        # Get a map from valid contigs (based on MIN_CONTIG_SIZE) to number of reads in that contig.
        valid_contigs_with_read_counts = self.get_valid_contigs_with_read_counts(contig_stats_json, min_contig_size)

        # Get a map from accession name to hits (reads and contigs)
        # Also get a map from taxons to accession names.
        (accession_data, taxons_to_accessions) = self.generate_accession_data(hit_summary, valid_contigs_with_read_counts)
        # Remove taxons with no contigs in any of their accessions.
        coverage_utils.remove_taxons_with_no_contigs(accession_data, taxons_to_accessions)

        # Add the total_length and name of the accession to the accession data.
        self.augment_accession_data_with_info(info_dict, accession_data)

        # Get unassigned reads. Use a set for performance.
        unassigned_reads_set = coverage_utils.get_unassigned_reads_set(accession_data)

        # Extract information about contigs and reads.
        contig_data = self.generate_contig_data(blast_top_m8, valid_contigs_with_read_counts)
        read_data = self.generate_read_data(gsnap_deduped_m8, unassigned_reads_set)

        # Add coverage to the contig data.
        self.augment_contig_data_with_coverage(contig_coverage_json, contig_data)

        # Add byteranges to the contig data.
        self.augment_contig_data_with_byteranges(contigs_fasta, contig_data)

        # Before selecting the best accessions and removing the rest, get the total accession count for each taxon.
        taxons_to_total_accession_count = coverage_utils.get_taxons_to_total_accession_count(taxons_to_accessions)

        # Select the best accessions for each taxon.
        (taxons_to_accessions, accession_data) = coverage_utils.select_best_accessions_per_taxon(
            taxons_to_accessions, accession_data, contig_data, read_data, num_accessions_per_taxon
        )

        # For each accession, generate the JSON that will be sent to the coverage_viz.
        coverage_viz_obj = coverage_utils.generate_coverage_viz_json(accession_data, contig_data, read_data, max_num_bins_coverage)

        # Generate the summary JSON file which is initially loaded on the report page.
        coverage_viz_summary_json = coverage_utils.generate_coverage_viz_summary_json(taxons_to_accessions, accession_data, coverage_viz_obj, taxons_to_total_accession_count)

        # Write the summary JSON.
        with open(coverage_viz_summary, 'w') as cvs:
            json.dump(coverage_viz_summary_json, cvs)

        # Create a coverage viz JSON file for each accession
        coverage_viz_dir = os.path.join(self.output_dir_local, "coverage_viz")
        command.execute(f"mkdir -p {coverage_viz_dir}")
        for accession_id in coverage_viz_obj:
            upload_file = os.path.join(coverage_viz_dir, f"{accession_id}_coverage_viz.json")

            with open(upload_file, 'w') as uf:
                json.dump(coverage_viz_obj[accession_id], uf)

        self.additional_folders_to_upload.append(coverage_viz_dir)

    # Get a dict that maps valid contigs to their read count.
    # A contig is valid if it is larger than min_contig_size.
    @staticmethod
    def get_valid_contigs_with_read_counts(contig_stats_json, min_contig_size):
        valid_contigs_with_read_counts = {}
        with open(contig_stats_json, 'r') as csj:
            contig_read_counts = json.load(csj)

            for contig in contig_read_counts:
                if contig != "*" and contig_read_counts[contig] >= min_contig_size:
                    valid_contigs_with_read_counts[contig] = contig_read_counts[contig]

        return valid_contigs_with_read_counts

    # Generate a dict that maps accessions to the reads and contigs that were assigned to them.
    # Also generate a dict that maps taxons to accessions.
    @staticmethod
    def generate_accession_data(hit_summary, valid_contigs_with_read_counts):
        # Use a set for contigs, since multiple lines in the hitsummary can map to the same contig.
        accession_data = defaultdict(lambda: {'reads': [], 'contigs': set() })
        # Use a set for accessions, since multiple lines in the hitsummary can map to the same accession.
        taxons_to_accessions = defaultdict(set)

        line_count = 0
        with open(hit_summary, 'r') as hs:
            for line in hs:
                line_count += 1
                if line_count % 100000 == 0:
                    log.write(f"{line_count} lines in the hitsummary file processed.")

                values = line.rstrip().split("\t")

                # Only count the contig if the contig is valid, i.e. it is larger than min_contig_size.
                if len(values) == 12 and values[7] in valid_contigs_with_read_counts:
                    taxons_to_accessions[values[9]].add(values[8])
                    accession_data[values[8]]["contigs"].add(values[7])
                else:
                    taxons_to_accessions[values[4]].add(values[3])
                    accession_data[values[3]]["reads"].append(values[0])

        # Convert the contig set to a list.
        for accession_id, data in accession_data.items():
            accession_data[accession_id]["contigs"] = list(data["contigs"])

        return (accession_data, taxons_to_accessions)

    # Augment the accession data dictionary with accession data pulled from the info_db.
    @staticmethod
    def augment_accession_data_with_info(info_dict, accession_data):
        for accession_id in accession_data:
            entry = info_dict.get(accession_id)

            if entry:
                # The name of the accession.
                accession_data[accession_id]["name"] = entry[0]
                # The total length of the accession in base pairs.
                accession_data[accession_id]["total_length"] = int(entry[1])
            else:
                accession_data[accession_id]["name"] = "Unknown accession"
                accession_data[accession_id]["total_length"] = 0

    # Generate hit data from an m8 file.
    # Only include hits whose name appears in the valid_hits collection.
    @staticmethod
    def generate_hit_data_from_m8(m8_file, valid_hits):
        hits = {}

        # File is empty.
        if os.path.getsize(m8_file) < MIN_M8_FILE_SIZE:
            return hits

        # iterate_m8 automatically removes invalid hits.
        for (hit_id, accession_id, percent_id, alignment_length, num_mismatches, num_gaps,
            query_start, query_end, subject_start, subject_end, _e_value, _bitscore, line) in m8.iterate_m8(m8_file, full_line=True):

            if hit_id in valid_hits:
                # Map the hit_id to a dict of hit data.
                hits[hit_id] = {
                    "accession": accession_id,
                    "percent_id": percent_id,
                    "alignment_length": alignment_length,
                    "num_mismatches": num_mismatches,
                    "num_gaps": num_gaps,
                    "query_start": query_start,
                    "query_end": query_end,
                    "subject_start": subject_start,
                    "subject_end": subject_end,
                    "prop_mismatch": num_mismatches / alignment_length,
                }

        return hits

    # Generate contig data from blast_top_m8.
    @staticmethod
    def generate_contig_data(blast_top_m8, valid_contigs_with_read_counts):
        contigs = PipelineStepGenerateCoverageViz.generate_hit_data_from_m8(blast_top_m8, valid_contigs_with_read_counts)

        # Include some additional data.
        for contig_id, contig_obj in contigs.items():
            name_parts = contig_id.split("_")
            # Total length of the contig. We extract this from the contig name.
            contig_obj["total_length"] = int(name_parts[3])
            # The contig read count.
            contig_obj["num_reads"] = valid_contigs_with_read_counts[contig_id]

        return contigs

    # Generate read data from deduped_m8.
    # We process gsnap.deduped.m8 instead of gsnap.reassigned.m8 because we ignore contigs with read_count < 4.
    # However, these contigs still get reassigned in gsnap.reassigned.m8,
    # and overwrite the original read alignment to the accession, which we need. So we can't use gsnap.reassigned.m8.
    @staticmethod
    def generate_read_data(deduped_m8, unassigned_reads_set):
        return PipelineStepGenerateCoverageViz.generate_hit_data_from_m8(deduped_m8, unassigned_reads_set)

    # Augment contig data with the contig coverage array.
    @staticmethod
    def augment_contig_data_with_coverage(contig_coverage_json, contig_data):
        with open(contig_coverage_json, 'r') as ccj:
            contig_coverage = json.load(ccj)

            for contig_name in contig_coverage:
                if contig_name in contig_data:
                    contig_data[contig_name]["coverage"] = contig_coverage[contig_name]["coverage"]

    # Augment contig data with the byterange location in the contigs.fasta file for each contig.
    @staticmethod
    def augment_contig_data_with_byteranges(contigs_fasta, contig_data):
        with open(contigs_fasta, 'r') as cf:
            seq_offset = 0
            seq_len = 0
            contig_name = ""

            # Process each line in contigs.fasta.
            for line in cf:
                # If the line is a header file, process the contig we just traversed.
                if line[0] == '>':  # header line
                    if seq_len > 0 and contig_name in contig_data:
                        contig_data[contig_name]["byterange"] = [seq_offset, seq_len]

                    seq_offset = seq_offset + seq_len
                    seq_len = len(line)
                    contig_name = line[1:].strip()
                else:
                    seq_len += len(line)

            # Process the last contig once we reach the end of the file.
            if seq_len > 0 and contig_name in contig_data:
                contig_data[contig_name]["byterange"] = [seq_offset, seq_len]
