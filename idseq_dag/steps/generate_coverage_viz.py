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
        contigs_map = self.generate_contigs_map(blast_top_m8, valid_contigs_with_read_counts)
        reads_map = self.generate_reads_map(gsnap_deduped_m8, unassigned_reads_set)

        # Add coverage to the contig data.
        self.augment_contigs_map_with_coverage(contig_coverage_json, contigs_map)

        # Add byteranges to the contig data.
        self.augment_contigs_map_with_byteranges(contigs_fasta, contigs_map)

        # Before selecting the best accessions and removing the rest, get the total accession count for each taxon.
        taxons_to_total_accession_count = coverage_utils.get_taxons_to_total_accession_count(taxons_to_accessions)

        # Select the best accessions for each taxon.
        (taxons_to_accessions, accession_data) = coverage_utils.select_best_accessions_per_taxon(
            taxons_to_accessions, accession_data, contigs_map, reads_map, num_accessions_per_taxon
        )

        # For each accession, generate the JSON that will be sent to the coverage_viz.
        coverage_viz_obj = coverage_utils.generate_coverage_viz_json(accession_data, contigs_map, reads_map, max_num_bins_coverage)

        # Generate the summary JSON which will get to the report page. JSON contains a map of taxons to valid accessions.
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

            self.additional_files_to_upload.append(upload_file)


    def get_valid_contigs_with_read_counts(self, contig_stats_json, min_contig_size):
        valid_contigs_with_read_counts = {}
        with open(contig_stats_json, 'r') as csj:
            contig_read_counts = json.load(csj)

            for contig in contig_read_counts:
                if contig != "*" and contig_read_counts[contig] >= min_contig_size:
                    valid_contigs_with_read_counts[contig] = contig_read_counts[contig]

        return valid_contigs_with_read_counts


    def generate_accession_data(self, hit_summary, valid_contigs_with_read_counts):
        accession_data = defaultdict(lambda: {'reads': [], 'contigs': set() })
        # We only allow species taxons for now.
        taxons_to_accessions = defaultdict(set)

        line_count = 0
        with open(hit_summary, 'r') as hs:
            for line in hs:
                line_count += 1
                if line_count % 100000 == 0:
                    log.write(f"{line_count} lines in the hitsummary file processed.")

                values = line.rstrip().split("\t")

                # Only add contig if the contig is "valid", i.e. it has 4 or more reads.
                if len(values) == 12 and values[7] in valid_contigs_with_read_counts:
                    taxons_to_accessions[values[9]].add(values[8])
                    accession_data[values[8]]["contigs"].add(values[7])
                else:
                    taxons_to_accessions[values[4]].add(values[3])
                    accession_data[values[3]]["reads"].append(values[0])

        for accession_id, data in accession_data.items():
            accession_data[accession_id]["contigs"] = list(data["contigs"])

        return (accession_data, taxons_to_accessions)


    def augment_accession_data_with_info(self, info_dict, accession_data):
        for accession_id in accession_data:
            entry = info_dict.get(accession_id)

            if entry:
                accession_data[accession_id]["name"] = entry[0]
                accession_data[accession_id]["total_length"] = int(entry[1])
            else:
                accession_data[accession_id]["name"] = "Unknown accession"
                accession_data[accession_id]["total_length"] = 0


    def generate_contigs_map(self, blast_top_m8, valid_contigs_with_read_counts):
        contigs = {}

        # File is empty.
        if os.path.getsize(blast_top_m8) < MIN_M8_FILE_SIZE:
            return contigs

        # iterate_m8 automatically removes invalid hits.
        for _contig_id, _accession_id, _percent_id, _alignment_length, _e_value, _bitscore, line in m8.iterate_m8(blast_top_m8):
            parts = line.split("\t")
            name_parts = parts[0].split("_")

            if parts[0] in valid_contigs_with_read_counts:
                contigs[parts[0]] = {
                    "total_length": int(name_parts[3]),
                    "accession": parts[1],
                    "query_start": int(parts[6]),
                    "query_end": int(parts[7]),
                    "subject_start": int(parts[8]),
                    "subject_end": int(parts[9]),
                    "prop_mismatch": int(parts[4]) / int(parts[3]),
                    "percent_id": float(parts[2]),
                    "alignment_length": int(parts[3]),
                    "num_mismatches": int(parts[4]),
                    "num_gaps": int(parts[5]),
                    "num_reads": valid_contigs_with_read_counts[parts[0]]
                }

        return contigs


    # We process gsnap.deduped.m8 instead of gsnap.reassigned.m8 because we ignore contigs with read_count < 4.
    # However, these contigs still get reassigned in gsnap.reassigned.m8,
    # and overwrite the original read alignment to the accession, which we want.
    def generate_reads_map(self, deduped_m8, unassigned_reads_set):
        reads = {}

        # File is empty.
        if os.path.getsize(deduped_m8) < MIN_M8_FILE_SIZE:
            return contigs

        # iterate_m8 automatically removes invalid hits.
        for _read_id, _accession_id, _percent_id, _alignment_length, _e_value, _bitscore, line in m8.iterate_m8(deduped_m8):
            parts = line.split("\t")

            if parts[0] in unassigned_reads_set:
                reads[parts[0]] = {
                    # Extract the length from the contig name.
                    "accession": parts[1],
                    "query_start": int(parts[6]),
                    "query_end": int(parts[7]),
                    "subject_start": int(parts[8]),
                    "subject_end": int(parts[9]),
                    "prop_mismatch": int(parts[4]) / int(parts[3]),
                    "percent_id": float(parts[2]),
                    "alignment_length": int(parts[3]),
                    "num_mismatches": int(parts[4]),
                    "num_gaps": int(parts[5])
                }

        return reads


    def augment_contigs_map_with_coverage(self, contig_coverage_json, contigs_map):
        with open(contig_coverage_json, 'r') as ccj:
            contig_coverage = json.load(ccj)

            for contig_name in contig_coverage:
                if contig_name in contigs_map:
                    contigs_map[contig_name]["coverage"] = contig_coverage[contig_name]["coverage"]


    def augment_contigs_map_with_byteranges(self, contigs_fasta, contigs_map):
        with open(contigs_fasta, 'r') as cf:
            seq_offset = 0
            seq_len = 0
            contig_name = ""

            for line in cf:
                if line[0] == '>':  # header line
                    if seq_len > 0 and contig_name in contigs_map:
                        contigs_map[contig_name]["byterange"] = [seq_offset, seq_len]

                    seq_offset = seq_offset + seq_len
                    seq_len = len(line)
                    contig_name = line[1:].strip()
                else:
                    seq_len += len(line)

            if seq_len > 0 and contig_name in contigs_map:
                contigs_map[contig_name]["byterange"] = [seq_offset, seq_len]
