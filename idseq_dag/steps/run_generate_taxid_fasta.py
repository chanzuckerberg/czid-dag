from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.taxid_lineage as taxid_lineage
import shelve


class PipelineStepRunGenerateTaxidFasta(PipelineStep):
    """Generate taxid FASTA from hit summaries. Intermediate conversion step
    that includes handling of non-specific hits with artificial tax_ids.
    """

    def run(self):
        input_fasta_file = ""
        hit_summary_files = {}
        lineage_path = ""
        output_fasta_file = ""

        lineage_map = shelve.open(lineage_path)
        valid_hits = PipelineStepRunGenerateTaxidFasta.parse_hits(hit_summary_files)

        input_fasta_f = open(input_fasta_file, 'rb')
        output_fasta_f = open(output_fasta_file, 'wb')
        sequence_name = input_fasta_f.readline()
        sequence_data = input_fasta_f.readline()
        while len(sequence_name) > 0 and len(sequence_data) > 0:
            # Example read_id: "NR::NT:CP010376.2:NB501961:14:HM7TLBGX2:1:23109
            # :12720:8743/2"
            # Translate the read information into our custom format with fake
            # taxids.
            accession_annotated_read_id = sequence_name.decode(
                "utf-8").rstrip().lstrip('>')
            read_id = accession_annotated_read_id.split(":", 4)[-1]

            nr_taxid_species, nr_taxid_genus, nr_taxid_family = PipelineStepRunGenerateTaxidFasta.get_valid_lineage(
                valid_hits, lineage_map, read_id, 'NR')
            nt_taxid_species, nt_taxid_genus, nt_taxid_family = PipelineStepRunGenerateTaxidFasta.get_valid_lineage(
                valid_hits, lineage_map, read_id, 'NT')

            family_str = 'family_nr:' + nr_taxid_family + ':family_nt:' + nt_taxid_family
            genus_str = ':genus_nr:' + nr_taxid_genus + ':genus_nt:' + nt_taxid_genus
            species_str = ':species_nr:' + nr_taxid_species + ':species_nt:' + nt_taxid_species
            new_read_name = (family_str + genus_str + species_str + ':' +
                             accession_annotated_read_id)

            output_fasta_f.write(">%s\n" % new_read_name)
            output_fasta_f.write(sequence_data)
            sequence_name = input_fasta_f.readline()
            sequence_data = input_fasta_f.readline()
        input_fasta_f.close()
        output_fasta_f.close()

    @staticmethod
    def parse_hits(hit_summary_files):
        """Return map of {NT, NR} => read_id => (hit_taxid_str, hit_level_str)"""
        valid_hits = {}
        for count_type, summary_file in hit_summary_files.items():
            hits = {}
            with open(summary_file, "rb") as sf:
                for hit_line in sf:
                    hit_line_columns = hit_line.decode("utf-8").strip().split(
                        "\t")
                    if len(hit_line_columns) >= 3:
                        hit_read_id = hit_line_columns[0]
                        hit_level_str = hit_line_columns[1]
                        hit_taxid_str = hit_line_columns[2]
                        hits[hit_read_id] = (hit_taxid_str, hit_level_str)
            valid_hits[count_type] = hits
        return valid_hits

    @staticmethod
    def get_valid_lineage(valid_hits, lineage_map, read_id, count_type):
        # If the read aligned to something, then it would be present in the
        # summary file for count type, and correspondingly in valid_hits[
        # count_type], even if the hits disagree so much that the
        # "valid_hits" entry is just ("-1", "-1"). If the read didn't align
        # to anything, we also represent that with ("-1", "-1"). This ("-1",
        # "-1") gets translated to NULL_LINEAGE.
        hit_taxid_str, hit_level_str = valid_hits[count_type].get(
            read_id, ("-1", "-1"))
        hit_lineage = lineage_map.get(hit_taxid_str,
                                      taxid_lineage.NULL_LINEAGE)
        return taxid_lineage.validate_taxid_lineage(hit_lineage, hit_taxid_str,
                                                    hit_level_str)
