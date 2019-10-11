import json
import os
import multiprocessing
from collections import defaultdict
from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.trace_lock import TraceLock
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.s3 as s3
import idseq_dag.util.m8 as m8
import idseq_dag.util.log as log

from idseq_dag.steps.run_assembly import PipelineStepRunAssembly

MIN_REF_FASTA_SIZE = 25
MIN_ASSEMBLED_CONTIG_SIZE = 25

# When composing a query cover from non-overlapping fragments, consider fragments
# that overlap less than this fraction to be disjoint.
NT_MIN_OVERLAP_FRACTION = 0.1

# NT alginments with shorter length are associated with a high rate of false positives.
# NR doesn't have this problem because Rapsearch contains an equivalent filter.
# TODO: Nevertheless, it may be useful to re-filter blastx results.
NT_MIN_ALIGNMENT_LEN = 36

# Ignore NT local alignments (in blastn) with sequence similarity below 80%.
#
# Considerations:
#
#   This should be no higher than the equivalent threshold for GSNAP. That way, if
#   BLAST finds a better alignment for the contig than GSNAP found for some of the
#   reads, the BLAST result wouldn't be discarded in favor of an inferior GSNAP
#   result for those reads.
#
#   Conversely, this threshold should be pretty high, with experts agreeing that
#   a match with lesser than 80% quality here should likely be left for NR/blastx
#   to identify.
#
NT_MIN_PIDENT = 80


def intervals_overlap(p, q):
    '''Return True iff the intersection of p and q covers more than NT_MIN_OVERLAP_FRACTION of either p or q.'''
    p_len = p[1] - p[0]
    q_len = q[1] - q[0]
    shorter_len = min(p_len, q_len)
    intersection_len = max(0, min(p[1], q[1]) - max(p[0], q[0]))
    return intersection_len > (NT_MIN_OVERLAP_FRACTION * shorter_len)


# Note:  HSP is a BLAST term that stands for "highest-scoring pair", i.e., a local
# alignment with no gaps that achieves one of the highest alignment scores
# in a given search.  Each HSP spans an interval in the query sequence as well
# as an interval in the subject sequence.


def query_interval(row):
    # decode HSP query interval
    # depending on strand, qstart and qend may be reversed
    q = (row["qstart"], row["qend"])
    return min(*q), max(*q)


def hsp_overlap(hsp_1, hsp_2):
    # TODO:  Should we take into account subject_interval overlap?
    return intervals_overlap(query_interval(hsp_1), query_interval(hsp_2))


def intersects(needle, haystack):
    ''' Return True iff needle intersects haystack.  Ignore overlap < NT_MIN_OVERLAP_FRACTION. '''
    return any(hsp_overlap(needle, hay) for hay in haystack)


class Candidate:
    '''A list of blast highest-scoring pairs (HSPs) from the same query to the same subject sequence (candidate), ordered by decreasing bitscore.'''

    def __init__(self, hsps):
        self.hsps = hsps
        self.optimal_cover = None
        # agscores are non-negative
        self.agscore = None

    def find_optimal_fragment_cover(self):
        # Find a subset of disjoint HSPs with maximum sum of bitscores.
        # Initial implementation:  Super greedy.  Takes advantage of the fact
        # that blast results are sorted by bitscore, highest first.
        self.optimal_cover = [self.hsps[0]]
        for next_hsp in self.hsps[1:]:
            if not intersects(next_hsp, self.optimal_cover):
                optimal_cover.append(next_hsp)
        self.agscore = sum(hsp["pident"] * hsp["length"] for hsp in self.optimal_cover)

    def summary_row(self):
        '''Optimal cover stats are used later in generate_coverage_viz.'''
        r = dict(self.optimal_cover[0])
        # aggregate across the optimal cover's HSPs
        r["length"] = sum(hsp["length"] for hsp in self.optimal_cover)
        if r["length"] == 0:
            # it's never zero, but...
            r["length"] = 1
        r["pident"] = sum(hsp["pident"] * hsp["length"] for hsp in self.optimal_cover) / r["length"]
        r["bitscore"] = sum(hsp["bitscore"] for hsp in self.optimal_cover)
        r["qstart"] = min(hsp["qstart"] for hsp in self.optimal_cover)
        r["qend"] = max(hsp["qend"] for hsp in self.optimal_cover)
        r["sstart"] = min(hsp["sstart"] for hsp in self.optimal_cover)
        r["send"] = max(hsp["send"] for hsp in self.optimal_cover)
        # add these two
        r["qcov"] = r["length"] / r["qlen"]
        r["hsp_count"] = len(self.optimal_cover)
        return r


class PipelineStepBlastContigs(PipelineStep):
    """ The BLAST step is run independently for the contigs. First against the NT-BLAST database
    constructed from putative taxa identified from short read alignments to NCBI NT using GSNAP.
    Then, against the NR-BLAST database constructed from putative taxa identified from short read
    alignments to NCBI NR with Rapsearch2.

    For NT:
    ```
    blast_command
    -query {assembled_contig}
    -db {blast_index_path}
    -out {blast_m8}
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'
    -evalue 1e-10
    -max_target_seqs 5000
    -num_threads 16
    ```

    For NR:

    ```
    blast_command
    -query {assembled_contig}
    -db {blast_index_path}
    -out {blast_m8}
    -outfmt 6
    -num_alignments 5
    -num_threads 16
    ```
    """

    # Opening lineage sqlite in parallel for NT and NR sometimes hangs.
    # In fact it's generally not helpful to run NT and NR in parallel in this step,
    # so we may broaden the scope of this mutex.
    cya_lock = multiprocessing.RLock()

    def run(self):
        '''
            1. summarize hits
            2. built blast index
            3. blast assembled contigs to the index
            4. update the summary
        '''
        (_align_m8, deduped_m8, hit_summary, orig_counts) = self.input_files_local[0]
        assembled_contig, _assembled_scaffold, bowtie_sam, _contig_stats = self.input_files_local[1]
        reference_fasta = self.input_files_local[2][0]

        (blast_m8, refined_m8, refined_hit_summary, refined_counts, contig_summary_json, blast_top_m8) = self.output_files_local()
        db_type = self.additional_attributes["db_type"]
        if os.path.getsize(assembled_contig) < MIN_ASSEMBLED_CONTIG_SIZE or \
           os.path.getsize(reference_fasta) < MIN_REF_FASTA_SIZE:
            # No assembled results or refseq fasta available.
            # Create empty output files.
            command.write_text_to_file(' ', blast_m8)
            command.write_text_to_file(' ', blast_top_m8)
            command.copy_file(deduped_m8, refined_m8)
            command.copy_file(hit_summary, refined_hit_summary)
            command.copy_file(orig_counts, refined_counts)
            command.write_text_to_file('[]', contig_summary_json)
            return #FIXME return in the middle of the function

        (read_dict, accession_dict, _selected_genera) = m8.summarize_hits(hit_summary)
        PipelineStepBlastContigs.run_blast(db_type, blast_m8, assembled_contig, reference_fasta, blast_top_m8)
        read2contig = {}
        contig_stats = defaultdict(int)
        PipelineStepRunAssembly.generate_info_from_sam(bowtie_sam, read2contig, contig_stats)

        (updated_read_dict, read2blastm8, contig2lineage, added_reads) = self.update_read_dict(
            read2contig, blast_top_m8, read_dict, accession_dict, db_type)
        self.generate_m8_and_hit_summary(updated_read_dict, added_reads, read2blastm8,
                                         hit_summary, deduped_m8,
                                         refined_hit_summary, refined_m8)

        # Generating taxon counts based on updated results
        lineage_db = s3.fetch_from_s3(
            self.additional_files["lineage_db"],
            self.ref_dir_local,
            allow_s3mi=False) # Too small to waste s3mi
        deuterostome_db = None
        evalue_type = 'raw'
        if self.additional_files.get("deuterostome_db"):
            deuterostome_db = s3.fetch_from_s3(self.additional_files["deuterostome_db"],
                                               self.ref_dir_local, allow_s3mi=False) # Too small for s3mi
        with TraceLock("PipelineStepBlastContigs-CYA", PipelineStepBlastContigs.cya_lock, debug=False):
            with log.log_context("PipelineStepBlastContigs", {"substep": "generate_taxon_count_json_from_m8", "db_type": db_type, "refined_counts": refined_counts}):
                m8.generate_taxon_count_json_from_m8(refined_m8, refined_hit_summary,
                                                     evalue_type, db_type.upper(),
                                                     lineage_db, deuterostome_db, refined_counts)

        # generate contig stats at genus/species level
        with log.log_context("PipelineStepBlastContigs", {"substep": "generate_taxon_summary"}):
            contig_taxon_summary = self.generate_taxon_summary(read2contig, contig2lineage, updated_read_dict, added_reads, db_type)

        with log.log_context("PipelineStepBlastContigs", {"substep": "generate_taxon_summary_json", "contig_summary_json": contig_summary_json}):
            with open(contig_summary_json, 'w') as contig_outf:
                json.dump(contig_taxon_summary, contig_outf)

        # Upload additional file
        contig2lineage_json = os.path.join(os.path.dirname(contig_summary_json), f"contig2lineage.{db_type}.json")
        with log.log_context("PipelineStepBlastContigs", {"substep": "contig2lineage_json", "contig2lineage_json": contig2lineage_json}):
            with open(contig2lineage_json, 'w') as c2lf:
                json.dump(contig2lineage, c2lf)

        self.additional_files_to_upload.append(contig2lineage_json)

    @staticmethod
    def run_blast(db_type, blast_m8, *args):
        blast_index_path = os.path.join(os.path.dirname(blast_m8), f"{db_type}_blastindex")
        if db_type == 'nt':
            PipelineStepBlastContigs.run_blast_nt(blast_index_path, blast_m8, *args)
        else:
            assert db_type == 'nr'
            PipelineStepBlastContigs.run_blast_nr(blast_index_path, blast_m8, *args)

    @staticmethod
    def generate_taxon_summary(read2contig, contig2lineage, read_dict, added_reads_dict, db_type):
        # Return an array with
        # { taxid: , tax_level:, contig_counts: { 'contig_name': <count>, .... } }
        genus_summary = defaultdict(lambda: defaultdict(int))
        species_summary = defaultdict(lambda: defaultdict(int))
        for read_id, read_info in read_dict.items():
            contig = read2contig.get(read_id)
            if contig and contig != '*' and contig2lineage.get(contig):
                # It's possible a contig doesn't blast to anything
                species_taxid, genus_taxid, _family_taxid = contig2lineage[contig]
                species_summary[species_taxid][contig] += 1
                genus_summary[genus_taxid][contig] += 1
            else:
                # not mapping to a contig
                species_taxid, genus_taxid = read_info[4:6]
                species_summary[species_taxid]['*'] += 1
                genus_summary[genus_taxid]['*'] += 1

        for read_id, read_info in added_reads_dict.items():
            contig = read2contig[read_id]
            species_taxid, genus_taxid, _family_taxid = contig2lineage[contig]
            species_summary[species_taxid][contig] += 1
            genus_summary[genus_taxid][contig] += 1
        # constructing the array for output
        output_array = []
        for idx, summary in enumerate([species_summary, genus_summary]):
            tax_level = idx + 1
            for taxid, contig_counts in summary.items():
                entry = {'taxid': taxid, 'tax_level': tax_level,
                         'count_type': db_type.upper(), 'contig_counts': contig_counts}
                output_array.append(entry)
        return output_array

    @staticmethod
    def generate_m8_and_hit_summary(consolidated_dict, added_reads, read2blastm8,
                                    hit_summary, deduped_m8,
                                    refined_hit_summary, refined_m8):
        ''' generate new m8 and hit_summary based on consolidated_dict and read2blastm8 '''
        # Generate new hit summary
        new_read_ids = added_reads.keys()
        with open(refined_hit_summary, 'w') as rhsf:
            with open(hit_summary, 'r', encoding='utf-8') as hsf:
                for line in hsf:
                    read_id = line.rstrip().split("\t")[0]
                    read = consolidated_dict[read_id]
                    output_str = "\t".join(read)
                    rhsf.write(output_str + "\n")
            # add the reads that are newly blasted
            for read_id in new_read_ids:
                read_info = added_reads[read_id]
                output_str = "\t".join(read_info)
                rhsf.write(output_str + "\n")
        # Generate new M8
        with open(refined_m8, 'w') as rmf:
            with open(deduped_m8, 'r', encoding='utf-8') as mf:
                for line in mf:
                    read_id = line.rstrip().split("\t")[0]
                    m8_line = read2blastm8.get(read_id)
                    if m8_line:
                        m8_fields = m8_line.split("\t")
                        m8_fields[0] = read_id
                        rmf.write("\t".join(m8_fields))
                    else:
                        rmf.write(line)
            # add the reads that are newly blasted
            for read_id in new_read_ids:
                m8_line = read2blastm8.get(read_id)
                m8_fields = m8_line.split("\t")
                m8_fields[0] = read_id
                rmf.write("\t".join(m8_fields))

    @staticmethod
    def update_read_dict(read2contig, blast_top_m8, read_dict, accession_dict, db_type):
        consolidated_dict = read_dict
        read2blastm8 = {}
        contig2accession = {}
        contig2lineage = {}
        added_reads = {}

        for row, raw_line in m8.parse_tsv(blast_top_m8, m8.RERANKED_BLAST_OUTPUT_SCHEMA[db_type], raw_lines=True):
            contig_id = row["qseqid"]
            accession_id = row["sseqid"]
            contig2accession[contig_id] = (accession_id, raw_line)
            contig2lineage[contig_id] = accession_dict[accession_id]

        for read_id, contig_id in read2contig.items():
            (accession, m8_line) = contig2accession.get(contig_id, (None, None))
            if accession:
                (species_taxid, genus_taxid, family_taxid) = accession_dict[accession]
                if consolidated_dict.get(read_id):
                    consolidated_dict[read_id] += [contig_id, accession, species_taxid, genus_taxid, family_taxid]
                    consolidated_dict[read_id][2] = species_taxid
                else:
                    added_reads[read_id] = [read_id, "1", species_taxid, accession, species_taxid, genus_taxid, family_taxid, contig_id, accession, species_taxid, genus_taxid, family_taxid, 'from_assembly']
            if m8_line:
                read2blastm8[read_id] = m8_line
        return (consolidated_dict, read2blastm8, contig2lineage, added_reads)


    @staticmethod
    def run_blast_nt(blast_index_path, blast_m8, assembled_contig, reference_fasta, blast_top_m8):
        blast_type = 'nucl'
        blast_command = 'blastn'
        min_alignment_length = NT_MIN_ALIGNMENT_LEN
        min_pident = NT_MIN_PIDENT
        command.execute(
            command_patterns.SingleCommand(
                cmd="makeblastdb",
                args=[
                    "-in",
                    reference_fasta,
                    "-dbtype",
                    blast_type,
                    "-out",
                    blast_index_path
                ]
            )
        )
        command.execute(
            command_patterns.SingleCommand(
                cmd=blast_command,
                args=[
                    "-query",
                    assembled_contig,
                    "-db",
                    blast_index_path,
                    "-out",
                    blast_m8,
                    "-outfmt",
                    '6 ' + ' '.join(m8.BLAST_OUTPUT_NT_SCHEMA.keys()),
                    '-evalue',
                    1e-10,
                    '-max_target_seqs',
                    5000,
                    "-num_threads",
                    16
                ]
            )
        )
        # further processing of getting the top m8 entry for each contig.
        PipelineStepBlastContigs.get_top_m8_nt(blast_m8, blast_top_m8, min_alignment_length, min_pident)


    @staticmethod
    def run_blast_nr(blast_index_path, blast_m8, assembled_contig, reference_fasta, blast_top_m8):
        blast_type = 'prot'
        blast_command = 'blastx'
        command.execute(
            command_patterns.SingleCommand(
                cmd="makeblastdb",
                args=[
                    "-in",
                    reference_fasta,
                    "-dbtype",
                    blast_type,
                    "-out",
                    blast_index_path
                ]
            )
        )
        command.execute(
            command_patterns.SingleCommand(
                cmd=blast_command,
                args=[
                    "-query",
                    assembled_contig,
                    "-db",
                    blast_index_path,
                    "-out",
                    blast_m8,
                    "-outfmt",
                    6,
                    "-num_alignments",
                    5,
                    "-num_threads",
                    16
                ]
            )
        )
        # further processing of getting the top m8 entry for each contig.
        PipelineStepBlastContigs.get_top_m8_nr(blast_m8, blast_top_m8)


    @staticmethod
    def get_top_m8_nr(orig_m8, blast_top_m8):
        ''' Get top m8 file entry for each read from orig_m8 and output to blast_top_m8 '''
        with open(blast_top_m8, 'w') as top_m8f:
            top_line = None
            top_bitscore = 0
            current_read_id = None
            for read_id, _accession_id, _percent_id, _alignment_length, _e_value, bitscore, line in m8.iterate_m8(orig_m8):
                # Get the top entry of each read_id based on the bitscore
                if read_id != current_read_id:
                    # Different batch start
                    if current_read_id: # Not the first line
                        top_m8f.write(top_line)
                    current_read_id = read_id
                    top_line = line
                    top_bitscore = bitscore
                elif bitscore > top_bitscore:
                    top_bitscore = bitscore
                    top_line = line
            if top_line is not None:
                top_m8f.write(top_line) # TODO: Unify reranked formats for NT and NR


    @staticmethod
    def get_top_m8_nt(blast_output, blast_top_m8, min_alignment_length, min_pident):
        '''
        For each contig Q (query) and reference S (subject), extend the highest-scoring
        fragment alignment of Q to S with other non-overlapping fragments as far as
        possible, forming a set of fragments called HSPs(Q, S) that maximizes cumulative
        bitscore while avoiding overlap in Q.

        For each Q, output the S with highest agscore(Q, S), which is defined as the
        number of matching base pairs in all fragments that belong to HSPs(Q, S).

        Note that agscore(Q, S) is NOT the sum of bitscores in HSPs(Q, S) because of
        concerns that cumulative bitscores might not be directly comparable across
        different reference sequences S.

        TODO:  Document and discuss these concerns and score choices.'''

        # HSP is a BLAST term that stands for "highest-scoring pair", i.e., a local
        # alignment with no gaps that achieves one of the highest alignment scores
        # in a given search.
        #
        # The output of BLAST consists of one HSP per line, where lines corresponding
        # to the same (query, subject) pair are ordered by decreasing bitscore.
        #
        # See http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
        # for documentation of Blast output.

        # Group blast output HSPs by (query_id, subject_id).
        HSPs = defaultdict(list)
        for hsp in m8.parse_tsv(blast_output, m8.BLAST_OUTPUT_NT_SCHEMA):
            # filter local alignment HSPs based on minimum length and sequence similarity
            if hsp["length"] < min_alignment_length:
                continue
            if hsp["pident"] < min_pident:
                continue
            candidate_id = (hsp["qseqid"], hsp["sseqid"])
            HSPs[candidate_id].append(hsp)

        # Identify each query's optimal hit through ranking candidate hits by agscore,
        # where a candidate hit's agscore is determined by its optimal fragment cover.
        argmax = {}
        for candidate_id, fragments in HSPs.items():
            candidate = Candidate(fragments)
            candidate.find_optimal_fragment_cover()
            (query_id, _) = candidate_id
            if query_id not in argmax or argmax[query_id].agscore < candidate.agscore:
                argmax[query_id] = candidate

        # Output the optimal hit for each query.
        m8.unparse_tsv(blast_top_m8, m8.RERANKED_BLAST_OUTPUT_NT_SCHEMA,
                       (optimal_hit.summary_row() for optimal_hit in argmax.values()))
