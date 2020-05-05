import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.fasta as fasta
import idseq_dag.util.m8 as m8

from idseq_dag.engine.pipeline_step import PipelineCountingStep
from idseq_dag.util.cdhit_clusters import parse_clusters_file


class PipelineStepGenerateAnnotatedFasta(PipelineCountingStep):
    '''
    generate annotated fasta
    '''

    def _annotated_fasta(self):
        return self.output_files_local()[0]

    def _unidentified_fasta(self):
        return self.output_files_local()[1]

    def run(self):
        ''' annotate fasta '''
        merged_fasta = self.input_files_local[0][-1]
        gsnap_m8 = self.input_files_local[1][1]
        rapsearch2_m8 = self.input_files_local[2][1]

        if len(self.input_files_local) == 4:
            cdhitdup_clusters, deduped_fasta = self.input_files_local[3]
        else:
            cdhitdup_clusters, deduped_fasta = None, None

        annotated_fasta = self._annotated_fasta()
        unidentified_fasta = self._unidentified_fasta()
        self.annotate_fasta_with_accessions(merged_fasta, gsnap_m8, rapsearch2_m8, annotated_fasta)
        self.generate_unidentified_fasta(annotated_fasta, unidentified_fasta, clusters_dict)

    def count_reads(self):
        # The webapp expects this count to be called "unidentified_fasta"
        super()._count_reads_work(
            cluster_key=PipelineStepGenerateAnnotatedFasta.old_read_name,
            counter_name="unidentified_fasta",
            fasta_files=[self._unidentified_fasta()]
        )

    @staticmethod
    def annotate_fasta_with_accessions(merged_input_fasta, nt_m8, nr_m8, output_fasta):
        def get_map(m8_file):
            return dict((read_id, accession_id)
                        for read_id, accession_id, _percent_id,
                        _alignment_length, _e_value, _bitscore, _line in m8.iterate_m8(
                            m8_file, 0, "annotate_fasta_with_accessions"))

        nt_map = get_map(nt_m8)
        nr_map = get_map(nr_m8)

        with open(merged_input_fasta, 'r', encoding='utf-8') as input_fasta_f:
            with open(output_fasta, 'w') as output_fasta_f:
                sequence_name = input_fasta_f.readline()
                sequence_data = input_fasta_f.readline()

                while sequence_name and sequence_data:
                    read_id = sequence_name.rstrip().lstrip('>')
                    # Need to annotate NR then NT in this order for alignment viz
                    new_read_name = "NR:{nr_accession}:NT:{nt_accession}:{read_id}".format(  # Its inverse is old_read_name()
                        nr_accession=nr_map.get(read_id, ''),
                        nt_accession=nt_map.get(read_id, ''),
                        read_id=read_id)
                    output_fasta_f.write(">%s\n" % new_read_name)
                    output_fasta_f.write(sequence_data)
                    sequence_name = input_fasta_f.readline()
                    sequence_data = input_fasta_f.readline()

    @staticmethod
    def old_read_name(new_read_name):
        # Inverse of new_read_name creation above.  Needed to cross-reference to original read_id
        # in order to identify all duplicate reads for this read_id.
        return new_read_name.split(":", 4)[-1]

    def generate_unidentified_fasta(self, input_fa, output_fa, clusters_dict):
        """
        Generates files with all unmapped reads. If COUNT_ALL, which was added
        in v4, then include non-unique reads extracted upstream by cdhitdup.
        """
        if READ_COUNTING_MODE == ReadCountingMode.COUNT_UNIQUE:
            # This is the old way of generating unmapped reads, kept here for consistency.
            # TODO  remove annotated fasta intermediate file and replace > with : below
            command.execute(
                command_patterns.ShellScriptCommand(
                    script=r'''grep -A 1 '>NR::NT::' "$1" | sed '/^--$/d' > "$2";''',
                    args=[
                        input_fa,
                        output_fa
                    ]
                )
            )
            return

        self._generate_unidentified_fasta(input_fa)

    def _generate_unidentified_fasta(
        self,
        input_fa: str,
        output_fa: str,
        clusters_dict
    ):
        assert READ_COUNTING_MODE == ReadCountingMode.COUNT_ALL
        # NOTE: this will load the set of all original read headers, which
        # could be several GBs in the worst case.
        clusters_dict = parse_clusters_file(
            # TODO: (gdingle): verify filenames
            self.input_files_local[2][0],
            self.input_files_local[3][0]
        )

        with fasta.iterator(fasta_file) as input_file, \
                open(output_fa, "w") as output_file:
            for read in input_file:
                if not read.header.startswith('>NR::NT::'):  # has accession
                    continue

                output_fa.write(read.header + "\n")
                output_fa.write(read.sequence + "\n")

                other_headers = clusters_dict[header][1:]
                for other_header in other_headers:
                    output_fa.write(other_header + "\n")
