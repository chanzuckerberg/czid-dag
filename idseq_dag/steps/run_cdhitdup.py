import os

import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.count as count
import idseq_dag.util.log as log

from idseq_dag.engine.pipeline_step import InputFileErrors, PipelineStep


class PipelineStepRunCDHitDup(PipelineStep):  # Deliberately not PipelineCountingStep
    """ Removes duplicate reads.

    ```
    cd-hit-dup
    -i {input_fasta}
    -o {output_fasta}
    -e 0.05
    -u 70
    ```

    Per the CDHit Documentation, available [here](https://github.com/weizhongli/cdhit/wiki),
    this command uses a 0.05 threshold for the number of mismatches - indicating that if
    two reads are > 95% similar they will be considered duplicates. It uses only the
    first/last 70 bases of each read to do the analysis on sequence similarity.

    cd-hit-dup will output three or four files. Two of them are the same as the
    output files of CD-HIT: one (named exactly the same as the file name
    specified by the “-o” option) is the cluster (or duplicate) representatives,
    the other is the clustering file (xxx.clstr) relating each duplicate to its
    representative. For paired end reads, another file by the “-o2” option is
    the cluster representatives for R2 reads. The last file (xxx2.clstr)
    contains the chimeric clusters. In this file, the description for each
    chimeric cluster contains cluster ids of its parent clusters from the
    clustering file xxx.clstr.
    """

    def validate_input_files(self):
        if not count.files_have_min_reads(self.input_files_local[0], 2):
            self.input_file_error = InputFileErrors.INSUFFICIENT_READS

    def run(self):
        ''' Invoking cd-hit-dup '''
        input_fas = self.input_files_local[0]

        output_files = self.output_files_local()
        assert len(output_files) == len(input_fas) + 2, f"Context: {input_fas} -> {output_files}."
        output_fas = output_files[:len(input_fas)]
        assert output_files[-2].endswith(".clstr"), str(output_files)
        assert output_files[-1].endswith(".tsv"), str(output_files)
        cdhit_clusters_path = output_files[-2]
        cdhit_cluster_sizes_path = output_files[-1]

        # Require exact match (-e 0) on first 70 nucleotides (-u 70).
        # (Only 70 because sequencer errors increase toward the end of a read;
        # for example, see https://insilicoseq.readthedocs.io/en/latest/iss/model.html)
        #
        # Note that cdhit runs without -f here so it would not identify chimeric reads,
        # and therefore the second .clstr file output would be empty, always.
        cdhitdup_params = [
            '-i', input_fas[0], '-o', output_fas[0],
            '-e', '0.0', '-u', '70'
        ]
        if len(input_fas) == 2:
            cdhitdup_params += ['-i2', input_fas[1], '-o2', output_fas[1]]
        command.execute(
            command_patterns.SingleCommand(
                cmd='cd-hit-dup',
                args=cdhitdup_params
            )
        )
        self._add_clstr_files()

        # Emit cluster sizes.  One line per cluster.  Format "<cluster_size> <cluster_read_id>".
        # The rest of the pipeline operates on unique reads (i.e. on the representatives of these clusters).
        # In the end, cluster sizes are used to convert unique read counts to original read counts.
        with open(cdhit_cluster_sizes_path, "w") as tsv,\
             open(cdhit_clusters_path, "r") as clusters_file:  # noqa
            read_id = None
            cluster_size = 0
            for line in clusters_file:
                if line.startswith(">"):
                    continue
                # Example input lines that form a cluster:
                #
                #    "0       140nt, >M05295:357:000000000-CRPNR:1:2119:16143:8253... *"
                #    "1       140nt, >M05295:357:000000000-CRPNR:1:1101:22051:10534... at 1:140:1:140/+/100.00%"
                #    "2       140nt, >M05295:357:000000000-CRPNR:1:1102:15401:7483... at 1:140:1:140/+/100.00%"
                #    ...
                #    "2334    140nt, >M05295:357:000000000-CRPNR:1:1102:13405:3483... at 1:140:1:140/+/100.00%"
                #
                # Corresponding output line for that cluster:
                #
                #    "2335    140nt, >M05295:357:000000000-CRPNR:1:2119:16143:8253"
                #
                # Please note that "..." above does not indicate truncation. CD-HIT-DUP appens "..." to the read
                # IDs even if the read IDs have not been truncated.
                #
                # Per CD-HIT-DUP docs, when a "-d" argument is not specified, each read ID is obtained by
                # splitting the corresponding FASTA line on whitespace and taking the first item.  That is also
                # how we do it throughout this pipeline.  The code below assumes (and relies upon) the read IDs
                # being free from whitespace.
                #
                # Furthremore, we expect read IDs to be unique in a sequencing run.  Violating that assumption,
                # if it does not break cdhit itself, might produce slightly bogus numbers, because of how it
                # affects subsampling and per-read DCR correction.  Preserving this uniqueness is why we do not
                # specify a "-d" flag to cdhit, and allow it thus to use the entire read id.
                #
                parts = line.strip().split()
                serial = int(parts[0])
                if serial != 0:
                    assert cluster_size == serial
                    cluster_size += 1
                    continue
                if cluster_size != 0:
                    tsv.write(f"{cluster_size}\t{read_id}\n")
                # This is the first read for a new cluster.
                cluster_size = 1
                assert parts[2][0] == ">", line
                assert parts[2].endswith("..."), line
                assert parts[3] == "*", line
                read_id = parts[2][1:-3]
            if cluster_size != 0:
                tsv.write(f"{cluster_size}\t{read_id}\n")

    def count_reads(self):
        self.should_count_reads = True
        # Here we intentionally count unique reads.
        self.counts_dict[self.name] = count.reads_in_group(self.output_files_local()[0:2])

    def _add_clstr_files(self):
        output_fas = self.output_files_local()[0]
        clstr_file = output_fas + '.clstr'  # clusters
        clstr_file2 = output_fas + '2.clstr'  # chimeric clusters
        if os.path.isfile(clstr_file) and os.path.isfile(clstr_file2):
            self.additional_output_files_visible += [clstr_file, clstr_file2]
        else:
            log.write(
                f"WARNING: Files not found: {clstr_file} and {clstr_file2}")
