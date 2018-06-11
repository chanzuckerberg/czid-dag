import multiprocessing
import os
import threading

from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.s3 as s3
import idseq_dag.util.log as log


class PipelineStepRunStar(PipelineStep):
    MAX_INPUT_READS = 75 * 1000 * 1000

    def run(self):
        input_file_0 = self.input_files_local[0][0]
        input_file_1 = self.input_files_local[0][1]
        output_file_0 = self.output_files_local()[0]
        output_file_1 = self.output_files_local()[1]

        STAR_COUNTS_OUT = ""
        STAR_GENOME = ""
        REF_DIR = ""
        SCRATCH_DIR = ""
        total_counts_from_star = {}
        SAMPLE_S3_OUTPUT_PATH = ""

        def result_path(basename):
            return os.path.join(self.output_dir_local, basename)

        """Run STAR to filter out host reads."""
        # what to do about STAR_COUNTS_OUT?
        star_outputs = [STAR_COUNTS_OUT, output_file_0, output_file_1]
        fastq_files = [input_file_0, input_file_1]
        num_fastqs = len(fastq_files)
        gene_count_output = None

        def unmapped_files_in(some_dir):
            return [
                "%s/Unmapped.out.mate%d" % (some_dir, i + 1)
                for i in range(num_fastqs)
            ]

        # TODO: what to do about STAR_GENOME?
        genome_dir = s3.fetch_genome(STAR_GENOME, REF_DIR)
        assert genome_dir is not None

        # Check if parts.txt file exists. If so, use the new version of partitioned
        # indices. Otherwise, stay put.
        parts_file = os.path.join(genome_dir, "parts.txt")
        assert os.path.isfile(parts_file)
        with open(parts_file, 'rb') as parts_f:
            num_parts = int(parts_f.read())

        unmapped = fastq_files
        for part_idx in range(num_parts):
            tmp_result_dir = "%s/star-part-%d" % (SCRATCH_DIR, part_idx)
            genome_part = "%s/part-%d" % (genome_dir, part_idx)
            self.run_star_part(tmp_result_dir, genome_part, unmapped, part_idx == 0)

            unmapped = self.sync_pairs(unmapped_files_in(tmp_result_dir))

            # Run part 0 in gene-counting mode:
            # (a) ERCCs are doped into part 0 and we want their counts
            # (b) if there is only 1 part (e.g. human), the host gene counts also
            # make sense
            # (c) at part 0, we can also extract out total input reads and if the
            # total_counts is exactly the same as MAX_INPUT_READS then we know the
            # input file is truncated.
            if part_idx == 0:
                gene_count_file = os.path.join(tmp_result_dir,
                                               "ReadsPerGene.out.tab")
                self.extract_total_counts_from_star_output(tmp_result_dir, num_fastqs,
                                                           total_counts_from_star)
                if os.path.isfile(gene_count_file):
                    gene_count_output = gene_count_file

        result_files = [gene_count_output] + unmapped
        for i, f in enumerate(result_files):
            if f is not None:
                output_i = result_path(star_outputs[i])
                command.execute("mv %s %s;" % (f, output_i))
                # TODO: question: which one should be called for this
                s3.upload_with_retries(output_i, SAMPLE_S3_OUTPUT_PATH + "/")

        # Cleanup
        command.execute("cd %s; rm -rf *" % SCRATCH_DIR)
        log.write("Finished running STAR.")

    def run_star_part(self, output_dir, genome_dir, fastq_files, count_genes=False):
        command.execute("mkdir -p %s" % output_dir)
        star_command_params = [
            'cd', output_dir, ';', 'STAR', '--outFilterMultimapNmax', '99999',
            '--outFilterScoreMinOverLread', '0.5', '--outFilterMatchNminOverLread',
            '0.5', '--outReadsUnmapped', 'Fastx', '--outFilterMismatchNmax', '999',
            '--outSAMmode', 'None', '--clip3pNbases', '0', '--runThreadN',
            str(multiprocessing.cpu_count()), '--genomeDir', genome_dir,
            '--readFilesIn', " ".join(fastq_files)
        ]
        if fastq_files[0][-3:] == '.gz':
            # Create a custom decompressor which does "zcat $input_file | head -
            # ${max_lines}"
            cmd = "echo 'zcat ${2} | head -${1}' > %s/gzhead; " % genome_dir
            command.execute(cmd)
            max_lines = self.max_input_lines(fastq_files[0])
            star_command_params += [
                '--readFilesCommand',
                '"sh %s/gzhead %d"' % (genome_dir, max_lines)
            ]
        path = "%s/sjdbList.fromGTF.out.tab" % genome_dir
        if count_genes and os.path.isfile(path):
            star_command_params += ['--quantMode', 'GeneCounts']
        cmd = " ".join(star_command_params), os.path.join(output_dir,
                                                          "Log.progress.out")
        command.execute(cmd)

    def write_lines(self, of, lines):
        for l in lines:
            of.write(l)

    def handle_outstanding_read(self, r0, r0id, outstanding_r0, outstanding_r1, of0, of1,
                                mem, max_mem):
        # If read r0 completes an outstanding r1, output the pair (r0, r1).
        # Else r0 becomes outstanding, so in future some r1 may complete it.
        if r0id:
            if r0id in outstanding_r1:
                self.write_lines(of0, r0)
                self.write_lines(of1, outstanding_r1.pop(r0id))
                mem -= 1
            else:
                outstanding_r0[r0id] = r0
                mem += 1
                if mem > max_mem:
                    max_mem = mem
        return mem, max_mem

    def sync_pairs_work(self, of0, of1, if0, if1):
        # TODO: Use this as a template for merging fasta?
        outstanding_r0 = {}
        outstanding_r1 = {}
        mem = 0
        max_mem = 0
        while True:
            r0, r0id = self.get_read(if0)
            r1, r1id = self.get_read(if1)
            if not r0 and not r1:
                break
            if r0id == r1id:
                # If the input pairs are already synchronized, we take this branch
                # on every iteration.
                self.write_lines(of0, r0)
                self.write_lines(of1, r1)
            else:
                mem, max_mem = self.handle_outstanding_read(r0, r0id, outstanding_r0,
                                                            outstanding_r1, of0, of1,
                                                            mem, max_mem)
                mem, max_mem = self.handle_outstanding_read(r1, r1id, outstanding_r1,
                                                            outstanding_r0, of1, of0,
                                                            mem, max_mem)
        return outstanding_r0, outstanding_r1, max_mem

    def sync_pairs(self, fastq_files, max_discrepancies=0):
        """The given fastq_files contain the same read IDs but in different order.
        Output the same data in synchronized order. Omit up to max_discrepancies
        if necessary. If more must be suppressed, raise assertion.
        """
        if len(fastq_files) != 2:
            return fastq_files

        output_fnames = [ifn + ".synchronized_pairs.fq" for ifn in fastq_files]
        with open(fastq_files[0], "rb") as if_0:
            with open(fastq_files[1], "rb") as if_1:
                with open(output_fnames[0], "wb") as of_0:
                    with open(output_fnames[1], "wb") as of_1:
                        outstanding_r0, outstanding_r1, max_mem = self.sync_pairs_work(
                            of_0, of_1, if_0, if_1)
        if max_mem:
            # This will be printed if some pairs were out of order.
            warning_message = "WARNING: Pair order out of sync in {fqf}. " \
                              "Synchronized using RAM for {max_mem} pairs.".format(
                fqf=fastq_files, max_mem=max_mem)
            log.write(warning_message)

        discrepancies_count = len(outstanding_r0) + len(outstanding_r1)
        if discrepancies_count:
            warning_message = "WARNING: Found {dc} broken pairs in {fqf}, e.g., " \
                              "{example}.".format(
                dc=discrepancies_count,
                fqf=fastq_files,
                example=(outstanding_r0 or outstanding_r1).popitem()[0])
            log.write(warning_message)
            assert discrepancies_count <= max_discrepancies, warning_message
        return output_fnames

    def extract_total_counts_from_star_output(self, result_dir, num_fastqs,
                                              total_counts_from_star):
        """Grab the total reads from the Log.final.out file."""
        log_file = os.path.join(result_dir, "Log.final.out")
        cmd = "grep 'Number of input reads' %s" % log_file
        total_reads = command.execute_with_output(cmd).split("\t")[1]
        total_reads = int(total_reads)
        # If it's exactly the same, it must have been truncated.
        if total_reads == self.MAX_INPUT_READS:
            total_counts_from_star['truncated'] = 1
        total_counts_from_star['total_reads'] = total_reads * num_fastqs

    def max_input_lines(self, input_file):
        """Returning number of lines corresponding to MAX_INPUT_READS based on file
        type.
        """
        if "fasta" in input_file:
            return self.MAX_INPUT_READS * 2
        # Assume it's FASTQ
        return self.MAX_INPUT_READS * 4

    def get_read(self, f):
        # The FASTQ format specifies that each read consists of 4 lines,
        # the first of which begins with @ followed by read ID.
        read, rid = [], None
        line = f.readline()
        if line:
            assert line[0] == "@"
            rid = line.split("\t", 1)[0].strip()
            read.append(line)
            read.append(f.readline())
            read.append(f.readline())
            read.append(f.readline())
        return read, rid
