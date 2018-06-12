import multiprocessing
import os

from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.s3 as s3
import idseq_dag.util.log as log


class PipelineStepRunStar(PipelineStep):
    MAX_INPUT_READS = 75 * 1000 * 1000

    def run(self):
        """Run STAR to filter out host reads."""
        print("Input files: " + str(self.input_files))
        print("Output files: " + str(self.output_files))
        print("Output dir local: " + self.output_dir_local)
        print("Output dir s3: " + self.output_dir_s3)
        print("Ref dir local: " + self.ref_dir_local)
        print("Additional files: " + str(self.additional_files))
        print("Additional attributes: " + str(self.additional_attributes))
        print("Output files local: " + str(self.output_files_local()))
        print("Input files local: " + str(self.input_files_local))

        # Setup
        input_files = self.input_files_local[0][0:2]
        star_genome = self.additional_files["star_genome"]
        ref_dir = self.ref_dir_local
        scratch_dir = os.path.join(self.output_dir_local, "scratch")
        total_counts_from_star = {}

        num_fastqs = len(self.input_files[0])

        def unmapped_files_in(some_dir):
            return [
                f"{some_dir}/Unmapped.out.mate{i+1}" for i in range(num_fastqs)
            ]

        genome_dir = s3.fetch_from_s3(star_genome, ref_dir, True, True, True)
        assert genome_dir is not None

        # Check if parts.txt file exists. If so, use the new version of
        # partitioned indices. Otherwise, stay put.
        parts_file = os.path.join(genome_dir, "parts.txt")
        assert os.path.isfile(parts_file)
        with open(parts_file, 'rb') as parts_f:
            num_parts = int(parts_f.read())

        unmapped = input_files
        for part_idx in range(num_parts):
            tmp_result_dir = f"{scratch_dir}/star-part-{part_idx}"
            genome_part = f"{genome_dir}/part-{part_idx}"
            count_genes = part_idx == 0
            self.run_star_part(tmp_result_dir, genome_part, unmapped,
                               count_genes)

            unmapped = self.sync_pairs(unmapped_files_in(tmp_result_dir))

            # Run part 0 in gene-counting mode:
            # (a) ERCCs are doped into part 0 and we want their counts.
            # (b) If there is only 1 part (e.g. human), the host gene counts also
            # make sense.
            # (c) At part 0, we can also extract out total input reads and if the
            # total_counts is exactly the same as MAX_INPUT_READS then we know the
            # input file is truncated.
            if part_idx == 0:
                gene_count_file = os.path.join(tmp_result_dir,
                                               "ReadsPerGene.out.tab")
                self.extract_total_counts_from_star_output(
                    tmp_result_dir, num_fastqs, total_counts_from_star)
                if os.path.isfile(gene_count_file):
                    gene_count_output = gene_count_file
                    self.additional_files_to_upload.append(gene_count_output)

        # Cleanup
        command.execute("cd %s; rm -rf *" % scratch_dir)
        log.write("Finished running STAR.")

    def run_star_part(self,
                      output_dir,
                      genome_dir,
                      fastq_files,
                      count_genes=False):
        command.execute("mkdir -p %s" % output_dir)
        cpus = str(multiprocessing.cpu_count())
        star_command_params = [
            'cd', output_dir, ';', 'STAR', '--outFilterMultimapNmax', '99999',
            '--outFilterScoreMinOverLread', '0.5',
            '--outFilterMatchNminOverLread', '0.5', '--outReadsUnmapped',
            'Fastx', '--outFilterMismatchNmax', '999', '--outSAMmode', 'None',
            '--clip3pNbases', '0', '--runThreadN',
            cpus, '--genomeDir', genome_dir,
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
        cmd = " ".join(star_command_params), os.path.join(
            output_dir, "Log.progress.out")
        command.execute(cmd)

    def handle_outstanding_read(self, r0, r0id, outstanding_r0, outstanding_r1,
                                of0, of1, mem, max_mem):
        # If read r0 completes an outstanding r1, output the pair (r0, r1).
        # Else r0 becomes outstanding, so in future some r1 may complete it.
        if r0id:
            if r0id in outstanding_r1:
                write_lines(of0, r0)
                write_lines(of1, outstanding_r1.pop(r0id))
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
                # If the input pairs are already synchronized, we take this
                # branch on every iteration.
                write_lines(of0, r0)
                write_lines(of1, r1)
            else:
                mem, max_mem = self.handle_outstanding_read(
                    r0, r0id, outstanding_r0, outstanding_r1, of0, of1, mem,
                    max_mem)
                mem, max_mem = self.handle_outstanding_read(
                    r1, r1id, outstanding_r1, outstanding_r0, of1, of0, mem,
                    max_mem)
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
            msg = "WARNING: Pair order out of sync in {fqf}. " \
                              "Synchronized using RAM for {max_mem} pairs.".format(
                fqf=fastq_files, max_mem=max_mem)
            log.write(msg)

        discrepancies_count = len(outstanding_r0) + len(outstanding_r1)
        if discrepancies_count:
            msg = "WARNING: Found {dc} broken pairs in {fqf}, e.g., " \
                              "{example}.".format(
                dc=discrepancies_count,
                fqf=fastq_files,
                example=(outstanding_r0 or outstanding_r1).popitem()[0])
            log.write(msg)
            assert discrepancies_count <= max_discrepancies, msg
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
        """Return number of lines corresponding to MAX_INPUT_READS based on file
        type.
        """
        res = self.MAX_INPUT_READS * 2
        if "fasta" not in input_file:  # Assume it's FASTQ
            res *= 2
        return res

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


def write_lines(of, lines):
    for l in lines:
        of.write(l)
