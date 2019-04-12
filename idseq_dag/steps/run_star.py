import multiprocessing
import os
import json

from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.s3 as s3
import idseq_dag.util.count as count
import idseq_dag.util.validate_constants as vc

class PipelineStepRunStar(PipelineStep):
    """ Implements the step for running STAR.
    The attributes 'validated_input_counts_file' and 'sequence_input_files' should be
    initialized in the run() method of any children that are not for execution directly
    on the output of the 'run_validate_input' step.
    Note that these attributes cannot be set in the __init__ method because required
    file path information only becomes available once the wait_for_input_files() has run.
    """
    def __init__(self, *args):
        super().__init__(*args)
        self.sequence_input_files = None
        self.validated_input_counts_file = None

    def run(self):
        """Run STAR to filter out host reads."""
        # Setup
        if self.sequence_input_files is not None and self.validated_input_counts_file is not None:
            # if input and counts files are available, use them
            validated_input_counts_file = self.validated_input_counts_file
            input_files = self.sequence_input_files
        else:
            if self.additional_attributes.get('raw_fastq_input'):
                # make this step not depending on validation
                validated_input_counts_file = None
                input_files = self.input_files_local[0][0:2]
            else:
                validated_input_counts_file = self.input_files_local[0][0]
                input_files = self.input_files_local[0][1:3]

        # Decide if we want to run STARLong
        use_starlong = False
        if validated_input_counts_file:
            with open(validated_input_counts_file) as validated_input_counts_f:
                validated_input_counts = json.load(validated_input_counts_f)

            use_starlong = validated_input_counts[vc.BUCKET_LONG] > 1 or \
                           validated_input_counts[vc.BUCKET_TOO_LONG] > 1


        num_inputs = len(input_files)
        scratch_dir = os.path.join(self.output_dir_local, "scratch_star")

        output_files_local = self.output_files_local()
        output_gene_file = self.additional_attributes.get("output_gene_file")

        genome_dir = s3.fetch_from_s3(
            self.additional_files["star_genome"],
            self.ref_dir_local,
            allow_s3mi=True,
            auto_untar=True)

        # Check parts file for the number of partitioned indexes
        parts_file = os.path.join(genome_dir, "parts.txt")
        assert os.path.isfile(parts_file)
        with open(parts_file, 'rb') as parts_f:
            num_parts = int(parts_f.read())

        # Run STAR on each partition and save the unmapped read info
        unmapped = input_files


        for part_idx in range(num_parts):
            tmp = f"{scratch_dir}/star-part-{part_idx}"
            genome_part = f"{genome_dir}/part-{part_idx}"
            count_genes = part_idx == 0
            self.run_star_part(tmp, genome_part, unmapped, count_genes, use_starlong)

            unmapped = PipelineStepRunStar.sync_pairs(
                PipelineStepRunStar.unmapped_files_in(tmp, num_inputs))

            # Run part 0 in gene-counting mode:
            # (a) ERCCs are doped into part 0 and we want their counts.
            # (b) If there is only 1 part (e.g. human), the host gene counts also
            # make sense.
            if part_idx == 0:
                gene_count_file = os.path.join(tmp, "ReadsPerGene.out.tab")
                if os.path.isfile(gene_count_file) and output_gene_file:
                    moved = os.path.join(self.output_dir_local,
                                         output_gene_file)
                    command.execute(f"mv {gene_count_file} {moved}")
                    self.additional_files_to_upload.append(moved)

        # Cleanup
        for src, dst in zip(unmapped, output_files_local):
            command.execute(f"mv {src} {dst}")  # Move out of scratch dir
        command.execute("cd %s; rm -rf *" % scratch_dir)

    def count_reads(self):
        self.should_count_reads = True
        self.counts_dict[self.name] = count.reads_in_group(self.output_files_local()[0:2])

    def run_star_part(self,
                      output_dir,
                      genome_dir,
                      input_files,
                      count_genes,
                      use_starlong):
        command.execute("mkdir -p %s" % output_dir)
        cpus = str(multiprocessing.cpu_count())
        params = [
            'cd', output_dir, ';', 'STARlong' if use_starlong else 'STAR',
            '--outFilterMultimapNmax', '99999',
            '--outFilterScoreMinOverLread', '0.5',
            '--outFilterMatchNminOverLread', '0.5',
            '--outReadsUnmapped', 'Fastx',
            '--outFilterMismatchNmax', '999',
            '--outSAMmode', 'None',
            '--clip3pNbases', '0',
            '--runThreadN', cpus,
            '--genomeDir', genome_dir,
            '--readFilesIn', " ".join(input_files)
        ]
        if use_starlong:
            params+= [
                '--seedSearchStartLmax', '20',
                '--seedPerReadNmax', '100000',
                '--seedPerWindowNmax', '1000',
                '--alignTranscriptsPerReadNmax', '100000',
                '--alignTranscriptsPerWindowNmax', '10000']
        count_file = f"{genome_dir}/sjdbList.fromGTF.out.tab"

        if count_genes and os.path.isfile(count_file):
            params += ['--quantMode', 'GeneCounts']
        cmd = " ".join(params)
        command.execute(cmd)

    @staticmethod
    def handle_outstanding_read(r0, r0id, outstanding_r0, outstanding_r1, of0,
                                of1, mem, max_mem):
        # If read r0 completes an outstanding r1, output the pair (r0, r1).
        # Else r0 becomes outstanding, so in future some r1 may complete it.
        if r0id:
            if r0id in outstanding_r1:
                PipelineStepRunStar.write_lines(of0, r0)
                PipelineStepRunStar.write_lines(of1, outstanding_r1.pop(r0id))
                mem -= 1
            else:
                outstanding_r0[r0id] = r0
                mem += 1
                if mem > max_mem:
                    max_mem = mem
        return mem, max_mem

    @staticmethod
    def sync_pairs_work(of0, of1, if0, if1):
        # TODO: Use this as a template for merging fasta?
        outstanding_r0 = {}
        outstanding_r1 = {}
        mem = 0
        max_mem = 0
        while True:
            r0, r0id = PipelineStepRunStar.get_read(if0)
            r1, r1id = PipelineStepRunStar.get_read(if1)
            if not r0 and not r1:
                break
            if r0id == r1id:
                # If the input pairs are already synchronized, we take this
                # branch on every iteration.
                PipelineStepRunStar.write_lines(of0, r0)
                PipelineStepRunStar.write_lines(of1, r1)
            else:
                mem, max_mem = PipelineStepRunStar.handle_outstanding_read(
                    r0, r0id, outstanding_r0, outstanding_r1, of0, of1, mem,
                    max_mem)
                mem, max_mem = PipelineStepRunStar.handle_outstanding_read(
                    r1, r1id, outstanding_r1, outstanding_r0, of1, of0, mem,
                    max_mem)
        return outstanding_r0, outstanding_r1, max_mem

    @staticmethod
    def sync_pairs(fastq_files, max_discrepancies=0):
        """The given fastq_files contain the same read IDs but in different order.
        Output the same data in synchronized order. Omit up to max_discrepancies
        if necessary. If more must be suppressed, raise assertion.
        """
        if len(fastq_files) != 2:
            return fastq_files

        output_fnames = [ifn + ".synchronized_pairs.fq" for ifn in fastq_files]
        with open(fastq_files[0], "rb") as if_0, open(fastq_files[1],
                                                      "rb") as if_1:
            with open(output_fnames[0], "wb") as of_0, open(
                    output_fnames[1], "wb") as of_1:
                outstanding_r0, outstanding_r1, max_mem = PipelineStepRunStar.sync_pairs_work(
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

    @staticmethod
    def get_read(f):
        # The FASTQ/FASTA format specifies that each read consists of 4/2 lines,
        # the first of which begins with @/> followed by read ID.
        read, rid = [], None
        line = f.readline()
        if line:
            if line[0] == 64:  # Equivalent to '@', fastq format
                rid = line.decode('utf-8').split('\t', 1)[0].strip()
                read.append(line)
                for i in range(3):
                    read.append(f.readline())
            elif line[0] == 62: # Equivalent to '>', fasta format
                rid = line.decode('utf-8').split('\t', 1)[0].strip()
                read.append(line)
                read.append(f.readline())
            else:
                raise RuntimeError("sync pair failed. unknown ")
        return read, rid

    @staticmethod
    def write_lines(of, lines):
        for l in lines:
            of.write(l)

    @staticmethod
    def unmapped_files_in(folder, num_inputs):
        return [f"{folder}/Unmapped.out.mate{i+1}" for i in range(num_inputs)]
