import multiprocessing
import os
import random
import shlex
import shutil
import threading
import time
import traceback

import boto3

from idseq_dag.engine.pipeline_step import InputFileErrors, PipelineStep

import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.count as count
import idseq_dag.util.log as log
import idseq_dag.util.m8 as m8

from idseq_dag.util.s3 import fetch_from_s3, fetch_reference
from idseq_dag.util.trace_lock import TraceLock

from idseq_dag.util.lineage import DEFAULT_BLACKLIST_S3, DEFAULT_WHITELIST_S3
from idseq_dag.util.m8 import NT_MIN_ALIGNMENT_LEN

MAX_CHUNKS_IN_FLIGHT = 32
CHUNK_MAX_ATTEMPTS = 3
CHUNK_ATTEMPT_TIMEOUT = 60 * 60 * 12  # 12 hours
GSNAP_CHUNK_SIZE = 60000
RAPSEARCH_CHUNK_SIZE = 80000

class PipelineStepRunAlignmentRemotely(PipelineStep):
    """ Runs gsnap/rapsearch2 remotely.

    For GSNAP:
    ```
    gsnapl
    -A m8
    --batch=0
    --use-shared-memory=0
    --gmap-mode=none
    --npaths=100
    --ordered
    -t 48
    --max-mismatches=40
    -D {remote_index_dir}
    -d nt_k16
    {remote_input_files} > {multihit_remote_outfile}
    ```

    GSNAP documentation is available [here](http://research-pub.gene.com/gmap/).
    -t (threads): For example, r5d.24xlarge machines have 96 vCPUs. These are the machines used by the alignment batch compute environment.
    Each batch job is allotted 48 vcpus. Use 48 threads and each instance will be able to concurrently process 2 chunks

    For Rapsearch2:
    ```
    rapsearch
    -d {remote_index_dir}/nr_rapsearch
    -e -6
    -l 10
    -a T
    -b 0
    -v 50
    -z 24
    -q {remote_input_files}
    -o {multihit_remote_outfile}
    ```

    Rapsearch2 documentation is available [here](http://omics.informatics.indiana.edu/mg/RAPSearch2/).
    """

    @staticmethod
    def _alignment_algorithm_inputs(host_filter_outputs):
        # Destructuring-bind for host filter outputs.
        # The host filter outputs either 1 file for unpaired reads, let's call it R1;
        # or, 3 files for paired-end reads, which can be succinctly described
        # as [R1, R2, R1/1 + R2/2].  All are usually named gsnap_filter_*.
        try:
            # Unpaired?
            gsnap_filter_1, = host_filter_outputs
            return {
                "gsnap": [gsnap_filter_1],
                "rapsearch2": [gsnap_filter_1]
            }
        except:
            # Paired!
            gsnap_filter_1, gsnap_filter_2, gsnap_filter_merged = host_filter_outputs
            return {
                "gsnap": [gsnap_filter_1, gsnap_filter_2],
                "rapsearch2": [gsnap_filter_merged]
            }

    def validate_input_files(self):
        # first two files are gsnap_filter_1.fa and gsnap_filter_2.fa
        if not count.files_have_min_reads(self.input_files_local[0][:-1], 1):
            self.input_file_error = InputFileErrors.INSUFFICIENT_READS

    def __init__(self, *args, **kwrds):
        PipelineStep.__init__(self, *args, **kwrds)
        # TODO: (tmorse) remove service compatibility https://jira.czi.team/browse/IDSEQ-2568
        self.alignment_algorithm = self.additional_attributes.get("alignment_algorithm", self.additional_attributes.get("service"))
        assert self.alignment_algorithm in ("gsnap", "rapsearch2")
        self.chunks_in_flight = threading.Semaphore(MAX_CHUNKS_IN_FLIGHT)
        self.chunks_result_dir_local = os.path.join(self.output_dir_local, "chunks")
        self.chunks_result_dir_s3 = os.path.join(self.output_dir_s3, "chunks")
        command.make_dirs(self.chunks_result_dir_local)

    def count_reads(self):
        pass

    def run(self):
        ''' Run alignmment remotely '''

        alignment_algorithm_inputs = PipelineStepRunAlignmentRemotely._alignment_algorithm_inputs(self.input_files_local[0])
        cdhit_cluster_sizes_path, = self.input_files_local[1]
        output_m8, deduped_output_m8, output_hitsummary, output_counts_with_dcr_json = self.output_files_local()
        assert output_counts_with_dcr_json.endswith("_with_dcr.json"), self.output_files_local()

        self.run_remotely(alignment_algorithm_inputs[self.alignment_algorithm], output_m8)

        # get database
        lineage_db = fetch_reference(self.additional_files["lineage_db"], self.ref_dir_local)
        accession2taxid_db = fetch_reference(self.additional_files["accession2taxid_db"], self.ref_dir_local, allow_s3mi=True)

        min_alignment_length = NT_MIN_ALIGNMENT_LEN if self.alignment_algorithm == 'gsnap' else 0
        m8.call_hits_m8(output_m8, lineage_db, accession2taxid_db,
                        deduped_output_m8, output_hitsummary, min_alignment_length)

        db_type = 'NT' if self.alignment_algorithm == 'gsnap' else 'NR'
        evalue_type = 'log10' if self.alignment_algorithm == 'rapsearch2' else 'raw'

        deuterostome_db = None
        if self.additional_files.get("deuterostome_db"):
            deuterostome_db = fetch_reference(self.additional_files["deuterostome_db"],
                                              self.ref_dir_local, allow_s3mi=True)

        blacklist_s3_file = self.additional_attributes.get('taxon_blacklist', DEFAULT_BLACKLIST_S3)
        taxon_blacklist = fetch_reference(blacklist_s3_file, self.ref_dir_local)

        taxon_whitelist = None
        if self.additional_attributes.get("use_taxon_whitelist"):
            taxon_whitelist = fetch_reference(self.additional_files.get("taxon_whitelist", DEFAULT_WHITELIST_S3),
                                              self.ref_dir_local)

        m8.generate_taxon_count_json_from_m8(
            deduped_output_m8, output_hitsummary, evalue_type, db_type,
            lineage_db, deuterostome_db, taxon_whitelist, taxon_blacklist, cdhit_cluster_sizes_path,
            output_counts_with_dcr_json)

    def run_remotely(self, input_fas, output_m8):
        # Split files into chunks for performance
        chunk_size = GSNAP_CHUNK_SIZE if self.alignment_algorithm == "gsnap" else RAPSEARCH_CHUNK_SIZE
        part_suffix, input_chunks = self.chunk_input(input_fas, chunk_size)
        chunk_count = len(input_chunks)

        # Process chunks
        chunk_output_files = [None] * chunk_count
        chunk_threads = []
        mutex = TraceLock("run_remotely", threading.RLock())
        # Randomize execution order for performance
        randomized = list(enumerate(input_chunks))
        random.shuffle(randomized)

        try:
            for n, chunk_input_files in randomized:
                self.chunks_in_flight.acquire()
                self.check_for_errors(mutex, chunk_output_files, input_chunks, self.alignment_algorithm)
                t = threading.Thread(
                    target=PipelineStepRunAlignmentRemotely.run_chunk_wrapper,
                    args=[
                        self.chunks_in_flight, chunk_output_files, n, mutex, self.run_chunk,
                        [
                            part_suffix, chunk_input_files, chunk_count, True
                        ]
                    ])
                t.start()
                chunk_threads.append(t)

        finally:
            # Check chunk completion
            for ct in chunk_threads:
                ct.join()

        self.check_for_errors(mutex, chunk_output_files, input_chunks, self.alignment_algorithm)

        assert None not in chunk_output_files
        # Concatenate the pieces and upload results
        self.concatenate_files(chunk_output_files, output_m8)

    def chunk_input(self, input_files, chunksize):
        """Chunk input files into pieces for performance and parallelism."""
        part_lists = []  # Lists of partial files
        known_nlines = None
        part_suffix = ""
        chunk_nlines = chunksize * 2

        for input_file in input_files:
            # Count number of lines in the file
            cmd_output = command.execute_with_output(
                command_patterns.SingleCommand(
                    cmd="wc",
                    args=[
                        "-l",
                        input_file
                    ]
                )
            )
            nlines = int(cmd_output.strip().split()[0])
            # Number of lines should be the same in paired files
            if known_nlines is not None:
                msg = "Mismatched line counts in supposedly paired files: {}".format(
                    input_files)
                assert nlines == known_nlines, msg
            known_nlines = nlines

            # Set number of pieces and names
            numparts = (nlines + chunk_nlines - 1) // chunk_nlines
            ndigits = len(str(numparts - 1))
            part_suffix = "-chunksize-%d-numparts-%d-part-" % (chunksize, numparts)
            out_prefix_base = os.path.basename(input_file) + part_suffix
            out_prefix = os.path.join(self.chunks_result_dir_local, out_prefix_base)

            # Split large file into smaller named pieces
            command.execute(
                command_patterns.SingleCommand(
                    cmd="split",
                    args=[
                        "-a",
                        ndigits,
                        "--numeric-suffixes",
                        "-l",
                        chunk_nlines,
                        input_file,
                        out_prefix
                    ]
                )
            )
            command.execute_with_retries(
                command_patterns.SingleCommand(
                    cmd="aws",
                    args=[
                        "s3",
                        "sync",
                        "--only-show-errors",
                        os.path.join(self.chunks_result_dir_local, ""),
                        os.path.join(self.chunks_result_dir_s3, ""),
                        "--exclude",
                        "*",
                        "--include",
                        out_prefix_base + "*"
                    ]
                )
            )

            # Get the partial file names
            partial_files = []
            paths = command.glob(
                glob_pattern=out_prefix + "*",
                strip_folder_names=True
            )
            partial_files.extend(paths)

            # Check that the partial files match our expected chunking pattern
            pattern = "{:0%dd}" % ndigits
            expected_partial_files = [(out_prefix_base + pattern.format(i))
                                      for i in range(numparts)]
            msg = "something went wrong with chunking: {} != {}".format(
                partial_files, expected_partial_files)
            assert expected_partial_files == partial_files, msg
            part_lists.append(partial_files)

        # Ex: [["input_R1.fasta-part-1", "input_R2.fasta-part-1"],
        # ["input_R1.fasta-part-2", "input_R2.fasta-part-2"],
        # ["input_R1.fasta-part-3", "input_R2.fasta-part-3"],...]
        input_chunks = [list(part) for part in zip(*part_lists)]
        return part_suffix, input_chunks

    @staticmethod
    def check_for_errors(mutex, chunk_output_files, input_chunks, alignment_algorithm):
        with mutex:
            if "error" in chunk_output_files:
                # We already have per-chunk retries to address transient (non-deterministic) system issues.
                # If a chunk fails on all retries, it must be due to a deterministic & reproducible problem
                # with the chunk data or command options, so we should not even attempt the other chunks.
                err_i = chunk_output_files.index("error")
                raise RuntimeError("All retries failed for {} chunk {}.".format(
                    alignment_algorithm, input_chunks[err_i]))

    @staticmethod
    def concatenate_files(chunk_output_files, output_m8):
        with log.log_context("run_alignment_remote.concatenate_files", {"chunk_output_files": chunk_output_files}):
            with open(output_m8, 'wb') as outf:
                for f in chunk_output_files:
                    with log.log_context("run_alignment_remote.concatenate_files#chunk", {"f": f}):
                        with open(f, 'rb') as fd:
                            shutil.copyfileobj(fd, outf)

    @staticmethod
    def run_chunk_wrapper(chunks_in_flight, chunk_output_files, n, mutex, target, args):
        result = "error"
        try:
            result = target(*args)
        except:
            with mutex:
                log.write(traceback.format_exc())
        finally:
            with mutex:
                chunk_output_files[n] = result
            chunks_in_flight.release()

    @staticmethod
    def _get_job_status(client, job_id):
        response = client.describe_jobs(jobs=[job_id])
        return response['jobs'][0]['status']

    def run_chunk(self, part_suffix, input_files, chunk_count, lazy_run):
        """
        Dispatch a chunk to worker machines for distributed GSNAP or RAPSearch
        group machines and handle their execution.
        """

        chunk_id = int(input_files[0].split(part_suffix)[-1])
        multihit_basename = f"multihit-{self.alignment_algorithm}-out{part_suffix}{chunk_id}.m8"
        multihit_local_outfile = os.path.join(self.chunks_result_dir_local, multihit_basename)
        multihit_s3_outfile = os.path.join(self.chunks_result_dir_s3, multihit_basename)

        if lazy_run and fetch_from_s3(multihit_s3_outfile, multihit_local_outfile, okay_if_missing=True, allow_s3mi=False):
            log.write(f"finished alignment for chunk {chunk_id} with {self.alignment_algorithm} by lazily fetching last result")
            return multihit_local_outfile

        # TODO: (tmorse) remove compat hack https://jira.czi.team/browse/IDSEQ-2568
        deployment_environment = os.environ.get("DEPLOYMENT_ENVIRONMENT", self.additional_attributes.get("environment"))
        priority_name = os.environ.get("PRIORITY_NAME", "normal")
        provisioning_model = os.environ.get("PROVISIONING_MODEL", "EC2")

        index_dir_suffix = self.additional_attributes["index_dir_suffix"]

        job_name = f"idseq-{deployment_environment}-{self.alignment_algorithm}"
        job_queue = f"idseq-{deployment_environment}-{self.alignment_algorithm}-{provisioning_model}-{index_dir_suffix}-{priority_name}"
        job_definition = f"idseq-{deployment_environment}-{self.alignment_algorithm}"

        environment = [{
            'name': f"INPUT_PATH_{i}",
            'value': os.path.join(self.chunks_result_dir_s3, input_file),
        } for i, input_file in enumerate(input_files)] + [{
            'name': "OUTPUT_PATH",
            'value': multihit_s3_outfile,
        }]

        session = boto3.session.Session()
        client = session.client("batch")
        response = client.submit_job(
            jobName=job_name,
            jobQueue=job_queue,
            jobDefinition=job_definition,
            containerOverrides={"environment": environment},
            retryStrategy={"attempts": CHUNK_MAX_ATTEMPTS},
            timeout={"attemptDurationSeconds": CHUNK_ATTEMPT_TIMEOUT}
        )
        job_id = response["jobId"]

        log.write(f"waiting for {self.alignment_algorithm} batch queue for chunk {chunk_id}.")
        total_timeout = CHUNK_MAX_ATTEMPTS * CHUNK_ATTEMPT_TIMEOUT
        start = time.time()
        end = start + total_timeout
        delay = 5 * chunk_count # ~1 chunk every 5 seconds to avoid throttling
        while time.time() < end:
            status = PipelineStepRunAlignmentRemotely._get_job_status(client, job_id)
            if status == "SUCCEEDED":
                break
            if status == "FAILED":
                log.log_event("alignment_batch_error", values={"job_id": job_id})
                raise Exception("chunk alignment failed")
            time.sleep(delay)

        fetch_from_s3(multihit_s3_outfile, multihit_local_outfile, okay_if_missing=True, allow_s3mi=False)

        log.write(f"finished alignment for chunk {chunk_id} on {self.alignment_algorithm}")

        return multihit_local_outfile

    def step_description(self, require_docstrings=False):
        if (self.alignment_algorithm == "gsnap"):
            return """
                Runs gsnap remotely.

                ```
                gsnapl
                -A m8
                --batch=0
                --use-shared-memory=0
                --gmap-mode=none
                --npaths=100
                --ordered
                -t 48
                --max-mismatches=40
                -D {remote_index_dir}
                -d nt_k16
                {remote_input_files} > {multihit_remote_outfile}
                ```

                GSNAP documentation is available [here](http://research-pub.gene.com/gmap/).
            """
        elif self.alignment_algorithm == "rapsearch2":
            return """
                Runs rapsearch2 remotely.

                ```
                rapsearch
                -d {remote_index_dir}/nr_rapsearch
                -e -6
                -l 10
                -a T
                -b 0
                -v 50
                -z 24
                -q {remote_input_files}
                -o {multihit_remote_outfile}
                ```

                Rapsearch2 documentation is available [here](http://omics.informatics.indiana.edu/mg/RAPSearch2/).
            """
        # If neither, then return default step_description method.
        return super(PipelineStepRunAlignmentRemotely, self).step_description()
