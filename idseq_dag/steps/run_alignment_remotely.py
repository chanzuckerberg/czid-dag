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
from idseq_dag.util.server import ASGInstance, ChunkStatus, chunk_status_tracker
from idseq_dag.util.trace_lock import TraceLock

from idseq_dag.util.lineage import DEFAULT_BLACKLIST_S3, DEFAULT_WHITELIST_S3
from idseq_dag.util.m8 import NT_MIN_ALIGNMENT_LEN

MAX_CONCURRENT_CHUNK_UPLOADS = 4
MAX_CHUNKS_IN_FLIGHT = 32
CORRECT_NUMBER_OF_OUTPUT_COLUMNS = 12
CHUNK_MAX_ATTEMPTS = 3
CHUNK_ATTEMPT_TIMEOUT = 60 * 60 * 12  # 12 hours
CHUNK_COMPLETE_CHECK_DELAY = 30 # 30 seconds

# Please override this with gsnap_chunk_timeout or rapsearch_chunk_timeout in DAG json.
# Default is several sigmas beyond the pale and indicates the data has to be QC-ed better.
# Note(2020-01-10): Raised to 3 hrs to mitigate Rapsearch chunk timeouts after recent index update.
DEFAULT_CHUNK_TIMEOUT = 60 * 60 * 3

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
    -t (threads): r5d.metal machines have 96 vCPUs. Use 48 threads and each process will be able to
    concurrently process 2 chunks (see attribute 'max_concurrent').

    For Rapsearch:
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
    def _service_inputs(host_filter_outputs):
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
        self.chunks_in_flight = threading.Semaphore(MAX_CHUNKS_IN_FLIGHT)
        self.chunks_result_dir_local = os.path.join(self.output_dir_local, "chunks")
        self.chunks_result_dir_s3 = os.path.join(self.output_dir_s3, "chunks")
        self.iostream_upload = multiprocessing.Semaphore(MAX_CONCURRENT_CHUNK_UPLOADS)
        command.make_dirs(self.chunks_result_dir_local)

    def count_reads(self):
        pass

    def run(self):
        ''' Run alignmment remotely '''

        service_inputs = PipelineStepRunAlignmentRemotely._service_inputs(self.input_files_local[0])
        cdhit_cluster_sizes_path, = self.input_files_local[1]
        output_m8, deduped_output_m8, output_hitsummary, output_counts_with_dcr_json = self.output_files_local()
        assert output_counts_with_dcr_json.endswith("_with_dcr.json"), self.output_files_local()

        service = self.additional_attributes["service"]
        self.run_remotely(service_inputs[service], output_m8, service)

        # get database
        lineage_db = fetch_reference(self.additional_files["lineage_db"], self.ref_dir_local)
        accession2taxid_db = fetch_reference(self.additional_files["accession2taxid_db"], self.ref_dir_local, allow_s3mi=True)

        min_alignment_length = NT_MIN_ALIGNMENT_LEN if service == 'gsnap' else 0
        m8.call_hits_m8(output_m8, lineage_db, accession2taxid_db,
                        deduped_output_m8, output_hitsummary, min_alignment_length)

        db_type = 'NT' if service == 'gsnap' else 'NR'
        evalue_type = 'log10' if service == 'rapsearch2' else 'raw'

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

    def run_remotely(self, input_fas, output_m8, service):
        sample_name = self.output_dir_s3.rstrip('/').replace('s3://', '').replace('/', '-')
        chunk_size = int(self.additional_attributes["chunk_size"])

        # Split files into chunks for performance
        part_suffix, input_chunks = self.chunk_input(input_fas, chunk_size)

        # Process chunks
        chunk_output_files = [None] * len(input_chunks)
        chunk_threads = []
        mutex = TraceLock("run_remotely", threading.RLock())
        # Randomize execution order for performance
        randomized = list(enumerate(input_chunks))
        random.shuffle(randomized)

        try:
            for n, chunk_input_files in randomized:
                self.chunks_in_flight.acquire()
                self.check_for_errors(mutex, chunk_output_files, input_chunks, service)
                t = threading.Thread(
                    target=PipelineStepRunAlignmentRemotely.run_chunk_wrapper,
                    args=[
                        self.chunks_in_flight, chunk_output_files, n, mutex, self.run_chunk,
                        [
                            part_suffix, chunk_input_files, service, True
                        ]
                    ])
                t.start()
                chunk_threads.append(t)

        finally:
            # Check chunk completion
            for ct in chunk_threads:
                ct.join()

        self.check_for_errors(mutex, chunk_output_files, input_chunks, service)

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
    def check_for_errors(mutex, chunk_output_files, input_chunks, service):
        with mutex:
            if "error" in chunk_output_files:
                # We already have per-chunk retries to address transient (non-deterministic) system issues.
                # If a chunk fails on all retries, it must be due to a deterministic & reproducible problem
                # with the chunk data or command options, so we should not even attempt the other chunks.
                err_i = chunk_output_files.index("error")
                raise RuntimeError("All retries failed for {} chunk {}.".format(
                    service, input_chunks[err_i]))

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


    def run_chunk(self, part_suffix, input_files, service, lazy_run):
        """
        Dispatch a chunk to worker machines for distributed GSNAP or RAPSearch
        group machines and handle their execution.
        """
        assert service in ("gsnap", "rapsearch2", "rapsearch")

        # TODO: (tmorse) standardize on rapsearch
        if service == "rapsearch2":
            service = "rapsearch"

        chunk_id = int(input_files[0].split(part_suffix)[-1])
        multihit_basename = f"multihit-{service}-out{part_suffix}{chunk_id}.m8"
        multihit_local_outfile = os.path.join(self.chunks_result_dir_local, multihit_basename)
        multihit_s3_outfile = os.path.join(self.chunks_result_dir_s3, multihit_basename)

        if lazy_run and fetch_from_s3(multihit_s3_outfile, multihit_local_outfile, okay_if_missing=True, allow_s3mi=False):
            log.write(f"finished alignment for chunk {chunk_id} with {service} by lazily fetching last result")
            return multihit_local_outfile

        # TODO: switch to environment variable
        environment = self.additional_attributes['environment']
        index_dir_suffix = self.additional_attributes['index_dir_suffix']
        # TODO: (tmorse) parameterize these
        priority_name = "normal"
        provisioning_model = "EC2"

        job_name = f"idseq-{service}-{environment}"
        job_queue = f"idseq-{service}-{environment}-{provisioning_model}-{index_dir_suffix}-{priority_name}"
        job_definition = f"idseq-{service}-{environment}"

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

        log.write(f"waiting for {service} batch queue for chunk {chunk_id}.")
        total_timeout = CHUNK_MAX_ATTEMPTS * CHUNK_ATTEMPT_TIMEOUT
        for _ in range((total_timeout // CHUNK_COMPLETE_CHECK_DELAY) + 1):
            status = PipelineStepRunAlignmentRemotely._get_job_status(client, job_id)
            if status == "SUCCEEDED":
                break
            if status == "FAILED":
                log.log_event("alignment_batch_error", values={"job_id": job_id})
                raise "chunk alignment failed"
            time.sleep(CHUNK_COMPLETE_CHECK_DELAY)

        fetch_from_s3(multihit_s3_outfile, multihit_local_outfile, okay_if_missing=True, allow_s3mi=False)

        log.write(f"finished alignment for chunk {chunk_id} on {service}")

        return multihit_local_outfile

    def step_description(self, require_docstrings=False):
        if (self.name == "gsnap_out"):
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
        elif (self.name == "rapsearch2_out"):
            return """
                Runs rapsearch remotely.

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
