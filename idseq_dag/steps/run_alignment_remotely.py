import os
import threading
import shutil
import random
import traceback

from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
from idseq_dag.util.s3 import fetch_from_s3

class PipelineStepRunAlignmentRemotely(PipelineStep):
    '''
    Run gsnap/rapsearch2 remotely
    '''

    def __init__(self):
        PipelineStep.__init__(self)
        self.chunks_in_flight = threading.Semaphore(self.additional_attributes['chunks_in_flight'])
        self.chunks_result_dir_local = os.path.join(self.output_dir_local, "chunks")
        self.chunks_result_dir_s3 = os.path.join(self.output_dir_s3, "chunks")
        command.execute("mkdir -p %s" % self.chunks_result_dir_local)

    def run(self):
        ''' Run alignmment remotely '''
        input_fas = self.get_input_fas()
        output_m8 = self.output_files_local()[0]
        service = self.additional_attributes["service"]
        # TODO: run the alignment remotely and make lazy_chunk=True, revisit this later
        self.run_remotely(input_fas, output_m8, service)

    def run_remotely(input_fas, output_m8, service):
        assert service in ("gsnap", "rapsearch2")
        key_path = self.fetch_key(os.environ['KEY_PATH_S3'])
        sample_name = self.output_dir_s3.rstrip('/').replace('s3://', '').replace('/', '-')
        chunk_size = self.additional_attributes["chunk_size"]
        if service == "gsnap":
            remote_username = "ubuntu"
            remote_home_dir = "/home/%s" % remote_username
            remote_index_dir = "%s/share" % remote_home_dir
        elif service == "rapsearch2":
            remote_username = "ec2-user"
            remote_home_dir = "/home/%s" % remote_username
            remote_index_dir = "%s/references/nr_rapsearch" % remote_home_dir

        remote_work_dir = "%s/batch-pipeline-workdir/%s" % (remote_home_dir, sample_name)

        # Split files into chunks for performance
        part_suffix, input_chunks = self.chunk_input(input_fas, chunk_size)

        # Process chunks
        chunk_output_files = [None] * len(input_chunks)
        chunk_threads = []
        mutex = threading.RLock()
        # Randomize execution order for performance
        randomized = list(enumerate(input_chunks))
        random.shuffle(randomized)

        for n, chunk_input_files in randomized:
            self.chunks_in_flight.acquire()
            self.check_for_errors(mutex, chunk_output_files, input_chunks, service)
            t = threading.Thread(
                target=self.run_chunk_wrapper,
                args=[
                      chunks_in_flight, chunk_output_files, n, mutex, self.run_chunk, [
                      part_suffix, remote_home_dir, remote_index_dir,
                      remote_work_dir, remote_username, chunk_input_files,
                      key_path, service, True
                    ]
                ])
            t.start()
            chunk_threads.append(t)

        # Check chunk completion
        for ct in chunk_threads:
            ct.join()
            self.check_for_errors(mutex, chunk_output_files, input_chunks, service)
        assert None not in chunk_output_files
        # Concatenate the pieces and upload results
        self.concatenate_files(chunk_output_files, output_m8)

    def get_input_fas(self):
        service = self.additional_attributes["service"]
        if service == 'gsnap':
            return self.input_files_local[0][0:2]
        elif service = 'rapsearch2'
            return [self.input_files_local[0][2]]
        return None


    def fetch_key(self, key_path_s3):
        key_path = fetch_from_s3(key_path_s3, self.output_dir_local)
        command.execute("chmod 400 %s" % key_path)
        return key_path

    def chunk_input(self, input_files, chunksize):
        """Chunk input files into pieces for performance and parallelism."""
        part_lists = []  # Lists of partial files
        known_nlines = None
        part_suffix = ""
        chunk_nlines = chunksize * 2

        for input_file in input_files:
            # Count number of lines in the file
            nlines = int(command.execute_with_output("wc -l %s" % input_file)
                                .decode('utf-8').strip().split()[0])
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
            command.execute("split -a %d --numeric-suffixes -l %d %s %s" %
                            (ndigits, chunk_nlines, input_file, out_prefix))
            command.execute("aws s3 cp --quiet %s/ %s/ --recursive --exclude '*' --include '%s*'" %
                            (self.chunks_result_dir_local, self.chunks_result_dir_s3, out_prefix_base))

            # Get the partial file names
            partial_files = []
            paths = command.execute_with_output("ls %s*" % out_prefix)
                           .decode('utf-8').rstrip().split("\n")
            for pf in paths:
                partial_files.append(os.path.basename(pf))

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
                    task, input_chunks[err_i]))
    @staticmethod
    def concatenate_files(chunk_output_files, output_m8):
        with open(output_file, 'wb') as outf:
            for f in file_list:
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, outf)

    def run_chunk_wrapper(self, chunks_in_flight, chunk_output_files, n, mutex, target, args):
        result = "error"
        try:
            result = target(*args)
        except:
            with mutex:
                traceback.print_exc()
        finally:
            with mutex:
                chunk_output_files[n] = result
            chunks_in_flight.release()

    def run_chunk(self, part_suffix, remote_home_dir, remote_index_dir, remote_work_dir,
                  remote_username, input_files, key_path, service, lazy_run):
        """Dispatch a chunk to worker machines for distributed GSNAP or RAPSearch
        group machines and handle their execution.
        """
        assert service in ("gsnap", "rapsearch2")

        chunk_id = input_files[0].split(part_suffix)[-1]
        # TODO: Switch to python 3.6 which supports interpolation in string
        # formatting, and we will half the number of lines below.
        multihit_basename = "multihit-{service}-out{part_suffix}{chunk_id}.m8".format(
            service=service,
            part_suffix=part_suffix,
            chunk_id=chunk_id,
        )
        multihit_local_outfile = os.path.join(self.chunks_result_dir_local, multihit_basename)
        multihit_remote_outfile = os.path.join(remote_work_dir, multihit_basename)
        multihit_s3_outfile = os.path.join(self.chunks_result_dir_s3, multihit_basename)

        base_str = "aws s3 cp --quiet {s3_path}/{input_fa} {remote_work_dir}/{input_fa} "
        download_input_from_s3 = " ; ".join(
            base_str.format(
                s3_path=self.chunks_result_dir_s3,
                input_fa=input_fa,
                remote_work_dir=remote_work_dir) for input_fa in input_files)

        base_str = "mkdir -p {remote_work_dir} ; {download_input_from_s3} ; "
        if service == "gsnap":
            commands = base_str + "{remote_home_dir}/bin/gsnapl -A m8 --batch=0 --use-shared-memory=0 --gmap-mode=none --npaths=100 --ordered -t 32 --maxsearch=1000 --max-mismatches=40 -D {remote_index_dir} -d nt_k16 {remote_input_files} > {multihit_remote_outfile}"
        else:
            commands = base_str + "/usr/local/bin/rapsearch -d {remote_index_dir}/nr_rapsearch -e -6 -l 10 -a T -b 0 -v 50 -z 24 -q {remote_input_files} -o {multihit_remote_outfile}"

        commands = commands.format(
            remote_work_dir=remote_work_dir,
            download_input_from_s3=download_input_from_s3,
            remote_home_dir=remote_home_dir,
            remote_index_dir=remote_index_dir,
            remote_input_files=" ".join(
                remote_work_dir + "/" + input_fa for input_fa in input_files),
            multihit_remote_outfile=multihit_remote_outfile
            if service == "gsnap" else multihit_remote_outfile[:-3]
            # Strip the .m8 for RAPSearch as it adds that
        )

        if not lazy_run or not s3.fetch_from_s3(multihit_s3_outfile,
                                                multihit_local_outfile):
            correct_number_of_output_columns = 12
            min_column_number = 0
            max_tries = 2
            try_number = 1
            instance_ip = ""

            # Check if every row has correct number of columns (12) in the output
            # file on the remote machine
            while min_column_number != correct_number_of_output_columns \
                    and try_number <= max_tries:
                log.write("waiting for {} server for chunk {}".format(
                    service, chunk_id))
                max_concurrent = self.additional_attributes["max_concurrent"]
                environment = self.additional_attributes["environment"]

                instance_ip = self.wait_for_server_ip(service, key_path,
                                                 remote_username, environment,
                                                 max_concurrent, chunk_id)
                log.write("starting alignment for chunk %s on %s server %s" %
                             (chunk_id, service, instance_ip))
                command.execute(command.remote(commands, key_path, remote_username, instance_ip))

                if service == "gsnap":
                    verification_command = "cat %s" % multihit_remote_outfile
                else:
                    # For rapsearch, first remove header lines starting with '#'
                    verification_command = "grep -v '^#' %s" % multihit_remote_outfile
                verification_command += " | awk '{print NF}' | sort -nu | head -n 1"
                min_column_number_string = command.execute_with_output(
                    command.remote(verification_command, key_path, remote_username, instance_ip)).decode('utf-8')
                min_column_number = interpret_min_column_number_string(
                    min_column_number_string, correct_number_of_output_columns,
                    try_number)
                try_number += 1

            # Move output from remote machine to local machine
            msg = "Chunk %s output corrupt; not copying to S3. Re-start pipeline " \
                  "to try again." % chunk_id
            assert min_column_number == correct_number_of_output_columns, msg

            with iostream_uploads:  # Limit concurrent uploads so as not to stall the pipeline.
                with iostream:  # Still counts toward the general semaphore.
                    execute_command(
                        scp(key_path, remote_username, instance_ip,
                            multihit_remote_outfile, multihit_local_outfile))
                    execute_command("aws s3 cp --quiet %s %s/" %
                                    (multihit_local_outfile,
                                     SAMPLE_S3_OUTPUT_CHUNKS_PATH))
            write_to_log("finished alignment for chunk %s on %s server %s" %
                         (chunk_id, service, instance_ip))
        return multihit_local_outfile





