import json
import sys
import os
import threading
import time
import idseq_dag.util.command as command
import idseq_dag.util.log as log
from abc import abstractmethod
from enum import Enum, IntEnum

import idseq_dag.util.count as count
import idseq_dag.util.s3

class StepStatus(IntEnum):
    INITIALIZED = 0
    STARTED = 1 # step.start() called
    FINISHED = 2 # step.run() finished
    UPLOADED = 3 # all results uploaded to s3

class InputErrorType(Enum):
    ''' This error type will be used by the front-end to display a user-friendly error message '''
    NO_ERROR =  "NO_ERROR" # Placeholder value.
    INSUFFICIENT_READS = "INSUFFICIENT_READS"

class PipelineStep(object):
    ''' Each Pipeline Run Step i.e. run_star, run_bowtie2, etc '''
    def __init__(self, name, input_files, output_files,
                 output_dir_local, output_dir_s3, ref_dir_local,
                 additional_files, additional_attributes):
        ''' Set up all the input_files and output_files here '''
        self.name = name
        self.input_files = input_files # list of list files
        self.output_files = output_files # s3 location
        self.output_dir_local = output_dir_local
        self.output_dir_s3 = output_dir_s3.rstrip('/')
        self.ref_dir_local = ref_dir_local
        self.create_local_dirs()

        self.additional_files = additional_files
        self.additional_attributes = additional_attributes

        self.status = StepStatus.INITIALIZED
        self.exec_thread = None
        self.upload_thread = None
        self.input_files_local = []
        self.additional_files_to_upload = []
        self.optional_files_to_upload = []
        self.additional_folders_to_upload = []
        self.counts_dict = {}
        self.should_terminate = False
        self.should_count_reads = False

    @abstractmethod
    def run(self):
        ''' implement what is actually being run in this step '''

    @abstractmethod
    def count_reads(self):
        ''' count reads '''

    def stop_waiting(self):
        ''' stop waiting to run '''
        self.should_terminate = True

    def save_counts(self):
        if self.counts_dict:
            count_file_name = "%s/%s.count" % (self.output_dir_local, self.name)
            with open(count_file_name, 'w') as count_file:
                json.dump(self.counts_dict, count_file)
            self.additional_files_to_upload.append(count_file_name)

    def output_files_local(self):
        ''' Get list of output files on local folder '''
        return [os.path.join(self.output_dir_local, f) for f in self.output_files]

    def create_local_dirs(self):
        ''' make sure proper local directories are created for files with subdirs '''
        for f in self.output_files_local():
            command.execute("mkdir -p %s" % os.path.dirname(f))

    def uploading_results(self):
        ''' Upload output files to s3 '''
        files_to_upload = self.output_files_local() + self.additional_files_to_upload + [f for f in self.optional_files_to_upload if os.path.isfile(f)]
        for f in files_to_upload:
            # upload to S3 - TODO(Boris): parallelize the following with better calls
            s3_path = self.s3_path(f)
            idseq_dag.util.s3.upload_with_retries(f, s3_path)
        for f in self.additional_folders_to_upload:
            idseq_dag.util.s3.upload_folder_with_retries(f, self.s3_path(f))
        self.status = StepStatus.UPLOADED

    def s3_path(self, local_path):
        relative_path = os.path.relpath(local_path, self.output_dir_local)
        s3_path = os.path.join(self.output_dir_s3, relative_path)
        return s3_path

    @staticmethod
    def done_file(filename):
        ''' get the done file for a particular local file '''
        return "%s.done" % filename

    def wait_for_input_files(self):
        ''' wait for all the input files to be available and update input_files_local '''
        for fl in self.input_files:
            flist = []
            for f in fl:
                local_file = os.path.join(self.output_dir_local, f)
                while True:
                    if os.path.exists(local_file) and os.path.exists(self.done_file(local_file)):
                        flist.append(local_file)
                        break
                    else:
                        if self.should_terminate:
                            # If the step is not supposed to be run any more.
                            raise RuntimeError("Step %s being terminated" % self.name)
                        time.sleep(5)
            self.input_files_local.append(flist)

    def validate_input_files(self):
        ''' Validate input files before running the step '''
        errors = self.get_input_file_validation_errors()
        if errors and "error_type" in errors and errors["error_type"] != InputErrorType.NO_ERROR:
            log.write("Invalid input detected for step %s" % self.name)
            self.write_input_errors_json({
                "step": self.name,
                "errors": errors["errors"] if "errors" in errors else [],
                "error_type": errors["error_type"].name
            })
            raise RuntimeError("input files for step %s were not valid" % self.name)

    def get_input_file_validation_errors(self):
        ''' Return any input file validation errors'''
        # Steps should overwrite this method if necessary.
        return {
            # Step-specific error messages.
            "errors": [],
            # The broad type of error that occurred.
            "error_type": InputErrorType.NO_ERROR
        }

    def write_input_errors_json(self, error_json):
        ''' Write an input_file_errors.json file to the output_dir_s3 for this step, which can be detected by other services like idseq-web. '''
        log.write("Writing input_file_errors.json for step %s" % self.name)
        input_errors_file_basename = "input_file_errors.json"
        local_input_errors_file = "%s/%s" % (self.output_dir_local, input_errors_file_basename)
        s3_input_errors_file = "%s/%s" % (self.output_dir_s3, input_errors_file_basename)

        with open(local_input_errors_file, 'w') as input_errors_file:
            json.dump(error_json, input_errors_file)

        idseq_dag.util.s3.upload_with_retries(local_input_errors_file, s3_input_errors_file)

    @staticmethod
    def validate_input_files_min_reads(input_files, min_reads):
        '''
            Checks whether input files (fa, fasta, fq, fastq) have the minimum number of reads.
            Useful helper method that pipeline steps can use for validation.
        '''
        errors = []

        for index, input_file in enumerate(input_files):
            num_reads = count.reads(input_file)
            if num_reads < min_reads:
                errors.append("Input file %s requires at least %s reads (%s found)" % (index + 1, min_reads, num_reads))

        return errors

    def save_progress(self):
        ''' save progress after step run '''
        # save stats
        # save other info
        # TO BE IMPLEMENTED
        pass

    def validate(self):
        ''' Make sure all the output files are generated. '''
        for f in self.output_files_local():
            if not os.path.exists(f):
                raise RuntimeError("output file %s should be generated after run" % f)
            # Tag the done files
            done_file = self.done_file(f)
            command.execute("date > %s" % done_file)
        self.count_reads()

    def wait_until_finished(self):
        self.exec_thread.join()
        if self.status < StepStatus.FINISHED:
            raise RuntimeError("step %s run failed" % self.name)

    def wait_until_all_done(self):
        self.wait_until_finished()
        # run finished
        self.upload_thread.join()
        if self.status < StepStatus.UPLOADED:
            raise RuntimeError("step %s uploading failed" % self.name)

    def thread_run(self):
        ''' Actually running the step '''
        self.status = StepStatus.STARTED
        v = {"step": self.name}
        with log.log_context("dag_step", v):
            with log.log_context("substep_wait_for_input_files", v):
                self.wait_for_input_files()
            with log.log_context("validate_input_files", v):
                self.validate_input_files()
            with log.log_context("substep_run", v):
                self.run()
            with log.log_context("substep_validate", v):
                self.validate()
            with log.log_context("substep_save_progress", v):
                self.save_progress()
            with log.log_context("substep_save_counts", v):
                self.save_counts()
        self.upload_thread = threading.Thread(target=self.uploading_results)
        self.upload_thread.start()
        self.status = StepStatus.FINISHED

    def start(self):
        ''' function to be called after instantiation to start running the step '''
        self.exec_thread = threading.Thread(target=self.thread_run)
        self.exec_thread.start()

    def s3_path(self, local_path):
        relative_path = os.path.relpath(local_path, self.output_dir_local)
        s3_path = os.path.join(self.output_dir_s3, relative_path)
        return s3_path
