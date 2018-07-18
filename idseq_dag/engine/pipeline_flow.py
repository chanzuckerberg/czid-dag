import importlib
import json
import sys
import os
import threading
import traceback

import idseq_dag
import idseq_dag.util.s3
import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.count as count
from idseq_dag.engine.pipeline_step import PipelineStep

DEFAULT_OUTPUT_DIR_LOCAL = '/mnt/idseq/results/%d' % os.getpid()
DEFAULT_REF_DIR_LOCAL = '/mnt/idseq/ref'

class PipelineFlow(object):
    def __init__(self, lazy_run, dag_json):
        '''
            See examples/example_dag.json and
                idseq_dag.main.validate_dag_json for more details.
        '''
        self.lazy_run = lazy_run
        dag = PipelineFlow.parse_and_validate_conf(dag_json)
        self.targets = dag["targets"]
        self.steps = dag["steps"]
        self.given_targets = dag["given_targets"]
        self.output_dir_s3 = os.path.join(dag["output_dir_s3"],
                                          self.parse_output_version(idseq_dag.__version__))
        self.output_dir_local = dag.get("output_dir_local", DEFAULT_OUTPUT_DIR_LOCAL).rstrip('/')
        self.ref_dir_local = dag.get("ref_dir_local", DEFAULT_REF_DIR_LOCAL)
        self.large_file_list = []

        command.execute("mkdir -p %s %s" % (self.output_dir_local, self.ref_dir_local))

    @staticmethod
    def parse_output_version(version):
        return ".".join(version.split(".")[0:2])

    def prefetch_large_files(self):
        log.write("downloading large files: %s" % ",".join(self.large_file_list))
        for f in self.large_file_list:
            log.write("downloading %s" % f)
            idseq_dag.util.s3.fetch_from_s3(f, self.ref_dir_local, allow_s3mi=True, auto_untar=True)

    @staticmethod
    def parse_and_validate_conf(dag_json):
        '''
        Validate the json format. see examples/*.json.
        Required fields are:
          "output_dir_s3": base results folder. a pipeline version number will be appended for real output folder.
          "targets": lists of files that are given or would be generated
          "steps": steps that species actions to generate input and output
          "given_targets": input files that are given

        '''
        dag = json.loads(open(dag_json).read())
        output_dir = dag["output_dir_s3"]
        targets = dag["targets"]
        steps = dag["steps"]
        given_targets = dag["given_targets"]
        covered_targets = set()
        for s in steps:
            # validate each step in/out are valid targets
            for itarget in s["in"]:
                if itarget not in targets:
                    raise ValueError("input %s doesn't exit for step %s" % (itarget, s["out"]))
            if s["out"] not in targets:
                raise ValueError("%s target doesn't exit" % s["out"])
            if s["out"] in covered_targets:
                raise ValueError("%s hasn't been generated in other steps" % s["out"])
            covered_targets.add(s["out"])
        for target_name, target_data in given_targets.items():
            # validate the given targets exist in s3
            s3_path = target_data["s3_dir"]
            covered_targets.add(target_name)
            for file_name in targets[target_name]:
                s3_file = os.path.join(s3_path, file_name)
                if not idseq_dag.util.s3.check_s3_presence(s3_file):
                    raise ValueError("%s file doesn't exist" % s3_file)
        # Check that all targets are covered
        # ALL Inputs Outputs VALIDATED
        for target_name in targets.keys():
            if target_name not in covered_targets:
                raise ValueError("%s couldn't be generated from the steps" % target_name)

        return dag

    def plan(self):
        '''
            Traverse through the targets and steps and calculate
            1. the large file download priority based on the how deep the step is
            2. if a step needs to be run based on the existence of output file and lazy run parameter
        '''
        covered_targets = {}
        large_file_download_list = []
        step_list = []
        for target_name in self.given_targets.keys():
            covered_targets[target_name] = { 'depth': 0, 'lazy_run': self.lazy_run, 's3_downloadable': True }
        steps_complete = set()
        while len(steps_complete) < len(self.steps):
            # run until all the steps can be run
            current_targets = {}
            for step in self.steps:
                if step["out"] not in steps_complete:
                    step_can_be_run = True
                    depth_max = 0
                    lazy_run = True
                    for target in step["in"]:
                        if target not in covered_targets:
                            step_can_be_run = False
                            break
                        else:
                            depth_max = max(covered_targets[target]['depth'], depth_max)
                            if covered_targets[target]['lazy_run'] == False:
                                lazy_run = False
                    if step_can_be_run: # All the input is satisfied
                        steps_complete.add(step["out"])
                        file_list= self.targets[step["out"]]
                        if lazy_run and idseq_dag.util.s3.check_s3_presence_for_file_list(self.output_dir_s3, file_list):
                            # output can be lazily generated. touch the output
                            #idseq_dag.util.s3.touch_s3_file_list(self.output_dir_s3, file_list)
                            s3_downloadable = True
                        else:
                            # steps need to be run
                            lazy_run = False
                            s3_downloadable = False
                            step_list.append(step)
                            # The following can be changed to append if we want to get the round information
                            large_file_download_list += step["additional_files"].values()
                        # update targets available for the next round
                        current_targets[step["out"]] = { 'depth': (depth_max + 1), 'lazy_run': lazy_run, 's3_downloadable': s3_downloadable}
            covered_targets.update(current_targets)
        return (step_list, large_file_download_list, covered_targets)

    @staticmethod
    def fetch_input_files_from_s3(input_files, input_dir_s3, result_dir_local):
        for f in input_files:
            s3_file = os.path.join(input_dir_s3, f)
            local_file = os.path.join(result_dir_local, f)
            local_dir = os.path.dirname(local_file)
            command.execute("mkdir -p %s" % local_dir)
            # copy the file over
            output_file = idseq_dag.util.s3.fetch_from_s3(s3_file, local_dir, allow_s3mi=True)
            if output_file:
                # write the done_file
                done_file = PipelineStep.done_file(local_file)
                command.execute("date > %s" % done_file)
            else:
                raise RuntimeError(f"{s3_file} likely doesn't exist")

    @staticmethod
    def count_input_reads(input_files, result_dir_local, result_dir_s3, target_name, max_fragments=None):
        local_input_files = [os.path.join(result_dir_local, f) for f in input_files[0:2]]
        count_file_basename = "%s.count" % target_name
        local_count_file = "%s/%s" % (result_dir_local, count_file_basename)
        s3_count_file = "%s/%s" % (result_dir_s3, count_file_basename)

        read_count = count.reads_in_group(local_input_files, max_fragments=max_fragments)
        counts_dict = { target_name: read_count }
        if read_count == len(local_input_files) * max_fragments:
            # If the number of reads is exactly equal to the maximum we specified,
            # it means that the input has been truncated.
            counts_dict["truncated"] = read_count

        with open(local_count_file, 'w') as count_file:
            json.dump(counts_dict, count_file)
        idseq_dag.util.s3.upload_with_retries(local_count_file, s3_count_file)


    def fetch_target_from_s3(self, target):
        ''' .done file should be written to the result dir when the download is complete '''
        log.write("Downloading target %s" % target)
        if target in self.given_targets:
            input_path_s3 = self.given_targets[target]["s3_dir"]
        else:
            input_path_s3 = self.output_dir_s3

        PipelineFlow.fetch_input_files_from_s3(input_files=self.targets[target],
                                               input_dir_s3=input_path_s3,
                                               result_dir_local=self.output_dir_local)
        if target in self.given_targets and self.given_targets[target].get("count_reads"):
            PipelineFlow.count_input_reads(input_files=self.targets[target],
                                           result_dir_local=self.output_dir_local,
                                           result_dir_s3=self.output_dir_s3,
                                           target_name=target,
                                           max_fragments=self.given_targets[target]["max_fragments"])


    def start(self):
        # Come up with the plan
        (step_list, self.large_file_list, covered_targets) = self.plan()

        for step in step_list: # download the files from s3 when necessary
            for target in step["in"]:
                target_info = covered_targets[target]
                if target_info['s3_downloadable']:
                    threading.Thread(target=self.fetch_target_from_s3, args=(target,)).start()

        # TODO(boris): check the following implementation
        threading.Thread(target=self.prefetch_large_files).start()


        # Start initializing all the steps and start running them and wait until all of them are done
        step_instances = []
        for step in step_list:
            log.write("Initializing step %s" % step["out"])
            StepClass = getattr(importlib.import_module(step["module"]), step["class"])
            step_output = self.targets[step["out"]]
            step_inputs = [self.targets[itarget] for itarget in step["in"]]
            step_instance = StepClass(step["out"], step_inputs, step_output,
                                      self.output_dir_local, self.output_dir_s3, self.ref_dir_local,
                                      step["additional_files"], step["additional_attributes"])
            step_instance.start()
            step_instances.append(step_instance)
        # Collecting stats files
        for step in step_instances:
            try:
                step.wait_until_all_done()
            except:
                # Some exception thrown by one of the steps
                traceback.print_exc()
                for s in step_instances:
                    # notify the waiting step instances to self destruct
                    s.stop_waiting()
        log.write("all steps are done")

