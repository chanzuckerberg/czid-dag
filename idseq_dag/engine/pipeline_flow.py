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
        self.nodes = dag["nodes"]
        self.steps = dag["steps"]
        self.head_nodes = dag["head_nodes"]
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
        for f in self.large_file_list:
            idseq_dag.util.s3.fetch_from_s3(f, self.ref_dir_local, allow_s3mi=True)

    @staticmethod
    def parse_and_validate_conf(dag_json):
        '''
        Validate the json format. see examples/*.json.
        Required fields are:
          "output_dir_s3": base results folder. a pipeline version number will be appended for real output folder.
          "nodes": lists of files that are given or would be generated
          "steps": steps that species actions to generate input and output
          "head_nodes": input files that are given

        '''
        dag = json.loads(open(dag_json).read())
        output_dir = dag["output_dir_s3"]
        nodes = dag["nodes"]
        steps = dag["steps"]
        head_nodes = dag["head_nodes"]
        covered_nodes = set()
        for s in steps:
            # validate each step in/out are valid nodes
            for inode in s["in"]:
                if inode not in nodes:
                    raise ValueError("input %s doesn't exit for step %s" % (inode, s["out"]))
            if s["out"] not in nodes:
                raise ValueError("%s node doesn't exit" % s["out"])
            if s["out"] in covered_nodes:
                raise ValueError("%s hasn't been generated in other steps" % s["out"])
            covered_nodes.add(s["out"])
        for node_name, node_data in head_nodes.items():
            # validate the head nodes exist in s3
            s3_path = node_data["s3_dir"]
            covered_nodes.add(node_name)
            for file_name in nodes[node_name]:
                s3_file = os.path.join(s3_path, file_name)
                if not idseq_dag.util.s3.check_s3_presence(s3_file):
                    raise ValueError("%s file doesn't exist" % s3_file)
        # Check that all nodes are covered
        # ALL Inputs Outputs VALIDATED
        for node_name in nodes.keys():
            if node_name not in covered_nodes:
                raise ValueError("%s couldn't be generated from the steps" % node_name)

        return dag

    def plan(self):
        '''
            Traverse through the nodes and steps and calculate
            1. the large file download priority based on the how deep the step is
            2. if a step needs to be run based on the existence of output file and lazy run parameter
        '''
        covered_nodes = {}
        large_file_download_list = []
        step_list = []
        for node_name in self.head_nodes.keys():
            covered_nodes[node_name] = { 'depth': 0, 'lazy_run': self.lazy_run, 's3_downloadable': True }
        steps_complete = set()
        while len(steps_complete) < len(self.steps):
            # run until all the steps can be run
            current_nodes = {}
            for step in self.steps:
                if step["out"] not in steps_complete:
                    step_can_be_run = True
                    depth_max = 0
                    lazy_run = True
                    for node in step["in"]:
                        if node not in covered_nodes:
                            step_can_be_run = False
                            break
                        else:
                            depth_max = max(covered_nodes[node]['depth'], depth_max)
                            if covered_nodes[node]['lazy_run'] == False:
                                lazy_run = False
                    if step_can_be_run: # All the input is satisfied
                        steps_complete.add(step["out"])
                        file_list= self.nodes[step["out"]]
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
                        # update nodes available for the next round
                        current_nodes[step["out"]] = { 'depth': (depth_max + 1), 'lazy_run': lazy_run, 's3_downloadable': s3_downloadable}
            covered_nodes.update(current_nodes)
        return (step_list, large_file_download_list, covered_nodes)

    @staticmethod
    def fetch_input_files_from_s3(input_files, input_dir_s3, result_dir_local):
        for f in input_files:
            s3_file = os.path.join(input_dir_s3, f)
            local_file = os.path.join(result_dir_local, f)
            # copy the file over
            idseq_dag.util.s3.fetch_from_s3(s3_file, result_dir_local, allow_s3mi=True)
            # write the done_file
            done_file = PipelineStep.done_file(local_file)
            command.execute("date > %s" % done_file)

    def fetch_node_from_s3(self, node):
        ''' .done file should be written to the result dir when the download is complete '''
        log.write("Downloading node %s" % node)
        if node in self.head_nodes:
            input_path_s3 = self.head_nodes[node]["s3_dir"]
        else:
            input_path_s3 = self.output_dir_s3

        PipelineFlow.fetch_input_files_from_s3(input_files=self.nodes[node],
                                               input_dir_s3=input_path_s3,
                                               result_dir_local=self.output_dir_local)


    def start(self):
        # Come up with the plan
        (step_list, self.large_file_list, covered_nodes) = self.plan()

        for step in step_list: # download the files from s3 when necessary
            for node in step["in"]:
                node_info = covered_nodes[node]
                if node_info['s3_downloadable']:
                    threading.Thread(target=self.fetch_node_from_s3, args=(node,)).start()


        # TODO(boris): check the following implementation
        threading.Thread(target=self.prefetch_large_files)


        # Start initializing all the steps and start running them and wait until all of them are done
        step_instances = []
        for step in step_list:
            log.write("Starting step %s" % step["out"])
            StepClass = getattr(importlib.import_module(step["module"]), step["class"])
            step_output = self.nodes[step["out"]]
            step_inputs = [self.nodes[inode] for inode in step["in"]]
            step_instance = StepClass(step["out"], step_inputs, step_output,
                                      self.output_dir_local, self.output_dir_s3, self.ref_dir_local,
                                      step["additional_files"], step["additional_attributes"])
            step_instance.start()
            step_instances.append(step_instance)
        # Collecting stats files
        for step in step_instances:
            step.wait_until_all_done()
        log.write("all steps are done")

