import idseq_dag.util.command as command
from idseq_dag.engine.pipeline_step import PipelineStep

class CountLines(PipelineStep):
    def run(self):
        result_file = self.output_files_local()
        line_count = sum([int(command.execute_with_output("wc -l %s" % file).strip().split()[0]) for file in
                          self.input_files_local[0]])
        print ("RES " + str(result_file) + "\n")
        with open(result_file[0], 'w') as out:
            out.write(str(line_count) + "\n")
