import shutil

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.s3 import fetch_from_s3
import idseq_dag.util.command as command

import os
import glob

class PipelineStepRunSRST2(PipelineStep):
    '''
    srst2 is used to detect resistance genes.
    See: https://github.com/katholt/srst2
    '''
    def run(self):
        ''' Invoking srst2 '''
        print("AT START: CURR DIR")
        print(os.listdir()) 
        print("SELF.INPUT_FILES AT START")
        print((self.input_files))
        print("SELF.OUTPUT_DIR_LOCAL")
        print(self.output_dir_local)
        print(os.listdir(self.output_dir_local))
        print("SELF.GETCWD()")
        print(os.getcwd())
        # files = glob.glob(os.getcwd()+ '/*.txt')
        # for f in files:
        #    os.remove(f)
        # print("AFTER REMOVING ALL FILES: CURR DIR:")
        # print(os.listdir())
        db_file_path = fetch_from_s3(self.additional_files["resist_gene_db"], self.output_dir_local)
        min_cov = str(self.additional_attributes['min_cov'])
        n_threads = str(self.additional_attributes['n_threads'])   
        # print("BEFORE RUNNING TOOL: CURRENT DIRECTORY")
        # print(os.listdir())
        print("self.input_files_local")
        print(self.input_files_local)
        if len(self.input_files_local[0]) == 2:
            fwd_rd = self.input_files_local[0][0]
            rev_rd = self.input_files_local[0][1]
            # command.execute('cp '+ fwd_rd + ' R1_001.fastq.gz')
            # command.execute('cp ' + rev_rd + ' R2_001.fastq.gz')
            # print(os.listdir())
            srst2_params = [
                'srst2', '--input_pe', fwd_rd, rev_rd, '--forward', '001', '--reverse', '001',
                '--min_coverage', min_cov,'--threads', n_threads, '--output', self.output_dir_local+'/output',
                '--log', '--gene_db', db_file_path
            ]
        else:
            rd = self.input_files_local[0][0]
            shutil.copy2(rd, 'rd.fastq.gz')
            srst2_params = [
                'srst2', '--input_se', 'rd.fastq.gz', '--min_coverage', min_cov,'--threads', n_threads,
                '--output', self.output_dir_local, '--log', '--gene_db', db_file_path
            ]
        print("srst2 command: ")
        print(" ".join(srst2_params))
        command.execute(" ".join(srst2_params))
        # import os
        print("AFTER RUNNING TOOL: INPUT files")
        print(self.input_files_local)
        print("AFTER RUNNING TOOL: OUTPUT DIR")
        print(self.output_dir_local)
        print(os.listdir(self.output_dir_local))
        # command.execute('cp ' + output.log '+ self.output_files_local()[0])
        # command.execute('cp ' + self.output_dir_local + '/output__genes__ARGannot_r2__results.txt ' + self.output_files_local()[1])
    # Inherited method
    def count_reads(self):
        pass
