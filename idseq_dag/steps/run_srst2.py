import os

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.s3 import fetch_from_s3
import idseq_dag.util.command as command

class PipelineStepRunSRST2(PipelineStep):
    '''
    srst2 is used to detect resistance genes.
    See: https://github.com/katholt/srst2
    '''
    def run(self):
        ''' Invoking srst2 '''
        db_file_path = fetch_from_s3(self.additional_files["resist_gene_db"], self.output_dir_local)
        min_cov = str(self.additional_attributes['min_cov'])
        n_threads = str(self.additional_attributes['n_threads'])
        if len(self.input_files_local[0]) == 2:
            fwd_rd = self.input_files_local[0][0]
            rev_rd = self.input_files_local[0][1]
            command.execute('mv ' + fwd_rd + ' R1_001.fastq.gz')
            command.execute('mv ' + rev_rd + ' R2_001.fastq.gz')
            srst2_params = [
                'srst2', '--input_pe', 'R1_001.fastq.gz', 'R2_001.fastq.gz', '--forward', '001', '--reverse', '001',
                '--min_coverage', min_cov,'--threads', n_threads, '--output', self.output_dir_local+'/output',
                '--log', '--gene_db', db_file_path
            ]
        else:
            rd = self.input_files_local[0][0]
            command.execute('mv '+ rd + ' rd.fastq.gz')
            srst2_params = [
                'srst2', '--input_se', 'rd.fastq.gz', '--min_coverage', min_cov,'--threads', n_threads,
                '--output', self.output_dir_local+'/output', '--log', '--gene_db', db_file_path
            ]
        command.execute(" ".join(srst2_params))
        if 'output__fullgenes__ARGannot_r2__results.txt' not in os.listdir(self.output_dir_local): 
            # Cp to aws fails if output file is empty
            with open('output__fullgenes__ARGannot_r2__results.txt', 'a') as out:
                out.write('\n')
            command.execute('mv output__fullgenes__ARGannot_r2__results.txt ' + self.output_files_local()[2])
    
    # Inherited method
    def count_reads(self):
        pass
