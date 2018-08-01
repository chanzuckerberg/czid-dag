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
            # Rename to format srst2 expects
            command.execute(f"ln -sf {fwd_rd} _R1_001.fastq.gz")
            command.execute(f"ln -sf {rev_rd} _R2_001.fastq.gz")
            srst2_params = [
                'srst2', '--input_pe', '_R1_001.fastq.gz', '_R2_001.fastq.gz', '--forward', '_R1_001', '--reverse', '_R2_001',
                '--min_coverage', min_cov,'--threads', n_threads, '--output', os.path.join(self.output_dir_local, 'output'),
                '--log', '--gene_db', db_file_path
            ]
        else:
            rd = self.input_files_local[0][0]
            command.execute(f"ln -sf {rd} rd.fastq.gz")
            srst2_params = [
                'srst2', '--input_se', 'rd.fastq.gz', '--min_coverage', min_cov,'--threads', n_threads,
                '--output', self.output_dir_local+'/output', '--log', '--gene_db', db_file_path
            ]
        command.execute(" ".join(srst2_params))
        if not os.path.exists(os.path.join(self.output_dir_local, 'output__fullgenes__ARGannot_r2__results.txt')): 
            # Cp to aws fails if output file is empty
            with open('output__fullgenes__ARGannot_r2__results.txt', 'a') as out:
                out.write('\n')
            command.execute(f"mv output__fullgenes__ARGannot_r2__results.txt {self.output_files_local()[2]}")
    
    # Inherited method
    def count_reads(self):
        pass
