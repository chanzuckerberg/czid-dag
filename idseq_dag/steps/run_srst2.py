from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import os
import shutil

class PipelineStepRunSRST2(PipelineStep):
    '''
    srst2 is used to detect resistance genes.
    See: https://github.com/katholt/srst2
    '''
    def run(self):
        ''' Invoking srst2 '''
        print("-------SELF.INPUT FILES SRST2----")
        print(self.input_files_local)
        fwd_rd = self.input_files_local[0][0]
        rev_rd = self.input_files_local[0][1]
        # Note: this is different from db_name, e.g. ARGannot_r2. Consider string splicing db_path to get db_name.
        # print("---------SELF.ADDTIONAL_FILES-----")
        # print(self.additional_files)
        # db_file_path = self.additional_files['resist_gene_db'] 
        db_file_path = self.input_files_local[1][0]
        db_name = "ARGannot_r2" # hardcode for now
        shutil.copy2(db_file_path, db_name)
        print("----NEW DB FILE PATH-----")
        print(db_file_path)
        min_cov = str(self.additional_attributes['min_cov'])
        n_threads = str(self.additional_attributes['n_threads'])
        output_prefix = self.additional_attributes['output_prefix']
        # db_name = self.additional_attributes['db_name']   
        # Note that we pass in 001 assuming reads are in standard fastq.gz format.
        srst2_params = [
            'srst2', '--input_pe', fwd_rd, rev_rd, '--forward', '001', '--reverse', '001',
            '--min_coverage', min_cov,'--threads', n_threads, '--output', output_prefix,
            '--log', '--gene_db', db_file_path
        ]
        command.execute(" ".join(srst2_params))
        print(os.listdir())
        print('-----SELF.OUTPUT_FILES_LOCAL-----')
        print(self.output_files_local())
        print('------CWD-----')
        print(os.getcwd())
        shutil.copy2('output__genes__ARGannot_r2__results.txt', self.output_files_local()[0])
    # Inherited method
    def count_reads(self):
        pass
