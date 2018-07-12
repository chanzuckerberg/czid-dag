import shutil

from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command

class PipelineStepRunSRST2(PipelineStep):
    '''
    srst2 is used to detect resistance genes.
    See: https://github.com/katholt/srst2
    '''
    def run(self):
        ''' Invoking srst2 '''
        fwd_rd = self.input_files_local[0][0]
        rev_rd = self.input_files_local[0][1]
        # TODO: pass in database in additional_files field 
        db_file_path = self.input_files_local[1][0]
        db_name = "ARGannot_r2" # hardcode for now
        # srst2 assumes database file is in current directory
        shutil.copy2(db_file_path, db_name)
        min_cov = str(self.additional_attributes['min_cov'])
        n_threads = str(self.additional_attributes['n_threads'])
        output_prefix = self.additional_attributes['output_prefix']   
        # Note that we pass in 001 assuming reads are in standard fastq.gz format.
        srst2_params = [
            'srst2', '--input_pe', fwd_rd, rev_rd, '--forward', '001', '--reverse', '001',
            '--min_coverage', min_cov,'--threads', n_threads, '--output', output_prefix,
            '--log', '--gene_db', db_file_path
        ]
        command.execute(" ".join(srst2_params))
        shutil.copy2('output__genes__ARGannot_r2__results.txt', self.output_files_local()[0])

    # Inherited method
    def count_reads(self):
        pass
