import os
import pandas as pd

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
            # TODO: Handle Fasta files
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
        # TODO: Handle case where local output file is named the same as that in dag and mv fails with "same file" error
        command.execute(f"mv {os.path.join(self.output_dir_local, 'output.log')} {self.output_files_local()[0]}")
        command.execute(f"mv {os.path.join(self.output_dir_local, 'output__genes__ARGannot_r2__results.txt')} {self.output_files_local()[1]}")
        if not os.path.exists(os.path.join(self.output_dir_local, 'output__fullgenes__ARGannot_r2__results.txt')): 
            # Cp to aws fails if output file is empty
            with open('empty_file.txt', 'a') as out:
                out.write('\n')
            # TODO: Have more efficient command than cp
            command.execute(f"cp empty_file.txt {self.output_files_local()[2]}")
            command.execute(f"cp empty_file.txt {self.output_files_local()[3]}")
            command.execute(f"cp empty_file.txt {self.output_files_local()[4]}")
        else:
            # Post processing of amr data
            # TODO: Handle case where local output file is named the same as that in dag and mv fails with "same file" error
            command.execute(f"mv {os.path.join(self.output_dir_local, 'output__fullgenes__ARGannot_r2__results.txt')} {self.output_files_local()[2]}")
            PipelineStepRunSRST2.generate_processed_amr_full_results(self.output_files_local()[2], self.output_files_local()[3])
            PipelineStepRunSRST2.generate_processed_amr_summary_results(self.output_files_local()[2], self.output_files_local()[4])

    # Inherited method
    def count_reads(self):
        pass
   
    @staticmethod 
    def generate_processed_amr_full_results(amr_raw_path, results_path):
        amr_processed_results = PipelineStepRunSRST2._process_amr_results(amr_raw_path)
        amr_processed_results.to_csv(
            'amr_processed_results.csv',
            mode='w',
            index=False,
            encoding='utf-8')
        command.execute(f"mv amr_processed_results.csv {results_path}")        
    
    @staticmethod
    def generate_processed_amr_summary_results(amr_raw_path, results_path):
        amr_summary_results = PipelineStepRunSRST2._summarize_amr_gene_families(amr_raw_path)
        amr_summary_results.to_csv(
            'amr_summary_results.csv',
            mode='w',
            index=False,
            encoding='utf-8')
        command.execute(f"mv amr_summary_results.csv {results_path}")

    @staticmethod
    def _get_pre_proc_amr_results(amr_raw_path):
        """ Reads in raw amr results file outputted by srst2, and does initial processing of marking gene family,
        + adding in blank columns for relevant metrics (total_genes_hit, total_coverage, total_depth)."""
        amr_results = pd.read_csv(amr_raw_path, delimiter="\t")
        # Parse out gene family as substring after '_', e.g. Aph_AGly's gene family would be AGly
        amr_results['gene_family'] = amr_results.apply(lambda row: row.gene[row.gene.index('_')+1:], axis=1)
        # Add empty columns for later processing
        amr_results['total_genes_hit'], amr_results['total_coverage'], amr_results['total_depth'] = (
            [amr_results.apply(lambda _: '', axis=1)]) * 3
        return amr_results

    @staticmethod
    def _summarize_amr_gene_families(amr_raw_path):
        """Returns a gene family-level summary of total_genes_hit, total_coverage, and total_depth."""
        amr_results = PipelineStepRunSRST2._get_pre_proc_amr_results(amr_raw_path)
        amr_results_by_gene_fam = amr_results.groupby(['gene_family'])
        total_genes_hit, total_coverage, total_depth = (amr_results_by_gene_fam.size().reset_index(name='total_genes_hit'),
                                                        amr_results_by_gene_fam['coverage'].sum().reset_index(name='total_coverage'),
                                                        amr_results_by_gene_fam['depth'].sum().reset_index(name='total_depth')
                                                       )
        amr_summary = pd.merge(pd.merge(total_genes_hit, total_coverage, on='gene_family', how='left'),
                               total_depth, on='gene_family', how='left')
        return amr_summary

    @staticmethod
    def _get_loc(df, row_i, col_name):
        """Helper function to get a pandas df location given a row index and col name."""
        return df.iloc[row_i, df.columns.get_loc(col_name)] 

    @staticmethod
    def _assign_locations(df1, row_idx1, df2, row_idx2, *col_names):
        """ Assigns values of col_names in df2 to the same col_names in df1.
            Assumes df1 and df2 contain the same col_names. """
        for col_name in col_names:
            df1.iloc[row_idx1, df1.columns.get_loc(col_name)] = PipelineStepRunSRST2._get_loc(df2, row_idx2, col_name)
    
    @staticmethod
    def _process_amr_results(amr_results_path):
        """ Returns processed amr result table with total_genes_hit, total_coverage,
            and total_depth column values filled in for all genes. """
        amr_summary = PipelineStepRunSRST2._summarize_amr_gene_families(amr_results_path)
        amr_results = PipelineStepRunSRST2._get_pre_proc_amr_results(amr_results_path)
        # Sorting reduces run time when amalgamating amr_results and amr_summary
        sorted_amr = amr_results.sort_values(by=['gene_family'])
        i_srt_amr_row = 0
        for i_amr_smry_row in range(amr_summary.shape[0]):
            gene_fam = PipelineStepRunSRST2._get_loc(amr_summary, i_amr_smry_row, 'gene_family')
            while i_srt_amr_row < sorted_amr.shape[0] and (
                PipelineStepRunSRST2._get_loc(sorted_amr, i_srt_amr_row, 'gene_family') == gene_fam):
                PipelineStepRunSRST2._assign_locations(sorted_amr,
                                          i_srt_amr_row,
                                          amr_summary,
                                          i_amr_smry_row,
                                          'total_genes_hit',
                                          'total_coverage',
                                          'total_depth')
                i_srt_amr_row += 1
        return sorted_amr    


