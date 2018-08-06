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
        log = os.path.join(self.output_dir_local, 'output.log')
        log_dest = self.output_files_local()[0]
        results = os.path.join(self.output_dir_local, 'output__genes__ARGannot_r2__results.txt')
        results_dest = self.output_files_local()[1]
        os.rename(log, log_dest)
        os.rename(results, results_dest) 
        if not os.path.exists(os.path.join(self.output_dir_local, 'output__fullgenes__ARGannot_r2__results.txt')): 
            PipelineStepRunSRST2.fill_file_path(self.output_files_local()[2])
            PipelineStepRunSRST2.fill_file_path(self.output_files_local()[3])
            PipelineStepRunSRST2.fill_file_path(self.output_files_local()[4])
        else:
            # Post processing of amr data
            results_full = os.path.join(self.output_dir_local, 'output__fullgenes__ARGannot_r2__results.txt')
            results_full_dest = self.output_files_local()[2]
            PipelineStepRunSRST2.mv_to_dest(results_full, results_full_dest)
            self.process_amr_results(results_full_dest)

    # Inherited method
    def count_reads(self):
        pass

    @staticmethod
    def fill_file_path(file_path):
        """Helper function to open an "empty" file at a given file location.
           Note that aws s3 cannot upload a 0 byte file from local to s3;
           see  https://github.com/aws/aws-cli/issues/2403.
           Possible problems seem to be from python versions mismatching on aws cli and linux
           system. This doesn't seem to be the case on the staging machine, though.
           So the official recommendation to install aws-cli from pip does not seem to apply since
           the Python versions match up.
           I tried many suggestions on the link above + others -- for now, the following
           seems a reasonable workaround.
           We use os file functions to avoid overhead of using Python's File object
           functions to open an empty file  and using a cp command to move to destination.
           TODO: See if there are aws CLI installation errors that can be fixed."""
        fd = os.open(file_path, os.O_RDWR|os.O_CREAT)
        os.write(fd, b"\n")
        os.close(fd)
             
    @staticmethod
    def _get_pre_proc_amr_results(amr_raw_path):
        """ Reads in raw amr results file outputted by srst2, and does initial processing of marking gene family."""
        amr_results = pd.read_csv(amr_raw_path, delimiter="\t")
        # Parse out gene family as substring after '_', e.g. Aph_AGly's gene family would be AGly
        amr_results['gene_family'] = amr_results.apply(lambda row: row.gene.split('_', 1)[1], axis=1)
        return amr_results

    @staticmethod
    def _summarize_amr_gene_families(amr_results):
        """Returns a gene family-level summary of total_genes_hit, total_coverage, and total_depth."""
        amr_summary = amr_results.groupby(['gene_family']).agg({'gene_family': ['size'],'coverage':['sum'], 'depth':['sum']})
        amr_summary.columns = [' '.join(col) for col in amr_summary.columns]
        amr_summary = amr_summary.rename(columns={'gene_family size': 'total_gene_hits', 'coverage sum': 'total_coverage', 'depth sum': 'total_depth'}).reset_index()
        return amr_summary


    def process_amr_results(self, amr_results_path):
        """ Writes processed amr result table with total_genes_hit, total_coverage,
            and total_depth column values filled in for all genes to output files; and likewise with
            the gene-family level summary of amr results table. """
        amr_results = PipelineStepRunSRST2._get_pre_proc_amr_results(amr_results_path)
        amr_summary = PipelineStepRunSRST2._summarize_amr_gene_families(amr_results)
        amr_summary.to_csv(
            self.output_files_local()[4],
            mode='w',
            index=False,
            encoding='utf-8')
        sorted_amr = amr_results.sort_values(by=['gene_family'])
        proc_amr = pd.merge_ordered(sorted_amr, amr_summary, fill_method='ffill', left_by=['gene_family'])
        proc_amr.to_csv(
            self.output_files_local()[3],
            mode='w',
            index=False,
            encoding='utf-8')


