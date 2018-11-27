''' Generate Coverage Statistics '''
import json
import pysam

from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.count as count

class PipelineStepGenerateCoverageStats(PipelineStep):
    ''' Generate Coverage Statistics from Assembly Output '''
    def run(self):
        """
          1. extract contigs.fasta and read-contig.sam
          2. run pile up
        """
        _contigs, _scaffolds, read_contig_sam, _stats = self.input_files_local[0]
        coverage_json, coverage_summary_csv = self.output_files_local()

        # generate bam files
        bam_file = read_contig_sam.replace(".sam", ".bam")
        command.execute(f"samtools view -S -b  {read_contig_sam} | samtools sort - -o {bam_file}")
        command.execute(f"samtools index {bam_file}")

        # run coverage info
        contig2coverage = {}
        with pysam.AlignmentFile(bam_file, "rb") as f:
            contig_names = f.references
            for c in contig_names:
                coverage = []
                for pileup_column in f.pileup(contig=c):
                    coverage.append(pileup_column.get_num_aligned())
                sorted_coverage = sorted(coverage)
                min_cov = sorted_coverage[:10]
                max_cov = sorted_coverage[-10:]
                contig2coverage[c] = {
                    "coverage": coverage,
                    "min_cov": min_cov,
                    "max_cov": max_cov
                }

        with open(coverage_json, 'w') as csf:
            json.dump(contig2coverage, csf)

        with open(coverage_summary_csv, 'w') as csc:
            csc.write("contig_name,min_cov_1,min_cov_2,min_cov_3,max_cov_1,maxx_cov_2,max_cov3\n")
            for contig, stats in contig2coverage.items():
                output_row = [contig] + stats['min_cov'][:3] + stats['max_cov'][-1:-4:-1]
                output_str = ','.join(str(e) for e in output_row)
                csc.write(output_str + "\n")

    def count_reads(self):
        ''' Count reads '''
        pass


