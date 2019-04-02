''' Generate loc db  '''
import dbm
import re
from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.command import run_in_subprocess
import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.count as count
from idseq_dag.util.dict import IdSeqDict, IdSeqDictValue
from idseq_dag.steps.generate_gsnap_index import PipelineStepGenerateGsnapIndex
BATCH_INSERT_SIZE = 300

class PipelineStepGenerateRRnaDB(PipelineStep):
    ''' Generate rRNA database from NT and build gsnap index  '''
    def run(self):
        """
        Main function
        """
        nt_file = self.input_files_local[0][0]
        output_fasta = self.output_files_local()[0]
        output_gsnap_index = self.output_files_local()[1]
        k = self.additional_attributes.get("k", 16) # kmer k
        self.generate_rrna_db(nt_file, output_fasta, k, output_gsnap_index)

    @run_in_subprocess
    def generate_rrna_db(self, nt_file, output_fasta, k, output_gsnap_index):
        outf = open(output_fasta, 'w', encoding='utf-8')
        with open(nt_file) as dbf:
            lines = 0
            output_current_seq = False
            for line in dbf:
                lines += 1
                if lines % 100000 == 0:
                    log.write(f"{lines/1000000.0}M lines")
                if line[0] == '>':  # header line
                    output_current_seq = False
                    if re.search("rrna", line, flags=re.IGNORECASE):
                        outf.write(line)
                        output_current_seq = True
                else:
                    if output_current_seq:
                        outf.write(line)
        outf.close()
        PipelineStepGenerateGsnapIndex.generate_gsnap_index(output_fasta, k, output_gsnap_index)



    def count_reads(self):
        ''' Count reads '''
        pass

