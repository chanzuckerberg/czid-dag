from idseq_dag.engine.pipeline_step import PipelineStep, InputErrorType
import idseq_dag.util.command as command
import idseq_dag.util.count as count

class PipelineStepRunCDHitDup(PipelineStep):
    '''
    CD-HIT-DUP is used to identify duplicates from single or paired reads.
    Two FASTA inputs means paired reads.
    See: http://weizhongli-lab.org/cd-hit/
    '''
    def get_input_file_validation_errors(self):
        # Return errors if either input file has less than 2 reads.
        errors = PipelineStep.validate_input_files_min_reads(self.input_files_local[0], 2)

        if errors:
            return {
                "errors": errors,
                "error_type": InputErrorType.INSUFFICIENT_READS
            }

        return None

    def run(self):
        ''' Invoking cd-hit-dup '''
        input_fas = self.input_files_local[0]
        output_fas = self.output_files_local()
        cdhitdup_params = [
            'cd-hit-dup', '-i', input_fas[0], '-o', output_fas[0],
            '-e', '0.05', '-u', '70'
        ]
        if len(input_fas) == 2:
            cdhitdup_params += ['-i2', input_fas[1], '-o2', output_fas[1]]
        command.execute(" ".join(cdhitdup_params))

    def count_reads(self):
        self.should_count_reads = True
        self.counts_dict[self.name] = count.reads_in_group(self.output_files_local()[0:2])

