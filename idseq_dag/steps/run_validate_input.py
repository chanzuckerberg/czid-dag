import json

from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.count as count

class PipelineStepRunValidateInput(PipelineStep):
    def run(self):
        # Setup
        input_files = self.input_files_local[0][0:2]
        num_inputs = len(input_files)
        assert num_inputs == 1 or num_inputs == 2, 'Invalid number of input files'
        output_files = self.output_files_local()[1:3]
        summary_file = self.output_files_local()[0]
        max_fragments = self.additional_attributes["truncate_fragments_to"]

        file_ext = self.additional_attributes.get("file_ext")
        assert file_ext == 'fastq' or file_ext == 'fasta', 'Invalid file extension'

        for i in range(num_inputs):
            file = input_files[i]

            # unzip if .gz file
            if file[-3:] == '.gz':
                cmd = f"gunzip {file}"
                command.execute(cmd)
                input_files[i] = file = file[:-3]

        # keep a dictionary of the distribution of read lengths in the files
        self.summary_dict = {'<50':0, '50-500': 0, '500-10000': 0, '10000+': 0}

        quick_check_passed = 
            self.quick_check_file(input_files[0], file_ext == 'fastq') and
            (num_inputs == 1 or self.quick_check_file(input_files[1], file_ext == 'fastq'))

        for infile, outfile in zip(input_files, output_files):
            if quick_check_passed:
                self.truncate_file(infile, outfile, file_ext == 'fastq', max_fragments)
            else:
                self.full_check_and_truncate_file(file, outfile, file_ext == 'fastq', max_fragments)

        with open(summary_file, 'w') as summary_f:
            json.dump(self.summary_dict, summary_f)

        return

    # quick_check_file returns:
    #   True if the first 100 fragments all have the same length of reads and
    #   are well-formated single-line FASTA / FASTQ entries. 
    #
    #   False if the entries are not formatted simply or the read lengths are not
    #   all identical or if there is another possibly recoverable error
    #
    #   Throws an exception in the case of an unrecoverable abnormality
    def quick_check_file(self, file, is_fastq, max_fragments_to_check = 100):
        num_fragments = 0
        fragment_length = 0

        with open(file, 'r', encoding='utf-8') as input_f:
            while True:
                num_fragments += 1
                if num_fragments > max_fragments_to_check:
                    break

                identifier_l = input_f.readline()
                if len(identifier_l) == 0: # EOF
                    break

                read_l = input_f.readline()
                if len(read_l) == 0: # unexpected EOF
                    raise RuntimeError(f"Invalid input file, unexpected EOF: {file}")

                if is_fastq:
                    identifier2_l = input_f.readline()
                    if len(identifier2_l) == 0:
                        raise RuntimeError(f"Invalid FASTQ file, unexpected EOF: {file}")

                    quality_l = input_f.readline()
                    if len(quality_l) == 0:
                        raise RuntimeError(f"Invalid FASTQ file, unexpected EOF: {file}")

                if is_fastq:
                    if identifier_l[0] != '@' or identifier2_l[0] != '+':
                        # may be FASTQ file with multi-line reads, requires full check
                        return False
                else:
                    if identifier_l[0] != '>':
                         # may be FASTA file with multi-line reads, requires full check
                        return False

                if fragment_length == 0:
                    fragment_length = len(identifier_l)
                    if fragment_length < 50 or fragment_length > 500:
                         # non-standard fragment lengths require more detailed examination
                        return False

                if fragment_length != len(identifier_l) or fragment_length != len(identifier2_l):
                     # file does not meet "quick check" requirements since fragments/quality
                     # scores are not all of same length
                    return False

        return True

    def truncate_file(self, infile, outfile, is_fastq, max_fragments):
        if is_fastq:
            num_lines = max_fragments * 4
        else:
            num_lines = max_fragments * 2
        command.execute(f"head -n {num_lines} {infile} > {outfile}")
        self.summary_dict['50-500'] = count.reads_in_group(self.output_files_local()[0:2])
        return

    # full_check_and_truncate_file does an exhaustive check of the input file, up to
    # max_fragments, and reformats the output to conform to what the rest of the 
    # computational pipeline expects (single-line reads of max-length 10,000). After
    # viewing max_fragments reads or encountering EOF, the function returns.
    #
    # Throws an exception in the case of an unrecoverable abnormality
    def full_check_and_truncate_file(self, infile, outfile, is_fastq, max_fragments):
        num_fragments = 0
        fragment_length = 0

        with open(infile, 'r', encoding='utf-8') as input_f, open(outfile, 'w') as output_f:
            next_line = input_f.readline()
            while True:
                num_fragments += 1
                if num_fragments > max_fragments:
                    break

                identifier_l = next_line
                if len(identifier_l) == 0: # EOF
                    break

                read_l = input_f.readline()
                if len(read_l) == 0:
                    raise RuntimeError(f"Invalid input file, unexpected EOF: {file}")

                read_l = read_l.rstrip()
                next_line = input_f.readline()
                while len(next_line) > 0 and next_line[0] not in ['>', '@', '+']:
                    read_l += next_line.rstrip()
                    next_line = input_f.readline()
                    
                if is_fastq:
                    identifier2_l = next_line
                    if len(identifier2_l) == 0:
                        raise RuntimeError(f"Invalid FASTQ file, unexpected EOF: {file}")
                    
                    quality_l = input_f.readline()
                    if len(quality_l) == 0:
                        raise RuntimeError(f"Invalid FASTQ file, unexpected EOF: {file}")    
                    
                    quality_l = quality_l.rstrip()
                    next_line = input_f.readline()
                    while len(next_line) > 0 and next_line[0] not in ['>', '@', '+']:
                        quality_l += next_line.rstrip()
                        next_line = input_f.readline()                    

                if is_fastq:
                    if identifier_l[0] != '@':
                        raise RuntimeError(f"Invalid FASTQ file: {infile}, \
                            invalid identifier: {identifier_l}")
                    if identifier2_l[0] != '+':
                        raise RuntimeError(f"Invalid FASTQ file: {infile}, \
                            invalid identifier: {identifier2_l}")
                else:
                    if identifier_l[0] != '>':
                        raise RuntimeError(f"Invalid FASTA file: {infile}, \
                            invalid identifier: {identifier_l}")

                # At this point, identifier_l and identifier2_l end in a newline and
                # read_l and quality_l do not end in a newline
                read_len = len(read_l)
                
                # Force read and quality lengths to be identical, either by padding quality
                # with the last quality score or truncating quality score
                if is_fastq:
                    if read_len > len(quality_l):
                        quality_l += (quality_l[-1] * (read_len - len(quality_l)))
                    elif read_len < len(quality_l):
                        quality_l = quality_l[0:read_len]

                if read_len < 50:
                    self.summary_dict['<50'] += 1
                    continue
                elif read_len < 500:
                    self.summary_dict['50-500'] += 1
                elif read_len < 10000:
                    self.summary_dict['500-10000'] += 1
                else:
                    self.summary_dict['10000+'] += 1
                    read_l = read_l[0:10000]
                    if is_fastq:
                        quality_l = quality_l[0:10000]

                output_f.write(identifier_l + read_l + "\n")
                if is_fastq:
                    output_f.write(identifier2_l + quality_l + "\n")

        return
    
    def count_reads(self):
        self.should_count_reads = True
        self.counts_dict[self.name] = self.summary_dict['50-500'] +
                                      self.summary_dict['500-10000'] +
                                      self.summary_dict['10000+']