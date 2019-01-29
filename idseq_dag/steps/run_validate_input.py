import json

from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.count as count
import idseq_dag.util.validate_constants as vc

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
        self.summary_dict = {vc.BUCKET_TOO_SHORT:0,
                             vc.BUCKET_NORMAL: 0,
                             vc.BUCKET_LONG: 0,
                             vc.BUCKET_TOO_LONG: 0}

        try:
            quick_check_passed = \
                self.quick_check_file(input_files[0], file_ext == 'fastq') and \
                (num_inputs == 1 or self.quick_check_file(input_files[1], file_ext == 'fastq'))

            all_fragments = []

            for infile, outfile in zip(input_files, output_files):
                if quick_check_passed:
                    num_fragments = self.truncate_file(infile, outfile, file_ext == 'fastq', max_fragments)
                else:
                    num_fragments = self.full_check_and_truncate_file(infile, outfile, file_ext == 'fastq', max_fragments)
                all_fragments.append(num_fragments)

            if len(all_fragments) == 2 and abs(all_fragments[1]-all_fragments[0]) > 1000:
                raise RuntimeError(f"Paired input files need to contain the same number of reads")

            with open(summary_file, 'w') as summary_f:
                json.dump(self.summary_dict, summary_f)
        except Exception as e:
            with open(summary_file, 'w') as summary_f:
                json.dump({'Validation error': str(e)}, summary_f)
                s3_path = self.s3_path(summary_f)
                idseq_dag.util.s3.upload_with_retries(summary_f, s3_path)
            raise e

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
                    fragment_length = len(read_l)
                    if fragment_length < vc.READ_LEN_CUTOFF_LOW or fragment_length > vc.READ_LEN_CUTOFF_MID:
                         # non-standard fragment lengths require more detailed examination
                        return False

                if fragment_length != len(read_l) or fragment_length != len(quality_l):
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
        num_fragments = count.reads(outfile)
        self.summary_dict[vc.BUCKET_NORMAL] += num_fragments
        return num_fragments

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
                    raise RuntimeError(f"Invalid input file, unexpected EOF: {infile}")

                read_l = read_l.rstrip()
                next_line = input_f.readline()
                while len(next_line) > 0 and next_line[0] not in ['>', '@', '+']:
                    read_l += next_line.rstrip()
                    next_line = input_f.readline()

                if is_fastq:
                    identifier2_l = next_line
                    if len(identifier2_l) == 0:
                        raise RuntimeError(f"Invalid FASTQ file, unexpected EOF: {infile}")

                    quality_l = input_f.readline()
                    if len(quality_l) == 0:
                        raise RuntimeError(f"Invalid FASTQ file, unexpected EOF: {infile}")

                    quality_l = quality_l.rstrip()
                    next_line = input_f.readline()
                    while len(next_line) > 0 and next_line[0] not in ['>', '@', '+']:
                        quality_l += next_line.rstrip()
                        next_line = input_f.readline()

                if is_fastq:
                    if identifier_l[0] != '@':
                        raise RuntimeError(f"Invalid FASTQ file: {infile}, " + \
                            f"invalid identifier: {identifier_l}")
                    if identifier2_l[0] != '+':
                        raise RuntimeError(f"Invalid FASTQ file: {infile}, " + \
                            f"invalid identifier: {identifier2_l}")
                else:
                    if identifier_l[0] != '>':
                        raise RuntimeError(f"Invalid FASTA file: {infile}, " + \
                            f"invalid identifier: {identifier_l}")

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

                if read_len < vc.READ_LEN_CUTOFF_LOW:
                    self.summary_dict[vc.BUCKET_TOO_SHORT] += 1
                    continue
                elif read_len < vc.READ_LEN_CUTOFF_MID:
                    self.summary_dict[vc.BUCKET_NORMAL] += 1
                elif read_len < vc.READ_LEN_CUTOFF_HIGH:
                    self.summary_dict[vc.BUCKET_LONG] += 1
                else:
                    self.summary_dict[vc.BUCKET_TOO_LONG] += 1
                    read_l = read_l[0:10000]
                    if is_fastq:
                        quality_l = quality_l[0:10000]

                output_f.write(identifier_l + read_l + "\n")
                if is_fastq:
                    output_f.write(identifier2_l + quality_l + "\n")

        return num_fragments

    def count_reads(self):
        self.should_count_reads = True
        self.counts_dict[self.name] = self.summary_dict[vc.BUCKET_NORMAL] + \
                                      self.summary_dict[vc.BUCKET_LONG] + \
                                      self.summary_dict[vc.BUCKET_TOO_LONG]
