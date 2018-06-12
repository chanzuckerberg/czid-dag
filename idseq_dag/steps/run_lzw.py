from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command

def lzw_fraction(sequence):
    if sequence == "":
        return 0.0
    sequence = sequence.upper()
    dict_size = 0
    dictionary = {}
    # Initialize dictionary with single char
    for c in sequence:
        dict_size += 1
        dictionary[c] = dict_size

    word = ""
    results = []
    for c in sequence:
        wc = word + c
        if dictionary.get(wc):
            word = wc
        else:
            results.append(dictionary[word])
            dict_size += 1
            dictionary[wc] = dict_size
            word = c
    if word != "":
        results.append(dictionary[word])
    return float(len(results)) / len(sequence)



def generate_lzw_filtered_single(fasta_file, out_file, cutoff_fractions):
    out_read_1 = open(out_prefix + '.1.fasta', 'wb')
    for cutoff_frac in cutoff_fractions:
        read_1 = open(fasta_file, 'rb')
        count = 0
        filtered = 0
        while True:
            line_r1_header = read_1.readline()
            line_r1_sequence = read_1.readline()
            if line_r1_header and line_r1_sequence:
                fraction_1 = lzw_fraction(line_r1_sequence.rstrip())
                count += 1
                if fraction_1 > cutoff_frac:
                    out_read_1.write(line_r1_header)
                    out_read_1.write(line_r1_sequence)
                else:
                    filtered += 1
            else:
                break
        msg = "LZW filter: cutoff_frac: %f, total reads: %d, filtered reads: %d, " \
              "kept ratio: %f" % (cutoff_frac, count, filtered, 1 - float(filtered) / count)
        write_to_log(msg)
        if count != filtered:
            break
    out_read_1.close()

def generate_lzw_filtered_paired(fasta_file_1, fasta_file_2, out_prefix,
                                 cutoff_fractions):
    out_read_1 = open(out_prefix + '.1.fasta', 'wb')
    out_read_2 = open(out_prefix + '.2.fasta', 'wb')
    for cutoff_frac in cutoff_fractions:
        read_1 = open(fasta_file_1, 'rb')
        read_2 = open(fasta_file_2, 'rb')
        count = 0
        filtered = 0
        while True:
            line_r1_header = read_1.readline()
            line_r1_seq = read_1.readline()
            line_r2_header = read_2.readline()
            line_r2_seq = read_2.readline()
            if line_r1_header and line_r1_seq and line_r2_header and line_r2_seq:
                fraction_1 = lzw_fraction(line_r1_seq.rstrip())
                fraction_2 = lzw_fraction(line_r2_seq.rstrip())
                count += 1
                if fraction_1 > cutoff_frac and fraction_2 > cutoff_frac:
                    out_read_1.write(line_r1_header)
                    out_read_1.write(line_r1_seq)
                    out_read_2.write(line_r2_header)
                    out_read_2.write(line_r2_seq)
                else:
                    filtered += 1
            else:
                break
        msg = "LZW filter: cutoff_frac: %f, total reads: %d, filtered reads: %d, " \
              "kept ratio: %f" % (cutoff_frac, count, filtered, 1 - float(filtered) / count)
        write_to_log(msg)
        if count != filtered:
            break
    out_read_1.close()
    out_read_2.close()

class PipelineStepRunLZW(PipelineStep):

    def run(self):
        input_fas = self.input_files[0]
        output_fas = self.output_files_local()
        for f in self.output_files_local():
            command.execute("date > %s" % f)
