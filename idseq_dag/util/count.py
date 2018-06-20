import idseq_dag.util.command as command

def reads_in_group(file_group, max_fragments=None):
    '''
    Count reads in a group of matching files, up to a maximum number of fragments.
    The term "fragment" refers to the physical DNA fragments the reads derive from.
    If the input is single, then 1 fragment == 1 read.
    If the input is paired, then 1 fragment == 2 reads.
    '''
    num_files = len(file_group)
    if max_fragments:
        max_reads = num_files * max_fragments
    else:
        max_reads = None
    reads_in_first_file = reads(file_group[0], max_reads)
    return num_files * reads_in_first_file

def reads(local_file_path, max_reads=None):
    '''
    Count reads in a local file based on file format inferred from extension,
    up to a maximum of max_reads.
    '''
    if max_reads:
        max_lines = reads2lines(max_reads)
        assert max_lines is not None, "Could not convert max_reads to max_lines"
        truncate_pipe = "| head -n {}".format(max_lines)
    else:
        truncate_pipe = ""
    if local_file_path.endswith(".gz"):
        count_cmd = "zcat {} {} | wc -l".format(local_file_path, truncate_pipe)
        file_format = local_file_path.split(".")[-2]
    else:
        count_cmd = "cat {} {} | wc -l".format(local_file_path, truncate_pipe)
        file_format = local_file_path.split(".")[-1]
    line_count = int(command.execute_with_output(count_cmd))
    return lines2reads(line_count, file_format)

def lines2reads(line_count, file_format):
    '''
    Convert line count to read count based on file format.
    Supports fastq and fasta formats.
    TODO: add support for m8 files here once the relevant steps are added in the pipeline engine.
    '''
    read_count = None
    if file_format in ["fq", "fastq"]:
        # Each read consists of exactly 4 lines, by definition.
        assert line_count % 4 == 0, "File does not follow fastq format"
        read_count = line_count // 4
    elif file_format in ["fa", "fasta"]:
        # Each read consists of exactly 2 lines (identifier and sequence).
        # ASSUMES the format is SINGLE-LINE fasta files, where each read's sequence
        # is on a single line rather than being spread over multiple lines.
        # This is fine for now, because we only count fasta files we generated
        # ourselves, following the single-line format.
        assert line_count % 2 == 0, "File does not follow single-line fasta format"
        read_count = line_count // 2
    else: 
        read_count = line_count
    return read_count

def reads2lines(read_count, file_format):
    '''
    Convert read count to line count based on file format.
    Currently supports fastq only.
    '''
    line_count = None
    if file_format in ["fq", "fastq"]:
        line_count = 4 * read_count
    return line_count
