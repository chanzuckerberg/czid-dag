import idseq_dag.util.command as command

def reads_in_group(file_group):
    ''' Count reads in a group of matching files '''
    num_files = len(file_group)
    reads_in_first_file = count.reads(file_group[0])
    return num_files * reads_in_first_file

def reads(local_file_path):
    ''' Count reads in a local file '''
    line_count = int(command.execute_with_output("wc -l < {}".format(local_file_path)))
    file_format = local_file_path.split(".")[-1]
    return lines2reads(line_count, file_format)

def lines2reads(line_count, file_format):
    ''' Convert line count to read count based on file format '''
    read_count = None
    if file_format in ["fq", "fastq"]:
        assert line_count % 4 == 0, "File does not follow fastq format"
        read_count = line_count // 4
    elif file_format in ["fa", "fasta"]:
        # Support only single-line fasta files, where each record's sequence
        # is on a single line rather than being spread over multiple lines.
        # This is fine for now, because we only count fasta files we generated
        # ourselves, following the single-line format.
        assert line_count % 2 == 0, "File does not follow single-line fasta format"
        read_count = line_count // 2
    else: 
        read_count = line_count
    return read_count
