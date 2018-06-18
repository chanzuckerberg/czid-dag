import idseq_dag.util.command as command

def reads(local_file_path):
    ''' Count reads in a local file '''
    line_count = int(command.execute_with_output("wc -l < {}".format(local_file_path)))
    file_format = local_file_path.split(".")[-1]
    assert file_format in ["fq", "fastq"], "Currently only support fq/fastq"
    return lines2reads(line_count, file_format)

def lines2reads(line_count, file_format):
    ''' Convert line count to read count based on file format '''
    read_count = None
    if file_format in ["fq", "fastq"]:
        read_count = 0.25 * line_count
    return read_count
