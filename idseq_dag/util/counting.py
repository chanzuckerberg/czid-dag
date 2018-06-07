def lines2reads(line_count, file_format):
    ''' Convert line count to read count based on file format '''
    read_count = None
    if file_format in ["fq", "fastq"]:
        read_count = 0.25 * line_count
    return read_count

def reads2lines(read_count, file_format):
    ''' Convert read count to line count based on file format '''
    line_count = None
    if file_format in ["fq", "fastq"]:
        line_count = 4 * read_count
    return line_count
