#!/bin/bash
DEFAULT_MAX_LINE_LENGTH=10000
while [ $# -gt 0 ]; do
  case "$1" in
    --max_line_length)
      shift
      max_line_length=$1
      ;;
    *)
      echo "*** Invalid argument $1" >>  "/dev/stderr"
      echo "Usage:" >>  "/dev/stderr"
      echo "   cat somefile.txt | ./fastq_fasta_line_validation.sh --max_line_length 20" >>  "/dev/stderr"
      exit 1
      ;;
  esac
  shift
done
max_line_length=${max_line_length:-$DEFAULT_MAX_LINE_LENGTH}
read -d "" awk_script << 'AWK_SCRIPT_END'
{
  if ($0 ~ /[^\\x20-\\x7F\\x01\\x09]/) {
    print "PARSE ERROR: stdin is not an ascii file. Line", NR, "contains non-ascii characters." > "/dev/stderr";
    exit 1;
  }
  if (length > max_line_length) {
    print "PARSE ERROR: stdin exceeds max line size of", max_line_length, ". Line", NR, "Current length:", length > "/dev/stderr";
    exit 1;
  }
  print $0;
}
AWK_SCRIPT_END

awk -v max_line_length="$max_line_length" "$awk_script"
