# Validates line length and checks if all characters are part of ascii set
# Required input variables:
#    max_line_length: This is the max allowed length for each line
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
