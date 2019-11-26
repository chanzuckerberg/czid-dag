# Description
*The intended change, motivation and consequences.*

# Version
[ ] I have increased the appropriate version number in https://github.com/chanzuckerberg/idseq-dag/blob/master/idseq_dag/__init__.py
[ ] I have added release notes for my new version to https://github.com/chanzuckerberg/idseq-dag/blob/master/README.md

Version numbers for this repo take the form X.Y.Z.
- We increase Z for a change that does not add or change results in any way. Example: adding a log statement.
- We increase Y for a change that adds results, changes results, or has the potential to change results. Example: changing a parameter in the GSNAP command.
- We increase X for a paradigm shift in how the pipeline is conceived. Example: adding a de-novo assembly step and then reassigning hits based on the assembled contigs.
Changes to X or Y force recomputation of all results when a sample is rerun using idseq-web. Changes to Z do not force recomputation when the sample is rerun - the pipeline will lazily reuse existing S3 outputs.

# Notes

*Optional observations, comments or explanations.*

# Tests

[ ] I have verified in IDseq staging that the pipeline still completes successfully:
    [ ] for single-end inputs
    [ ] for paired-end inputs
    [ ] for FASTQ inputs
    [ ] for FASTA inputs.
[ ] I have validated that my change does not introduce any correctness bugs to existing output types.
