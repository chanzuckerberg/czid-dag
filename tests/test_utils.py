import os
import re

import idseq_dag.util.command as command
import idseq_dag.util.log as log
from tests.idseq_step_setup import IdseqStepSetup


def should_match_exactly(expected, actual):
    contents = [None, None]
    for i, path in enumerate([expected, actual]):
        if path.startswith("s3://"):
            contents[i] = s3_file_contents(path)
        else:
            contents[i] = cat_file_contents(path)

    if contents[0] == contents[1]:
        log.write(f"File {expected} is the same as {actual}")
    else:
        raise ValueError(f"{expected} not equal to {actual}")


def should_match_sam(expected, actual):
    contents = [None, None]
    for i, path in enumerate([expected, actual]):
        if path.startswith("s3://"):
            contents[i] = s3_file_contents(path)
        else:
            contents[i] = cat_file_contents(path)
        # Remove lines starting with @PG from SAM file comparison because these
        # will differ each time. Ex: @PG ID:bowtie2 PN:bowtie2 VN:2.3.2
        # CL:"/usr/local/bin...
        contents[i] = re.sub("@PG.*\n", "", contents[i])

    if contents[0] == contents[1]:
        log.write(f"File {expected} is semantically the same as {actual}")
    else:
        raise ValueError(f"{expected} does not match {actual}")


def s3_file_contents(path):
    return command.execute_with_output(f"aws s3 cp {path} -")


def cat_file_contents(path):
    return command.execute_with_output("cat " + path)


def run_step_and_match_outputs(step_class,
                               step_name,
                               dag_file,
                               test_bundle,
                               output_dir_s3,
                               with_counts=False, subfolder=None):
    step = IdseqStepSetup.get_test_step_object(step_class, step_name, dag_file,
                                               test_bundle, output_dir_s3)
    step.start()
    step.wait_until_finished()

    if with_counts:
        step.output_files.append(f"{step_name}.count")

    expected_files = [os.path.join(test_bundle, f) for f in step.output_files]
    actual_files = step.output_files_local()

    if step.additional_files_to_upload:
        actual_files += step.additional_files_to_upload

        for f in step.additional_files_to_upload:
            basename = os.path.basename(f)
            if subfolder and subfolder in f:
                to_append = os.path.join(test_bundle, subfolder, basename)
            else:
                to_append = os.path.join(test_bundle, basename)
            expected_files.append(to_append)

    for expected, actual in zip(expected_files, actual_files):
        if expected.endswith(".sam"):
            should_match_sam(expected, actual)
        else:
            should_match_exactly(expected, actual)
