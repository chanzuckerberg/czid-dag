import idseq_dag.util.command as command
import idseq_dag.util.log as log


def should_match_exactly(expected, actual):
    if expected.startswith("s3://"):
        expected_content = s3_file_contents(expected)
    else:
        expected_content = cat_file_contents(expected)

    if actual.startswith("s3://"):
        actual_content = s3_file_contents(actual)
    else:
        actual_content = cat_file_contents(actual)

    if expected_content == actual_content:
        log.write(f"File {expected} is the same as {actual}")
    else:
        raise RuntimeError(f"{expected} not equal to {actual}")


def s3_file_contents(path):
    return command.execute_with_output(f"aws s3 cp {path} -")


def cat_file_contents(path):
    return command.execute_with_output("cat " + path)
