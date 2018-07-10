import idseq_dag.util.command as command
import idseq_dag.util.log as log


def should_match_exactly(expected, actual):
    expected_content = command.execute_with_output(f"aws s3 cp {expected} -")
    actual_content = command.execute_with_output(f"cat {actual}")
    if expected_content == actual_content:
        log.write(f"File {expected} is the same as {actual}")
    else:
        raise RuntimeError(f"{actual} not equal to {expected}")
