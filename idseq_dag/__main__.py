#!/usr/bin/env python3

import argparse
import os
import idseq_dag.util.log as log
from idseq_dag.engine.pipeline_flow import PipelineFlow
from idseq_dag import __version__

log.configure_logger()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument('dag_json', help='pipeline dag in json file format')
    parser.add_argument('--no-lazy-run', dest='lazy_run', action='store_false')
    parser.add_argument('--key-path-s3', dest='key_path_s3', help='ssh key')
    parser.set_defaults(lazy_run=True)
    args = parser.parse_args()
    if args.key_path_s3:
        os.environ["KEY_PATH_S3"] = args.key_path_s3
    try:
        flow = PipelineFlow(lazy_run=args.lazy_run,
                            dag_json=args.dag_json)
        log.write("everything is awesome. idseq dag is valid~")
    except:
        parser.print_help()
        raise
    log.write("start executing the dag")
    flow.start()
    log.write("all steps are done")
    log.write("Please contact us at idseqhelp@chanzuckerberg.com or Slack"
              " #idseqhelp within CZ Biohub for assistance.", user_friendly=True)
    log.write("============= START USER FRIENDLY LOG =============")
    log.write(f"\n{log.user_friendly_stream.getvalue()}")
    log.write("============= END USER FRIENDLY LOG =============")


if __name__ == "__main__":
    main()
