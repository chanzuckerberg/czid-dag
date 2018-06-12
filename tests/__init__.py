import unittest

from .run_bowtie2 import RunBowtie2Test
from .run_cdhitdup import RunCDHitDupTest
from .run_lzw import RunLZWTest
import idseq_dag.util.log as log

log.configure_logger()

if __name__ == '__main__':
  unittest.main(verbosity=2)


