import unittest

from .run_bowtie2 import RunBowtie2Test
from .run_cdhitdup import RunCDHitDupTest
from .run_lzw import RunLZWTest
from .run_priceseq import RunPriceSeqTest
import idseq_dag.util.log as log

if __name__ == '__main__':
    log.configure_logger()
    unittest.main(verbosity=2)
