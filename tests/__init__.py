import unittest

from .run_bowtie2 import RunBowtie2Test
from .run_cdhitdup import RunCDHitDupTest
from .run_lzw import RunLZWTest
from .run_priceseq import RunPriceSeqTest
from .run_star import RunStarTest
import idseq_dag.util.log as log
log.configure_logger()

if __name__ == '__main__':
    print("Configuring the log")
    log.set_up_stdout()
    log.configure_logger()
    unittest.main(verbosity=2)

# print("Configuring the log")
# log.set_up_stdout()
# log.configure_logger()
