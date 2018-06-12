import unittest

from .basic import MyFirstTest
from .advanced import MySecondTest
from .run_bowtie2 import RunBowtie2Test
from .run_star import RunStarTest
from .run_priceseq import RunPriceSeqTest
import idseq_dag.util.log as log

if __name__ == '__main__':
    unittest.main(verbosity=2)

log.configure_logger()
