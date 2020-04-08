import unittest
import tempfile
import os
import shelve
import random
import dbm

from idseq_dag.util.m8 import call_hits_m8

from tests.unit.unittest_helpers import file_contents, relative_file_path

class TestLog(unittest.TestCase):
    def test_call_hits_m8(self):
        input_m8 = relative_file_path(__file__, 'm8-test/gsnap.m8')

        lineages = relative_file_path(__file__, 'm8-test/taxid-lineages.db')
        taxid2wgs_accession = relative_file_path(__file__, 'm8-test/taxid2wgs_accession.db')

        """
        # Create an empty lineages db
        lineages_db = shelve.open(lineages.replace('.db', ''), 'c')
        lineages_db.close()

        taxid2wgs_accession_db = shelve.open(taxid2wgs_accession.replace('.db', ''), 'c')
        taxid2wgs_accession_db.close()
        """

        output_m8 = relative_file_path(__file__, 'm8-test/test.m8')
        output_summary = relative_file_path(__file__, 'm8-test/test.hitsummary.tab')

        call_hits_m8(
            input_m8,
            lineages,
            taxid2wgs_accession,
            output_m8,
            output_summary,
            36,
        )

        in_size = os.stat(input_m8).st_size
        out_size = os.stat(output_m8).st_size

        # File should shrink due to deduping
        self.assertLessEqual(out_size, in_size)

        sample_m8 = relative_file_path(__file__, 'm8-test/sample.m8')
        sample_summary = relative_file_path(__file__, 'm8-test/sample.hitsummary.tab')

        self.assertEqual(file_contents(output_m8), file_contents(sample_m8))
        self.assertEqual(file_contents(output_summary), file_contents(sample_summary))

        os.remove(output_m8)
        os.remove(output_summary)
        os.remove(lineages)
        os.remove(taxid2wgs_accession)
