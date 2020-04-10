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
        accession2taxid = relative_file_path(__file__, '/m8-test/accession2taxid.db')


        accession MK468611 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MK468612 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MK468613 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MK468615 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MK468617 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MH124576 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MH124577 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MH124578 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MH124579 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MH124580 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MK286896 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MK370031 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MK370032 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MK370033 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MK468608 37124
        lineage 37124 ('37124', '11019', '11018')
        accession CP015500 573
        lineage 573 ('573', '570', '543')
        accession CP015822 573
        lineage 573 ('573', '570', '543')
        accession CP015990 573
        lineage 573 ('573', '570', '543')
        accession CP016813 573
        lineage 573 ('573', '570', '543')
        accession CP016814 573
        lineage 573 ('573', '570', '543')
        accession CP018140 573
        lineage 573 ('573', '570', '543')
        accession CP018337 573
        lineage 573 ('573', '570', '543')
        accession CP018352 573
        lineage 573 ('573', '570', '543')
        accession CP018356 573
        lineage 573 ('573', '570', '543')
        accession CP018364 573
        lineage 573 ('573', '570', '543')
        accession MK468618 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MK468619 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MK468620 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MK468621 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MK468622 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MF740874 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MF773566 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MF774614 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MF774615 37124
        lineage 37124 ('37124', '11019', '11018')
        accession MF774616 37124
        lineage 37124 ('37124', '11019', '11018')
        accession CP010295 1280
        lineage 1280 ('1280', '1279', '90964')
        accession CP010296 1280
        lineage 1280 ('1280', '1279', '90964')
        accession CP010297 1280
        lineage 1280 ('1280', '1279', '90964')
        accession CP010298 1280
        lineage 1280 ('1280', '1279', '90964')
        accession CP010299 1280
        lineage 1280 ('1280', '1279', '90964')
        accession NC_038358 2065052
        lineage 2065052 ('2065052', '687333', '687329')
        accession CP017682 1280
        lineage 1280 ('1280', '1279', '90964')
        accession CP017804 1280
        lineage 1280 ('1280', '1279', '90964')
        accession AF325855 1280
        lineage 1280 ('1280', '1279', '90964')
        accession AM990992 523796
        lineage 523796 ('1280', '1279', '90964')
        accession AP014652 46170
        lineage 46170 ('1280', '1279', '90964')


        """
        # Create an empty lineages db
        lineages_db = shelve.open(lineages.replace('.db', ''), 'c')
        lineages_db
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
