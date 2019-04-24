import unittest

import idseq_dag.util.coverage as coverage_utils

class CalculateCoverageViz(unittest.TestCase):
    # Basic test where everything aligns nicely.
    def test_calculate_accession_coverage_basic(self):
        accession_data = {
            "total_length": 4,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contigs_map = {
            "CONTIG_1": {
                "total_length": 4,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 4,
                "subject_start": 1,
                "subject_end": 4,
                "coverage": [2] * 4,
            }
        }

        reads_map = {}

        (coverage, bin_size) = coverage_utils.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contigs_map, reads_map, 4
        )

        self.assertEqual(coverage, [
            [0, 2.0, 1, 1, 0],
            [1, 2.0, 1, 1, 0],
            [2, 2.0, 1, 1, 0],
            [3, 2.0, 1, 1, 0],
        ])
        self.assertEqual(bin_size, 1.0)

    # Test where contig is slightly smaller than accession.
    def test_calculate_accession_coverage_contig_small(self):
        accession_data = {
            "total_length": 10,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contigs_map = {
            "CONTIG_1": {
                "total_length": 9,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 9,
                "subject_start": 1,
                "subject_end": 10,
                "coverage": [2] * 5,
            }
        }

        reads_map = {}

        (coverage, bin_size) = coverage_utils.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contigs_map, reads_map, 10
        )

        self.assertEqual(coverage, [
            [0, 2.0, 1, 1, 0],
            [1, 2.0, 1, 1, 0],
            [2, 2.0, 1, 1, 0],
            [3, 2.0, 1, 1, 0],
            [4, 2.0, 1, 1, 0],
            [5, 2.0, 1, 1, 0],
            [6, 2.0, 1, 1, 0],
            [7, 2.0, 1, 1, 0],
            [8, 2.0, 1, 1, 0],
            [9, 2.0, 1, 1, 0],
        ])
        self.assertEqual(bin_size, 1.0)


    # Test where accession is slightly smaller than contig.
    def test_calculate_accession_coverage_accession_small(self):
        accession_data = {
            "total_length": 9,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contigs_map = {
            "CONTIG_1": {
                "total_length": 10,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 10,
                "subject_start": 1,
                "subject_end": 9,
                "coverage": [2] * 10,
            }
        }

        reads_map = {}

        (coverage, bin_size) = coverage_utils.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contigs_map, reads_map, 9
        )

        self.assertEqual(coverage, [
            [0, 2.0, 1, 1, 0],
            [1, 2.0, 1, 1, 0],
            [2, 2.0, 1, 1, 0],
            [3, 2.0, 1, 1, 0],
            [4, 2.0, 1, 1, 0],
            [5, 2.0, 1, 1, 0],
            [6, 2.0, 1, 1, 0],
            [7, 2.0, 1, 1, 0],
            [8, 2.0, 1, 1, 0],
        ])
        self.assertEqual(bin_size, 1.0)

    # Test where contig covers only part of accession, and coverage array is slightly small.
    def test_calculate_accession_coverage_partial_contig(self):
        accession_data = {
            "total_length": 10,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contigs_map = {
            "CONTIG_1": {
                "total_length": 7,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 7,
                "subject_start": 3,
                "subject_end": 9,
                "coverage": [2] * 6, # Coverage is 6 instead of 7.
            }
        }

        reads_map = {}

        (coverage, bin_size) = coverage_utils.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contigs_map, reads_map, 10
        )

        self.assertEqual(coverage, [
            [2, 2.0, 1, 1, 0],
            [3, 2.0, 1, 1, 0],
            [4, 2.0, 1, 1, 0],
            [5, 2.0, 1, 1, 0],
            [6, 2.0, 1, 1, 0],
            [7, 2.0, 1, 1, 0],
            [8, 2.0, 1, 1, 0],
        ])
        self.assertEqual(bin_size, 1.0)


    # Test where contig covers only part of accession and contig has varying coverage.
    def test_calculate_accession_coverage_offset_varying(self):
        accession_data = {
            "total_length": 10,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contigs_map = {
            "CONTIG_1": {
                "total_length": 5,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 5,
                "subject_start": 3,
                "subject_end": 9,
                "coverage": [1, 2, 3, 4, 5]
            }
        }

        reads_map = {}

        (coverage, bin_size) = coverage_utils.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contigs_map, reads_map, 10
        )

        self.assertEqual(coverage, [
            [2, 1.0, 1, 1, 0],
            [3, 1.5, 1, 1, 0],
            [4, 2.5, 1, 1, 0],
            [5, 3.0, 1, 1, 0],
            [6, 3.5, 1, 1, 0],
            [7, 4.5, 1, 1, 0],
            [8, 5.0, 1, 1, 0],
        ])
        self.assertEqual(bin_size, 1.0)

    # Previous test, but the alignment is reversed.
    def test_calculate_accession_coverage_offset_varying_reverse(self):
        accession_data = {
            "total_length": 10,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contigs_map = {
            "CONTIG_1": {
                "total_length": 5,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 5,
                "subject_start": 9,
                "subject_end": 3,
                "coverage": [1, 2, 3, 4, 5]
            }
        }

        reads_map = {}

        (coverage, bin_size) = coverage_utils.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contigs_map, reads_map, 10
        )

        self.assertEqual(coverage, [
            [2, 5.0, 1, 1, 0],
            [3, 4.5, 1, 1, 0],
            [4, 3.5, 1, 1, 0],
            [5, 3.0, 1, 1, 0],
            [6, 2.5, 1, 1, 0],
            [7, 1.5, 1, 1, 0],
            [8, 1.0, 1, 1, 0],
        ])
        self.assertEqual(bin_size, 1.0)

    # Test where accession size is slightly larger than max_bins.
    def test_calculate_accession_coverage_max_bins(self):
        accession_data = {
            "total_length": 10,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contigs_map = {
            "CONTIG_1": {
                "total_length": 5,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 5,
                "subject_start": 1,
                "subject_end": 10,
                "coverage": [2.0] * 5
            }
        }

        reads_map = {}

        (coverage, bin_size) = coverage_utils.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contigs_map, reads_map, 9
        )

        self.assertEqual(coverage, [
            [0, 2.0, 1, 1, 0],
            [1, 2.0, 1, 1, 0],
            [2, 2.0, 1, 1, 0],
            [3, 2.0, 1, 1, 0],
            [4, 2.0, 1, 1, 0],
            [5, 2.0, 1, 1, 0],
            [6, 2.0, 1, 1, 0],
            [7, 2.0, 1, 1, 0],
            [8, 2.0, 1, 1, 0],
        ])

        self.assertTrue((bin_size - 10 / 9) < 0.01)

    # Test where accession size is slightly larger than max_bins and contig has varying coverage.
    def test_calculate_accession_coverage_max_bins_varying(self):
        accession_data = {
            "total_length": 10,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contigs_map = {
            "CONTIG_1": {
                "total_length": 5,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 5,
                "subject_start": 1,
                "subject_end": 10,
                "coverage": [1, 2, 3, 4, 5]
            }
        }

        reads_map = {}

        (coverage, bin_size) = coverage_utils.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contigs_map, reads_map, 9
        )

        self.assertEqual(coverage, [
            [0, 1.0, 1, 1, 0],
            [1, 1.5, 1, 1, 0],
            [2, 2.0, 1, 1, 0],
            [3, 2.5, 1, 1, 0],
            [4, 3.0, 1, 1, 0],
            [5, 3.5, 1, 1, 0],
            [6, 4.0, 1, 1, 0],
            [7, 4.5, 1, 1, 0],
            [8, 5.0, 1, 1, 0],
        ])
        self.assertTrue((bin_size - 10 / 9) < 0.01)

    # Test where max bins is much smaller than accession.
    # This allows us to verify that the coverage is being reduced when the contig doesn't cover the whole bin.
    def test_calculate_accession_coverage_max_bins_small(self):
        accession_data = {
            "total_length": 10,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contigs_map = {
            "CONTIG_1": {
                "total_length": 6,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 4,
                "subject_start": 4,
                "subject_end": 7,
                "coverage": [1, 2, 3, 4, 5, 6]
            }
        }

        reads_map = {}

        (coverage, bin_size) = coverage_utils.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contigs_map, reads_map, 5
        )

        # The coverage is reduced (e.g. divided by 2) for bins with 0.5 and 2.0.
        # Since the contig only covers half of that bin.
        self.assertEqual(coverage, [
            [1, 0.5, 0.5, 1, 0],
            [2, 2.5, 1, 1, 0],
            [3, 2.0, 0.5, 1, 0],
        ])
        self.assertEqual(bin_size, 2)

    # Test multiple reads.
    def test_calculate_accession_coverage_multiple_reads(self):
        accession_data = {
            "total_length": 10,
            "contigs": [],
            "reads": ["READS_1", "READS_2"]
        }

        contigs_map = {}

        reads_map = {
            "READS_1": {
                "accession": "ACCESSION_1",
                "subject_start": 4,
                "subject_end": 7,
            },
            "READS_2": {
                "accession": "ACCESSION_1",
                "subject_start": 6,
                "subject_end": 9,
            }
        }

        (coverage, bin_size) = coverage_utils.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contigs_map, reads_map, 10
        )

        self.assertEqual(coverage, [
            [3, 1.0, 1, 0, 1],
            [4, 1.0, 1, 0, 1],
            [5, 2.0, 1, 0, 2],
            [6, 2.0, 1, 0, 2],
            [7, 1.0, 1, 0, 1],
            [8, 1.0, 1, 0, 1],
        ])
        self.assertEqual(bin_size, 1)

    # Test where max bins is much smaller than contig.
    # This allows us to verify that the coverage is being reduced when the reads and contigs don't cover the whole bin.
    def test_calculate_accession_coverage_max_bins_small_reads(self):
        accession_data = {
            "total_length": 10,
            "contigs": ["CONTIG_1"],
            "reads": ["READS_1", "READS_2"]
        }

        contigs_map = {
            "CONTIG_1": {
                "total_length": 6,
                "accession": "ACCESSION_1",
                "query_start": 3,
                "query_end": 4,
                "subject_start": 6,
                "subject_end": 7,
                "coverage": [1, 2, 3, 4, 5, 6]
            }
        }

        reads_map = {
            "READS_1": {
                "accession": "ACCESSION_1",
                "subject_start": 2,
                "subject_end": 3,
            },
            "READS_2": {
                "accession": "ACCESSION_1",
                "subject_start": 10,
                "subject_end": 10,
            }
        }

        (coverage, bin_size) = coverage_utils.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contigs_map, reads_map, 5
        )

        # [0.5, 0.5] from READS_1.
        # [1.5, 2.0] from CONTIG_1
        # last [0.5] from READS_2.
        self.assertEqual(coverage, [
            [0, 0.5, 0.5, 0, 1],
            [1, 0.5, 0.5, 0, 1],
            [2, 1.5, 0.5, 1, 0],
            [3, 2.0, 0.5, 1, 0],
            [4, 0.5, 0.5, 0, 1],
        ])
        self.assertEqual(bin_size, 2)
