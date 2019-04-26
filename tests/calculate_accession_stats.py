import unittest

import idseq_dag.util.coverage as coverage_utils

class CalculateAccessionStats(unittest.TestCase):
    def test_calculate_accession_stats_basic(self):
        accession_data = {
            "total_length": 10,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contig_data = {
            "CONTIG_1": {
                "total_length": 10,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 10,
                "subject_start": 1,
                "subject_end": 10,
                "coverage": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                "prop_mismatch": 0.1
            }
        }

        read_data = {}

        stats = coverage_utils.calculate_accession_stats(
          accession_data, contig_data, {}
        )

        self.assertEqual(stats["max_aligned_length"], 10)
        self.assertEqual(stats["coverage_depth"], 5.5)
        self.assertEqual(stats["coverage_breadth"], 1)
        self.assertEqual(stats["avg_prop_mismatch"], 0.1)


    # Add multiple contigs and reads
    def test_calculate_accession_stats_multiple(self):
        accession_data = {
            "total_length": 100,
            "contigs": ["CONTIG_1", "CONTIG_2"],
            "reads": ["READ_1", "READ_2"]
        }

        contig_data = {
            "CONTIG_1": {
                "total_length": 100,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 40,
                "subject_start": 1,
                "subject_end": 40,
                "coverage": [2] * 40,
                "prop_mismatch": 0.01
            },
            "CONTIG_2": {
                "total_length": 100,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 20,
                "subject_start": 100,
                "subject_end": 81,
                "coverage": [4] * 20,
                "prop_mismatch": 0.02
            }
        }

        read_data = {
            "READ_1": {
                "total_length": 20,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 20,
                "subject_start": 1,
                "subject_end": 20,
                "prop_mismatch": 0.03
            },
            "READ_2": {
                "total_length": 20,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 20,
                "subject_start": 100,
                "subject_end": 81,
                "prop_mismatch": 0.04
            }
        }

        stats = coverage_utils.calculate_accession_stats(
          accession_data, contig_data, read_data
        )

        self.assertEqual(stats["max_aligned_length"], 40)
        self.assertEqual(stats["coverage_depth"], 2)
        self.assertEqual(stats["coverage_breadth"], 0.6)
        self.assertEqual(stats["avg_prop_mismatch"], 0.025)
