import unittest

import idseq_dag.util.coverage as coverage_utils

class GenerateCoverageVizJson(unittest.TestCase):
    # Basic test of generate_coverage_viz_json
    def test_calculate_accession_stats_basic(self):
        accession_data = {
            "ACCESSION_1": {
                "total_length": 10,
                "name": "Test Name",
                "contigs": ["CONTIG_1", "CONTIG_2"],
                "reads": ["READS_1", "READS_2"]
            }
        }

        contigs_map = {
            "CONTIG_1": {
                "total_length": 6,
                "accession": "ACCESSION_1",
                "query_start": 3,
                "query_end": 4,
                "subject_start": 6,
                "subject_end": 7,
                "coverage": [1, 2, 3, 4, 5, 6],
                "prop_mismatch": 0.1,
                "num_reads": 4,
                "alignment_length": 2,
                "percent_id": 99,
                "num_mismatches": 2,
                "num_gaps": 1,
                "byterange": [1, 100],
            },
             "CONTIG_2": {
                "total_length": 6,
                "accession": "ACCESSION_1",
                "query_start": 4,
                "query_end": 4,
                "subject_start": 8,
                "subject_end": 8,
                "coverage": [2, 2, 2, 2, 2, 2],
                "prop_mismatch": 0.1,
                "num_reads": 5,
                "alignment_length": 1,
                "percent_id": 98,
                "num_mismatches": 4,
                "num_gaps": 6,
                "byterange": [100, 200],
            }
        }

        reads_map = {
            "READS_1": {
                "accession": "ACCESSION_1",
                "subject_start": 2,
                "subject_end": 3,
                "prop_mismatch": 0.3,
                "alignment_length": 2,
                "percent_id": 98,
                "num_mismatches": 1,
                "num_gaps": 0,
            },
            "READS_2": {
                "accession": "ACCESSION_1",
                "subject_start": 10,
                "subject_end": 10,
                "prop_mismatch": 0.3,
                "alignment_length": 1,
                "percent_id": 97,
                "num_mismatches": 0,
                "num_gaps": 1,
            }
        }

        coverage_viz_json = coverage_utils.generate_coverage_viz_json(
          accession_data, contigs_map, reads_map, 5
        )

        self.assertTrue("ACCESSION_1" in coverage_viz_json)

        accession_obj = coverage_viz_json["ACCESSION_1"]

        self.assertEqual(accession_obj["total_length"], 10)
        self.assertEqual(accession_obj["name"], "Test Name")

        # generate_read_group_json is tested further elsewhere.
        self.assertEqual(accession_obj["hit_groups"], [
            [0, 1, 0, 2, 3, 2, .98, 1, 0, 1, []],
            [1, 0, 4, 6, 7, 2, .99, 2, 1, 3, [[1, 100]]],
            [1, 0, 5, 8, 8, 1, .98, 4, 6, 3, [[100, 200]]],
            [0, 1, 0, 10, 10, 1, .97, 0, 1, 4, []]
        ])

        # calculate_accession_coverage is tested further elsewhere.
        self.assertEqual(accession_obj["coverage"], [
            [0, 0.5, 0.5, 0, 1],
            [1, 0.5, 0.5, 0, 1],
            [2, 1.5, 0.5, 1, 0],
            [3, 3.0, 1, 2, 0],
            [4, 0.5, 0.5, 0, 1],
        ])
        self.assertEqual(accession_obj["coverage_bin_size"], 2)

        self.assertEqual(accession_obj["max_aligned_length"], 2)
        self.assertEqual(accession_obj["coverage_depth"], 1.2)
        self.assertEqual(accession_obj["coverage_breadth"], 0.6)
        self.assertEqual(accession_obj["avg_prop_mismatch"], 0.2)
