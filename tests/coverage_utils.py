import unittest

import idseq_dag.util.coverage as coverage_utils

class CoverageUtils(unittest.TestCase):
    # Basic test of generate_read_group_json
    def test_generate_read_group_json_basic(self):
        accession_id = "ACCESSION_1"

        accession_data = {
            "total_length": 100,
            "name": "Test Name",
            "contigs": ["CONTIG_1", "CONTIG_2"],
            "reads": ["READS_1", "READS_2", "READS_3", "READS_4"]
        }

        contigs_map = {
            "CONTIG_1": {
                "total_length": 25,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 25,
                "subject_start": 40,
                "subject_end": 64,
                "coverage": [2] * 25,
                "prop_mismatch": 0.1,
                "num_reads": 22,
                "alignment_length": 25,
                "percent_id": 100,
                "num_mismatches": 6,
                "num_gaps": 3,
                "byterange": [1, 100],
            },
             "CONTIG_2": {
                "total_length": 9,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 9,
                "subject_start": 85,
                "subject_end": 93,
                "coverage": [2, 2, 2, 2, 2, 2, 2, 2, 2],
                "num_reads": 15,
                "alignment_length": 9,
                "percent_id": 88,
                "num_mismatches": 4,
                "num_gaps": 5,
                "byterange": [100, 200],
            }
        }

        reads_map = {
            "READS_1": {
                "accession": "ACCESSION_1",
                "subject_start": 1,
                "subject_end": 30,
                "alignment_length": 30,
                "percent_id": 98,
                "num_mismatches": 1,
                "num_gaps": 0,
            },
            "READS_2": {
                "accession": "ACCESSION_1",
                "subject_start": 69,
                "subject_end": 75,
                "alignment_length": 7,
                "percent_id": 97,
                "num_mismatches": 0,
                "num_gaps": 1,
            },
            "READS_3": {
                "accession": "ACCESSION_1",
                "subject_start": 83,
                "subject_end": 89,
                "alignment_length": 7,
                "percent_id": 90,
                "num_mismatches": 3,
                "num_gaps": 1,
            },
            "READS_4": {
                "accession": "ACCESSION_1",
                "subject_start": 90,
                "subject_end": 86,
                "alignment_length": 5,
                "percent_id": 95,
                "num_mismatches": 2,
                "num_gaps": 0,
            }
        }

        read_group_json = coverage_utils.generate_hit_group_json(
          accession_data, accession_id, contigs_map, reads_map, 10
        )

        self.assertEqual(len(read_group_json), 4)

        self.assertEqual(read_group_json, [
            [0, 1, 0, 1, 30, 30, .98, 1, 0, 1, []],
            [1, 0, 22, 40, 64, 25, 1, 6, 3, 5, [[1, 100]]],
            [0, 1, 0, 69, 75, 7, .97, 0, 1, 7, []],
            # CONTIG_2, READ_3, and READ_4 are grouped together and their stats are averaged.
            [1, 2, 15, 83, 93, 7, .91, 3, 2, 8, [[100, 200]]],
        ])

 # generate_read_group_json aggregates contigs
    def test_generate_read_group_json_contigs(self):
        accession_id = "ACCESSION_1"

        accession_data = {
           "total_length": 10,
            "name": "Test Name",
            "contigs": ["CONTIG_1", "CONTIG_2"],
            "reads": [],
        }


        contigs_map = {
            "CONTIG_1": {
                "total_length": 6,
                "accession": "ACCESSION_1",
                "query_start": 3,
                "query_end": 4,
                "subject_start": 7,
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

        read_group_json = coverage_utils.generate_hit_group_json(
          accession_data, accession_id, contigs_map, [], 5
        )

        self.assertEqual(len(read_group_json), 1)

        self.assertEqual(read_group_json, [
            [2, 0, 9, 7, 8, 1.5, .985, 3, 3.5, 3, [[1, 100], [100, 200]]],
        ])

    # Basic test of remove_taxons_with_no_contigs
    def test_remove_taxons_with_no_contigs_basic(self):
        accession_data = {
            "ACCESSION_1": {
                "reads": ["READS_1"],
                "contigs": ["CONTIGS_1"]
            },
            "ACCESSION_2": {
                "reads": ["READS_2"],
                "contigs": []
            },
            "ACCESSION_3": {
                "reads": ["READS_3"],
                "contigs": []
            },
            "ACCESSION_4": {
                "reads": ["READS_4"],
                "contigs": []
            }
        }

        taxons_to_accessions = {
            "TAXON_1": ["ACCESSION_1", "ACCESSION_2"],
            "TAXON_2": ["ACCESSION_3", "ACCESSION_4"]
        }

        coverage_utils.remove_taxons_with_no_contigs(
          accession_data, taxons_to_accessions
        )

        self.assertEqual(len(accession_data.keys()), 2)
        self.assertTrue("ACCESSION_1" in accession_data)
        self.assertTrue("ACCESSION_2" in accession_data)

        self.assertEqual(len(taxons_to_accessions.keys()), 1)
        self.assertTrue("TAXON_1" in taxons_to_accessions)


    # Basic test of get_unassigned_reads_set
    def test_get_unassigned_reads_set_basic(self):
        accession_data = {
            "ACCESSION_1": {
                "reads": ["READS_1", "READS_3", "READS_4"],
                "contigs": ["CONTIGS_1"]
            },
            "ACCESSION_2": {
                "reads": ["READS_2"],
                "contigs": []
            }
        }

        unassigned_reads_set = coverage_utils.get_unassigned_reads_set(
          accession_data
        )

        self.assertEqual(len(unassigned_reads_set), 4)
        self.assertTrue("READS_1" in unassigned_reads_set)
        self.assertTrue("READS_2" in unassigned_reads_set)
        self.assertTrue("READS_3" in unassigned_reads_set)
        self.assertTrue("READS_4" in unassigned_reads_set)


    def test_generate_coverage_viz_summary_json_basic(self):
        accession_data = {
            "ACCESSION_1": {
                "reads": ["READS_1"],
                "contigs": ["CONTIGS_1"],
                "name": "Test Accession 1",
                "score": 100,
            },
            "ACCESSION_2": {
                "reads": ["READS_2", "READS_3"],
                "contigs": [],
                "name": "Test Accession 2",
                "score": 200,
            }
        }

        taxons_to_accessions = {
            "TAXON_1": ["ACCESSION_1", "ACCESSION_2"],
        }

        coverage_viz_obj = {
            "ACCESSION_1": {
                "coverage_depth": 10,
            },
            "ACCESSION_2": {
                "coverage_depth": 20,
            }
        }

        taxons_to_total_accession_count = {
            "TAXON_1": 100
        }

        taxons_to_accessions_json = coverage_utils.generate_coverage_viz_summary_json(
            taxons_to_accessions, accession_data, coverage_viz_obj, taxons_to_total_accession_count
        )

        self.assertEqual(len(taxons_to_accessions_json.keys()), 1)
        self.assertEqual(taxons_to_accessions_json["TAXON_1"], {
            "best_accessions": [
                {
                    "id": "ACCESSION_1",
                    "name": "Test Accession 1",
                    "num_contigs": 1,
                    "num_reads": 1,
                    "score": 100,
                    "coverage_depth": 10
                },
                {
                    "id": "ACCESSION_2",
                    "name": "Test Accession 2",
                    "num_contigs": 0,
                    "num_reads": 2,
                    "score": 200,
                    "coverage_depth": 20
                }
            ],
            "num_accessions": 100
        })


    # Test that select best accessions correctly filters to top accessions.
    def test_select_best_accessions_per_taxon_basic(self):
        taxons_to_accessions = {
            "TAXON_1": ["ACCESSION_1", "ACCESSION_2", "ACCESSION_3"],
            "TAXON_2": ["ACCESSION_4", "ACCESSION_5", "ACCESSION_6"],
            "TAXON_3": ["ACCESSION_7"],
        }

        contigs_map = {
            "CONTIGS_1": { "alignment_length": 100 },
            "CONTIGS_2": { "alignment_length": 100 },
            "CONTIGS_3": { "alignment_length": 100 },
            "CONTIGS_4": { "alignment_length": 120 },
            "CONTIGS_5": { "alignment_length": 120 },
            "CONTIGS_6": { "alignment_length": 300 },
            "CONTIGS_7": { "alignment_length": 200 },
        }

        reads_map = {
            "READS_1": { "alignment_length": 150 },
            "READS_2": { "alignment_length": 150 },
            "READS_3": { "alignment_length": 150 },
            "READS_4": { "alignment_length": 150 },
        }

        accession_data = {
            "ACCESSION_1": {
                "reads": [],
                "contigs": ["CONTIGS_1", "CONTIGS_2", "CONTIGS_3"],
            },
            "ACCESSION_2": {
                "reads": [],
                "contigs": ["CONTIGS_4", "CONTIGS_5"],
            },
            "ACCESSION_3": {
                "reads": [],
                "contigs": ["CONTIGS_6"],
            },
            "ACCESSION_4": {
                "reads": ["READS_1"],
                "contigs": [],
            },
            "ACCESSION_5": {
                "reads": [],
                "contigs": ["CONTIGS_7"],
            },
            "ACCESSION_6": {
                "reads": ["READS_2", "READS_3"],
                "contigs": [],
            },
            "ACCESSION_7": {
                "reads": ["READS_4"],
                "contigs": [],
            },
        }

        (new_taxons_to_accessions, new_accession_data) = coverage_utils.select_best_accessions_per_taxon(
            taxons_to_accessions, accession_data, contigs_map, reads_map, 2
        )

        self.assertEqual(new_taxons_to_accessions, {
            "TAXON_1": ["ACCESSION_3", "ACCESSION_1", "ACCESSION_2"],
            "TAXON_2": ["ACCESSION_5", "ACCESSION_6"],
            "TAXON_3": ["ACCESSION_7"],
        })

        self.assertFalse("ACCESSION_4" in new_accession_data)

        self.assertTrue(new_accession_data["ACCESSION_3"]["score"] > new_accession_data["ACCESSION_1"]["score"])
        self.assertTrue(new_accession_data["ACCESSION_1"]["score"] > new_accession_data["ACCESSION_2"]["score"])
        self.assertTrue(new_accession_data["ACCESSION_5"]["score"] > new_accession_data["ACCESSION_6"]["score"])

    # Test that select best accessions correctly filters to top accessions.
    def test_select_best_accessions_per_taxon_high_read_count(self):
        taxons_to_accessions = {
            "TAXON_1": ["ACCESSION_1", "ACCESSION_2", "ACCESSION_3"],
            "TAXON_2": ["ACCESSION_4", "ACCESSION_5"],
        }

        contigs_map = {
            "CONTIGS_1": { "alignment_length": 1 },
            "CONTIGS_2": { "alignment_length": 1 },
            "CONTIGS_3": { "alignment_length": 1 },
            "CONTIGS_4": { "alignment_length": 1 },
        }

        reads_map = {
            "READS_1": { "alignment_length": 1 },
            "READS_2": { "alignment_length": 1 },
            "READS_3": { "alignment_length": 1 },
            "READS_4": { "alignment_length": 1 },
            "READS_5": { "alignment_length": 1 },
            "READS_6": { "alignment_length": 1 },
        }

        accession_data = {
            "ACCESSION_1": {
                "reads": ["READS_1", "READS_2", "READS_3"],
                "contigs": [],
            },
            "ACCESSION_2": {
                "reads": [],
                "contigs": ["CONTIGS_1"],
            },
            "ACCESSION_3": {
                "reads": [],
                "contigs": ["CONTIGS_2"],
            },
            "ACCESSION_4": {
                "reads": ["READS_4", "READS_5", "READS_6"],
                "contigs": ["CONTIGS_4"],
            },
            "ACCESSION_5": {
                "reads": [],
                "contigs": ["CONTIGS_3"],
            },
        }

        (new_taxons_to_accessions, new_accession_data) = coverage_utils.select_best_accessions_per_taxon(
            taxons_to_accessions, accession_data, contigs_map, reads_map, 2
        )

        self.assertEqual(new_taxons_to_accessions, {
            "TAXON_1": ["ACCESSION_2", "ACCESSION_3"],
            "TAXON_2": ["ACCESSION_4", "ACCESSION_5"],
        })

        self.assertFalse("ACCESSION_1" in new_accession_data)

        self.assertTrue(new_accession_data["ACCESSION_2"]["score"] == new_accession_data["ACCESSION_3"]["score"])
        self.assertTrue(new_accession_data["ACCESSION_4"]["score"] > new_accession_data["ACCESSION_5"]["score"])
