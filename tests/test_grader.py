"""Tests for grade_classification.

Run with:  python -m unittest discover -s tests -v
"""

import unittest

from classification_qc import CRITICAL, WARNING, grade_classification


def read(read_id, tax, seq="ACGT"):
    return {"read_id": read_id, "seq": seq, "tax_name": tax}


def types_of(flags):
    return {f.flag_type for f in flags}


class TestEmptyAndDegenerateInput(unittest.TestCase):
    def test_empty_read_list_returns_no_flags(self):
        self.assertEqual(grade_classification("gut", []), [])

    def test_empty_read_list_does_not_divide_by_zero(self):
        # Guards the unclassified-ratio and dominance calculations.
        try:
            grade_classification("gut", [])
        except ZeroDivisionError:
            self.fail("empty input should return cleanly, not divide by zero")

    def test_all_reads_unclassified_does_not_raise(self):
        # classified_total is zero here, so the dominance check must be
        # skipped rather than calling most_common() on an empty Counter.
        reads = [read(f"r{i}", "NA") for i in range(10)]
        flags = grade_classification("gut", reads)
        self.assertIn("HIGH_UNCLASSIFIED_RATE", types_of(flags))
        self.assertNotIn("DOMINANT_TAXON", types_of(flags))


class TestUnclassifiedRepresentations(unittest.TestCase):
    def test_na_string_counts_as_unclassified(self):
        reads = [read("r1", "NA")] + [read(f"r{i}", "Bacteroides fragilis") for i in range(2, 4)]
        flags = grade_classification("gut", reads)
        self.assertNotIn("UNEXPECTED_NICHE_TAXON", types_of(flags))

    def test_none_counts_as_unclassified(self):
        # A pipeline may emit either the literal string or a null.
        reads = [read(f"r{i}", None) for i in range(6)] + [read("r7", "Bacteroides fragilis")]
        flags = grade_classification("gut", reads)
        self.assertIn("HIGH_UNCLASSIFIED_RATE", types_of(flags))


class TestHostDnaLeak(unittest.TestCase):
    def test_human_read_is_flagged_critical(self):
        reads = [read("r1", "Homo sapiens")] + [
            read(f"r{i}", "Bacteroides fragilis") for i in range(2, 12)
        ]
        flags = grade_classification("gut", reads)
        host = [f for f in flags if f.flag_type == "HOST_DNA_LEAK"]
        self.assertEqual(len(host), 1)
        self.assertEqual(host[0].severity, CRITICAL)
        self.assertEqual(host[0].read_id, "r1")

    def test_host_check_short_circuits_niche_check(self):
        # A human read is out of niche for every sample type. Reporting both
        # flags would be true twice over but would bury the serious finding.
        reads = [read("r1", "Homo sapiens")] + [
            read(f"r{i}", "Bacteroides fragilis") for i in range(2, 12)
        ]
        r1_flags = [f for f in grade_classification("gut", reads) if f.read_id == "r1"]
        self.assertEqual(len(r1_flags), 1)
        self.assertEqual(r1_flags[0].flag_type, "HOST_DNA_LEAK")

    def test_host_flagged_regardless_of_sample_type(self):
        for sample_type in ("gut", "skin", "soil", "unknown-niche"):
            with self.subTest(sample_type=sample_type):
                flags = grade_classification(sample_type, [read("r1", "Homo sapiens")])
                self.assertIn("HOST_DNA_LEAK", types_of(flags))


class TestNichePlausibility(unittest.TestCase):
    def test_expected_genus_is_not_flagged(self):
        reads = [read(f"r{i}", "Bacteroides fragilis") for i in range(10)]
        self.assertNotIn("UNEXPECTED_NICHE_TAXON", types_of(grade_classification("gut", reads)))

    def test_unexpected_genus_is_flagged_warning(self):
        reads = [read("r1", "Streptomyces coelicolor")] + [
            read(f"r{i}", "Bacteroides fragilis") for i in range(2, 12)
        ]
        flags = grade_classification("gut", reads)
        niche = [f for f in flags if f.flag_type == "UNEXPECTED_NICHE_TAXON"]
        self.assertEqual(len(niche), 1)
        self.assertEqual(niche[0].severity, WARNING)

    def test_unknown_sample_type_disables_niche_check(self):
        # Better to run no niche check than to flag every read in a niche we
        # have no expectations for.
        reads = [read(f"r{i}", "Streptomyces coelicolor") for i in range(10)]
        flags = grade_classification("marine-sediment", reads)
        self.assertNotIn("UNEXPECTED_NICHE_TAXON", types_of(flags))

    def test_sample_type_is_case_insensitive(self):
        reads = [read("r1", "Streptomyces coelicolor")] + [
            read(f"r{i}", "Bacteroides fragilis") for i in range(2, 12)
        ]
        self.assertIn("UNEXPECTED_NICHE_TAXON", types_of(grade_classification("GUT", reads)))

    def test_custom_expectations_override_defaults(self):
        reads = [read(f"r{i}", "Thermus aquaticus") for i in range(10)]
        flags = grade_classification(
            "hot-spring", reads, niche_expectations={"hot-spring": {"Thermus"}}
        )
        self.assertNotIn("UNEXPECTED_NICHE_TAXON", types_of(flags))


class TestRunLevelChecks(unittest.TestCase):
    def test_high_unclassified_rate_flagged_at_run_level(self):
        reads = [read(f"u{i}", "NA") for i in range(6)] + [
            read(f"c{i}", "Bacteroides fragilis") for i in range(4)
        ]
        flags = grade_classification("gut", reads)
        hit = [f for f in flags if f.flag_type == "HIGH_UNCLASSIFIED_RATE"]
        self.assertEqual(len(hit), 1)
        self.assertTrue(hit[0].is_run_level)
        self.assertEqual(hit[0].severity, CRITICAL)

    def test_unclassified_rate_at_threshold_is_not_flagged(self):
        # Threshold is exclusive: exactly 50% does not trip it.
        reads = [read(f"u{i}", "NA") for i in range(5)] + [
            read(f"c{i}", "Bacteroides fragilis") for i in range(5)
        ]
        self.assertNotIn("HIGH_UNCLASSIFIED_RATE", types_of(grade_classification("gut", reads)))

    def test_dominant_taxon_flagged(self):
        reads = [read(f"r{i}", "Bacteroides fragilis") for i in range(9)] + [
            read("r9", "Bifidobacterium longum")
        ]
        flags = grade_classification("gut", reads)
        hit = [f for f in flags if f.flag_type == "DOMINANT_TAXON"]
        self.assertEqual(len(hit), 1)
        self.assertTrue(hit[0].is_run_level)

    def test_even_community_is_not_flagged_as_dominant(self):
        reads = [read(f"a{i}", "Bacteroides fragilis") for i in range(5)] + [
            read(f"b{i}", "Bifidobacterium longum") for i in range(5)
        ]
        self.assertNotIn("DOMINANT_TAXON", types_of(grade_classification("gut", reads)))

    def test_thresholds_are_configurable(self):
        reads = [read(f"a{i}", "Bacteroides fragilis") for i in range(7)] + [
            read(f"b{i}", "Bifidobacterium longum") for i in range(3)
        ]
        self.assertNotIn("DOMINANT_TAXON", types_of(grade_classification("gut", reads)))
        relaxed = grade_classification("gut", reads, dominance_threshold=0.5)
        self.assertIn("DOMINANT_TAXON", types_of(relaxed))


class TestCleanInput(unittest.TestCase):
    def test_plausible_diverse_classified_sample_raises_nothing(self):
        reads = (
            [read(f"a{i}", "Bacteroides fragilis") for i in range(4)]
            + [read(f"b{i}", "Bifidobacterium longum") for i in range(3)]
            + [read(f"c{i}", "Escherichia coli") for i in range(3)]
        )
        self.assertEqual(grade_classification("gut", reads), [])


if __name__ == "__main__":
    unittest.main()
