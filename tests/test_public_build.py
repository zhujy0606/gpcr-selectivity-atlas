import csv
import json
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PUBLIC = ROOT / "public_database"
SITE = PUBLIC / "site"


class PublicBuildTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sys
        sys.path.insert(0, str(PUBLIC))
        import build
        cls.summary = build.build()
        cls.compounds = json.loads((SITE / "data/compounds.json").read_text())
        cls.pairs = json.loads((SITE / "data/pairs.json").read_text())
        cls.receptors = json.loads((SITE / "data/receptors.json").read_text())

    def test_frozen_counts(self):
        counts = self.summary["counts"]
        expected = {
            "surface_receptors": 287,
            "surface_matrix_receptors": 286,
            "selected_receptor_pairs": 163,
            "hotspots": 489,
            "directed_tasks": 326,
            "unique_seed_molecules": 138,
            "seed_task_records": 178,
            "generated_compounds": 904,
            "admet_complete": 904,
            "cross_domain_shortlist": 438,
            "mmgbsa_seed_baseline_complete": 323,
            "mmgbsa_dual_endpoint_improved": 141,
            "final_selected": 111,
        }
        for key, value in expected.items():
            self.assertEqual(counts[key], value, key)

    def test_relations(self):
        self.assertEqual(len(self.receptors), 287)
        self.assertEqual(len(self.pairs), 163)
        self.assertTrue(all(len(pair["hotspots"]) == 3 for pair in self.pairs))
        self.assertEqual(len({compound["compound_id"] for compound in self.compounds}), 904)
        self.assertTrue(all(compound["seed_record_id"] and compound["pair_id"] for compound in self.compounds))

    def test_distance_matrix(self):
        with (SITE / "downloads/classA286_distance_matrix.csv").open(newline="") as handle:
            table = list(csv.reader(handle))
        self.assertEqual(len(table) - 1, 286)
        self.assertEqual(len(table[0]) - 2, 286)

    def test_no_sensitive_paths_or_credentials(self):
        forbidden = re.compile(r"/Users/|/home/zkyd/|10\.102\.106\.53|zkyd@12", re.I)
        for path in SITE.rglob("*"):
            if path.is_file() and path.suffix.lower() in {".json", ".csv", ".html", ".css", ".js"}:
                self.assertIsNone(forbidden.search(path.read_text(errors="ignore")), str(path))


if __name__ == "__main__":
    unittest.main()
