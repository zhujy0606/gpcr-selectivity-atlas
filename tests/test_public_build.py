import csv
import json
import re
import unittest
import zipfile
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
        cls.detail_mode_compounds = json.loads((SITE / "data/detail_mode_compounds.json").read_text())

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
            "final_candidates": 904,
            "zero_generation_receptor_pairs": 48,
            "detail_mode_selected_compounds_for_zero_generation_pairs": 87,
            "zero_generation_pairs_with_detail_mode_compounds": 48,
            "admet_complete": 904,
            "structure_bundles": 904,
            "ligand_sdf_files": 904,
            "complex_pdb_files": 1444,
        }
        for key, value in expected.items():
            self.assertEqual(counts[key], value, key)

    def test_creator_metadata(self):
        creator = self.summary["creators"][0]
        self.assertEqual(creator["name"], "Zhu, Jingyi")
        self.assertEqual(creator["affiliation"], "Shandong University")
        self.assertEqual(creator["orcid"], "0009-0003-8404-0455")
        citation = (PUBLIC / "CITATION.cff").read_text()
        self.assertIn('given-names: "Jingyi"', citation)
        self.assertIn('orcid: "https://orcid.org/0009-0003-8404-0455"', citation)

    def test_relations(self):
        self.assertEqual(len(self.receptors), 287)
        self.assertEqual(len(self.pairs), 163)
        self.assertTrue(all(len(pair["hotspots"]) == 3 for pair in self.pairs))
        self.assertEqual(len({compound["compound_id"] for compound in self.compounds}), 904)
        self.assertTrue(all(compound["seed_record_id"] and compound["pair_id"] for compound in self.compounds))
        self.assertTrue(all(compound["evidence"]["final_candidate_904"] for compound in self.compounds))
        self.assertTrue(all(compound["seed_zinc_id"].startswith("ZINC") for compound in self.compounds))
        self.assertEqual(sum(pair["final_candidate_count"] for pair in self.pairs), 904)
        self.assertTrue(all(set(compound["evidence"]) == {"final_candidate_904"} for compound in self.compounds))
        self.assertTrue(all("strict_final_selected_count" not in pair for pair in self.pairs))

        pair_ids = {pair["pair_id"] for pair in self.pairs}
        self.assertTrue(all(compound["pair_id"] in pair_ids for compound in self.compounds))
        linked_counts = {
            pair_id: sum(compound["pair_id"] == pair_id for compound in self.compounds)
            for pair_id in pair_ids
        }
        self.assertTrue(all(
            linked_counts[pair["pair_id"]] == pair["final_candidate_count"]
            for pair in self.pairs
        ))

        zero_generation_pairs = {
            pair["pair_id"] for pair in self.pairs
            if pair["pocketxmol_compound_count"] == 0
        }
        detail_pairs = {record["pair_id"] for record in self.detail_mode_compounds}
        self.assertEqual(len(zero_generation_pairs), 48)
        self.assertEqual(len(self.detail_mode_compounds), 87)
        self.assertEqual(detail_pairs, zero_generation_pairs)
        self.assertTrue(all(record["zinc_id"].startswith("ZINC") for record in self.detail_mode_compounds))
        self.assertTrue(all(record["pocketxmol_generated"] is False for record in self.detail_mode_compounds))
        self.assertTrue(all(1 <= pair["detail_mode_selected_compound_count"] <= 2 for pair in self.pairs if pair["pair_id"] in zero_generation_pairs))
        self.assertTrue(all(pair["detail_mode_selected_compound_count"] == 0 for pair in self.pairs if pair["pair_id"] not in zero_generation_pairs))

    def test_public_interface_is_two_module_table_catalog(self):
        html = (SITE / "index.html").read_text()
        javascript = (SITE / "app.js").read_text()
        pages = re.findall(r'<section id="([^"]+)" class="page', html)
        self.assertEqual(pages, ["overview", "receptors", "pairs"])
        self.assertIn("dMaSIF receptor surfaces", html)
        self.assertIn("Receptor-pair selectivity evidence", html)
        self.assertIn("Top 3 hotspots", javascript)
        self.assertIn("Input seeds", javascript)
        self.assertIn("Pocketxmol-generated compounds", html + javascript)
        self.assertIn("Detail-mode selected compounds", html + javascript)
        self.assertIn("detail_mode_compounds", javascript)
        self.assertIsNone(re.search(r"[\u4e00-\u9fff]", html + javascript))
        self.assertNotIn("PocketXMol", html + javascript)
        self.assertNotIn("robust", html.lower() + javascript.lower())
        self.assertNotIn("438", html + javascript)
        self.assertIn("structure.bundle_url", javascript)
        self.assertNotIn("pair.strict_final_selected_count", javascript)
        self.assertNotIn("evidenceLabel", javascript)
        self.assertIn("'Detail ΔE (kcal/mol)'", javascript)
        self.assertIn("'Worst ΔE (kcal/mol)'", javascript)
        self.assertIn("'ΔE change vs seed (kcal/mol)'", javascript)
        self.assertIn("ΔE = E<sub>target</sub> − E<sub>off-target</sub>", javascript)

    def test_structure_downloads(self):
        self.assertTrue(all(compound["structure_download"]["ligand_sdf_count"] == 1 for compound in self.compounds))
        self.assertGreater(sum(compound["structure_download"]["complex_pdb_count"] > 0 for compound in self.compounds), 0)
        self.assertEqual(sum(compound["structure_download"]["complex_pdb_count"] for compound in self.compounds), 1444)
        for compound in self.compounds:
            structure = compound["structure_download"]
            bundle = SITE / structure["bundle_url"]
            self.assertTrue(bundle.is_file(), bundle)
            self.assertEqual(bundle.name, f'{compound["compound_id"]}.zip')
        with_pdb = next(c for c in self.compounds if c["structure_download"]["complex_pdb_count"])
        without_pdb = next(c for c in self.compounds if not c["structure_download"]["complex_pdb_count"])
        for compound, expect_pdb in ((with_pdb, True), (without_pdb, False)):
            with zipfile.ZipFile(SITE / compound["structure_download"]["bundle_url"]) as archive:
                members = archive.namelist()
                self.assertIn(f'ligand/{compound["compound_id"]}.sdf', members)
                self.assertEqual(any(name.startswith("complexes/") and name.endswith(".pdb") for name in members), expect_pdb)

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
