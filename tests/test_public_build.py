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
            "admet_complete": 904,
            "cross_domain_shortlist": 438,
            "mmgbsa_seed_baseline_complete": 323,
            "mmgbsa_dual_endpoint_improved": 141,
            "final_selected": 111,
            "strict_final_selected": 111,
            "structure_bundles": 904,
            "ligand_sdf_files": 904,
            "compounds_with_complex_pdb": 438,
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
        self.assertEqual(sum(pair["strict_final_selected_count"] for pair in self.pairs), 111)

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

    def test_public_interface_is_two_module_table_catalog(self):
        html = (SITE / "index.html").read_text()
        javascript = (SITE / "app.js").read_text()
        pages = re.findall(r'<section id="([^"]+)" class="page', html)
        self.assertEqual(pages, ["overview", "receptors", "pairs"])
        self.assertIn("287受体 dMaSIF 数据", html)
        self.assertIn("163对受体与选择性分子", html)
        self.assertIn("Top 3热点", javascript)
        self.assertIn("输入种子", javascript)
        self.assertIn("robust生成分子", javascript)
        self.assertIn("structure.bundle_url", javascript)
        self.assertNotIn("严格精选", html + javascript)
        self.assertNotIn("pair.strict_final_selected_count", javascript)
        self.assertNotIn("evidenceLabel", javascript)
        self.assertIn("'ΔDD vs seed','结构下载'", javascript)

    def test_structure_downloads(self):
        self.assertTrue(all(compound["structure_download"]["ligand_sdf_count"] == 1 for compound in self.compounds))
        self.assertEqual(sum(compound["structure_download"]["complex_pdb_count"] > 0 for compound in self.compounds), 438)
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
