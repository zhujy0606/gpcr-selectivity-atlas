#!/usr/bin/env python3
"""Create public, per-compound structure bundles from frozen Agent6 evidence.

This script is intended to run at the project root on the evidence host.  It
never edits a frozen run.  Each output ZIP contains the candidate ligand SDF
and, when available, the MM/GBSA receptor-ligand complex PDB files from the
audited 438-candidate run.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import tarfile
import tempfile
import zipfile
from collections import defaultdict
from datetime import date
from pathlib import Path


MASTER_CANDIDATES = (
    Path("deliverables/agent6_bidirectional_final_20260824/sampling_robust_improved_candidates.csv"),
    Path("runs/agent6_pocketxmol_bidirectional_redocking_analysis_20260731/sampling_robust_improved_candidates.csv"),
)
ADMET = Path("runs/agent6_pocketxmol_bidirectional_filter244_admet_20260731/admet_predictions.csv")
LIGAND_ROOT_CANDIDATES = (
    Path("runs/agent6_pocketxmol_bidirectional_filter244_admet_20260731"),
    Path("runs/agent6_pocketxmol_bidirectional_docking_selection244_20260731"),
    Path("runs/agent6_pocketxmol_production_filter88_20260726"),
)
SHORTLIST_CANDIDATES = (
    Path("deliverables/agent6_bidirectional_final_20260824/cross_domain_shortlist.csv"),
    Path("runs/agent6_pocketxmol_bidirectional_production_pareto_20260731/cross_domain_shortlist.csv"),
)
COMPLEX_RUN = Path("runs/agent6_pocketxmol_bidirectional_candidate_mmgbsa438_scfrescue_topologyextract_contract_20260805")
COMPLEX_AUDIT = Path("runs/agent6_pocketxmol_bidirectional_candidate_mmgbsa438_scfrescue_topologyextract_contract_audit_20260809")
EXPECTED_COMPOUNDS = 904
EXPECTED_MMGBSA_COMPOUNDS = 438
EXPECTED_COMPLEXES = 1444


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def resolve_source(repo: Path, candidates: tuple[Path, ...]) -> Path:
    for candidate in candidates:
        if (repo / candidate).is_file():
            return candidate
    raise FileNotFoundError(f"None of the frozen source paths exists: {candidates}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def add_file(archive: zipfile.ZipFile, source: Path, member: str) -> dict[str, object]:
    archive.write(source, member, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    return {
        "member": member,
        "format": source.suffix.lstrip(".").upper(),
        "size_bytes": source.stat().st_size,
        "sha256": sha256(source),
    }


def parse_job(job_name: str) -> dict[str, object]:
    match = re.search(r"__(target|offtarget)__c(\d+)__([A-Z0-9]+)$", job_name)
    if not match:
        raise ValueError(f"Unrecognized MM/GBSA job name: {job_name}")
    return {"receptor_role": match.group(1), "cluster": int(match.group(2)), "receptor_uniprot": match.group(3)}


def build(repo: Path, output_tar: Path) -> dict[str, object]:
    repo = repo.resolve()
    master = resolve_source(repo, MASTER_CANDIDATES)
    shortlist = resolve_source(repo, SHORTLIST_CANDIDATES)
    master_rows = read_rows(repo / master)
    compound_ids = {row["prefilter_id"] for row in master_rows}
    if len(compound_ids) != EXPECTED_COMPOUNDS:
        raise RuntimeError(f"Expected {EXPECTED_COMPOUNDS} final candidates, found {len(compound_ids)}")

    shortlist_ids = {row["prefilter_id"] for row in read_rows(repo / shortlist)}
    if len(shortlist_ids) != EXPECTED_MMGBSA_COMPOUNDS or not shortlist_ids <= compound_ids:
        raise RuntimeError("Frozen 438-candidate shortlist does not match the 904-candidate master set")

    ligand_by_id: dict[str, tuple[Path, str]] = {}
    for row in read_rows(repo / ADMET):
        compound_id = row.get("prefilter_id", "")
        if compound_id not in compound_ids:
            continue
        expected = row["copied_sdf_sha256"]
        attempted = [repo / source_root / row["copied_sdf_path"] for source_root in LIGAND_ROOT_CANDIDATES]
        source = next((path for path in attempted if path.is_file() and sha256(path) == expected), None)
        if source is None:
            raise FileNotFoundError(f"No hash-matching ligand SDF for {compound_id}; checked {attempted}")
        observed = sha256(source)
        if observed != expected:
            raise RuntimeError(f"Ligand SDF hash mismatch for {compound_id}: {observed} != {expected}")
        ligand_by_id[compound_id] = (source, expected)
    if set(ligand_by_id) != compound_ids:
        missing = sorted(compound_ids - set(ligand_by_id))
        raise RuntimeError(f"Missing indexed ligand SDF files: {missing[:10]}")

    complexes_by_id: dict[str, list[Path]] = defaultdict(list)
    for path in sorted((repo / COMPLEX_RUN / "jobs").glob("*/complex_leap.pdb")):
        compound_id = path.parent.name.split("__", 1)[0]
        if compound_id in shortlist_ids:
            complexes_by_id[compound_id].append(path)
    complex_count = sum(len(paths) for paths in complexes_by_id.values())
    if set(complexes_by_id) != shortlist_ids or complex_count != EXPECTED_COMPLEXES:
        raise RuntimeError(
            f"Expected {EXPECTED_MMGBSA_COMPOUNDS}/{EXPECTED_COMPLEXES} complex coverage, "
            f"found {len(complexes_by_id)}/{complex_count}"
        )

    with tempfile.TemporaryDirectory(prefix="agent6_public_structures_") as temp_name:
        stage = Path(temp_name) / "structures"
        stage.mkdir()
        records = []
        for compound_id in sorted(compound_ids):
            zip_path = stage / f"{compound_id}.zip"
            with zipfile.ZipFile(zip_path, "w") as archive:
                ligand_source, ligand_hash = ligand_by_id[compound_id]
                ligand = add_file(archive, ligand_source, f"ligand/{compound_id}.sdf")
                if ligand["sha256"] != ligand_hash:
                    raise RuntimeError(f"Unexpected ligand hash drift for {compound_id}")
                complexes = []
                for pdb_path in sorted(complexes_by_id.get(compound_id, [])):
                    job_name = pdb_path.parent.name
                    item = add_file(archive, pdb_path, f"complexes/{job_name}.pdb")
                    item.update(parse_job(job_name))
                    complexes.append(item)
                readme = (
                    f"GPCR Selectivity Atlas structure bundle: {compound_id}\n"
                    f"Candidate ligand: ligand/{compound_id}.sdf\n"
                    f"Audited MM/GBSA complex PDB files: {len(complexes)}\n"
                    "PDB files are computational receptor-ligand complexes, not experimental structures.\n"
                )
                archive.writestr("README.txt", readme)
            records.append(
                {
                    "compound_id": compound_id,
                    "bundle_filename": zip_path.name,
                    "bundle_size_bytes": zip_path.stat().st_size,
                    "bundle_sha256": sha256(zip_path),
                    "ligand_sdf": ligand,
                    "complex_pdbs": complexes,
                }
            )

        manifest = {
            "schema_version": "1.0",
            "build_date": date.today().isoformat(),
            "scope": "904 candidate ligand SDF files; computational receptor-ligand complex PDB files when available",
            "counts": {
                "compounds": len(records),
                "ligand_sdf_files": len(records),
                "complex_pdb_files": sum(len(record["complex_pdbs"]) for record in records),
            },
            "sources": {
                "candidate_master": master.as_posix(),
                "ligand_index": ADMET.as_posix(),
            },
            "records": records,
        }
        (stage / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
        output_tar.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(output_tar, "w:gz", compresslevel=1) as archive:
            archive.add(stage, arcname="structures")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = build(args.repo, args.output)
    print(json.dumps({"output": str(args.output), **manifest["counts"]}, indent=2))


if __name__ == "__main__":
    main()
