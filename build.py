#!/usr/bin/env python3
"""Build the public, static GPCR Selectivity Atlas data bundle."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import shutil
import statistics
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
SITE = HERE / "site"
DATA = SITE / "data"
DOWNLOADS = SITE / "downloads"
STRUCTURE_DOWNLOADS = DOWNLOADS / "structures"
STRUCTURE_MANIFEST = STRUCTURE_DOWNLOADS / "manifest.json"

SOURCES = {
    "surface_assets": Path("data/masif/single_gpcr_outputs"),
    "surface_distance_matrix": Path("data/masif/classA286_distance_matrix.csv"),
    "surface_stats": Path("data/masif/_matrix_structure_stats.json"),
    "surface_mds": Path("data/masif/classA_mds_coords.csv"),
    "pairs_top3": Path("runs/batch_curated_masif_20260714/_top163_verified_regions.csv"),
    "seed_roster": Path("runs/agent6_bidirectional_seed_roster_post_wave4_304_20260730/seeds_304.csv"),
    "compounds_904": Path("deliverables/agent6_bidirectional_final_20260824/sampling_robust_improved_candidates.csv"),
    "admet": Path("runs/agent6_pocketxmol_bidirectional_filter244_admet_20260731/admet_predictions.csv"),
    "pose": Path("deliverables/agent6_bidirectional_final_20260824/pose_convergence.csv"),
    "pair_status": Path("deliverables/agent6_bidirectional_final_20260824/pairs_163_terminal_status.csv"),
    "direction_status": Path("deliverables/agent6_bidirectional_final_20260824/directions_326_terminal_status.csv"),
    "active_structure_refs": Path("result_analyse/slides_source_data/_allosteric_ref_gpcrdb_structures.csv"),
}

ADMET_FIELDS = {
    "ames": "admet_ai__AMES",
    "bbb": "admet_ai__BBB_Martins",
    "bioavailability": "admet_ai__Bioavailability_Ma",
    "cyp1a2_inhibition": "admet_ai__CYP1A2_Veith",
    "cyp2c19_inhibition": "admet_ai__CYP2C19_Veith",
    "cyp2c9_inhibition": "admet_ai__CYP2C9_Veith",
    "cyp2d6_inhibition": "admet_ai__CYP2D6_Veith",
    "cyp3a4_inhibition": "admet_ai__CYP3A4_Veith",
    "carcinogenicity": "admet_ai__Carcinogens_Lagunin",
    "clintox": "admet_ai__ClinTox",
    "dili": "admet_ai__DILI",
    "hia": "admet_ai__HIA_Hou",
    "pampa": "admet_ai__PAMPA_NCATS",
    "pgp_inhibition": "admet_ai__Pgp_Broccatelli",
    "herg": "admet_ai__hERG",
    "caco2": "admet_ai__Caco2_Wang",
    "hepatocyte_clearance": "admet_ai__Clearance_Hepatocyte_AZ",
    "microsome_clearance": "admet_ai__Clearance_Microsome_AZ",
    "half_life": "admet_ai__Half_Life_Obach",
    "solubility": "admet_ai__Solubility_AqSolDB",
    "vdss": "admet_ai__VDss_Lombardo",
}


def rows(rel: Path) -> list[dict[str, str]]:
    with (ROOT / rel).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def number(value: Any) -> float | int | None:
    if value is None or value == "":
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(parsed):
        return None
    return int(parsed) if parsed.is_integer() else parsed


def boolean(value: Any) -> bool | None:
    if value is None or value == "":
        return None
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    return None


def accession(part: str) -> str:
    match = re.match(r"([A-Z0-9]+)", part.strip())
    if not match:
        raise ValueError(f"Cannot parse UniProt accession from {part!r}")
    return match.group(1)


def directed_pair(pair_key: str) -> tuple[str, str]:
    left, right = pair_key.split("__", 1)
    return accession(left), accession(right)


def unordered_pair(a: str, b: str) -> str:
    return "__".join(sorted((a, b)))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dump_json(name: str, payload: Any) -> None:
    with (DATA / name).open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, separators=(",", ":"))


def write_csv(name: str, payload: Iterable[dict[str, Any]], fields: list[str]) -> None:
    with (DOWNLOADS / name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(payload)


def histogram(values: list[float], bins: int = 18) -> list[dict[str, float | int]]:
    if not values:
        return []
    low, high = min(values), max(values)
    width = (high - low) / bins if high > low else 1.0
    counts = [0] * bins
    for value in values:
        index = min(int((value - low) / width), bins - 1)
        counts[index] += 1
    return [
        {"start": round(low + i * width, 4), "end": round(low + (i + 1) * width, 4), "count": count}
        for i, count in enumerate(counts)
    ]


def build() -> dict[str, Any]:
    DATA.mkdir(parents=True, exist_ok=True)
    DOWNLOADS.mkdir(parents=True, exist_ok=True)

    compound_source = rows(SOURCES["compounds_904"])
    compound_ids = {row["prefilter_id"] for row in compound_source}
    if not STRUCTURE_MANIFEST.is_file():
        raise FileNotFoundError(
            f"Public structure manifest is missing: {STRUCTURE_MANIFEST}. "
            "Run prepare_structure_downloads.py against the frozen evidence host first."
        )
    structure_manifest = json.loads(STRUCTURE_MANIFEST.read_text(encoding="utf-8"))
    structure_by_id = {record["compound_id"]: record for record in structure_manifest["records"]}
    if set(structure_by_id) != compound_ids:
        raise RuntimeError("Public structure bundles do not exactly cover the 904-candidate master set")
    for compound_id, structure in structure_by_id.items():
        bundle_path = STRUCTURE_DOWNLOADS / structure["bundle_filename"]
        if not bundle_path.is_file() or sha256(bundle_path) != structure["bundle_sha256"]:
            raise RuntimeError(f"Public structure bundle is missing or hash-invalid: {compound_id}")
    pose_by_id = {row["prefilter_id"]: row for row in rows(SOURCES["pose"])}

    admet_by_id: dict[str, dict[str, str]] = {}
    with (ROOT / SOURCES["admet"]).open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            compound_id = row.get("prefilter_id", "")
            if compound_id in compound_ids:
                admet_by_id[compound_id] = row

    name_by_accession: dict[str, str] = {}
    top_rows = rows(SOURCES["pairs_top3"])
    for row in top_rows:
        name_by_accession[row["uniprot_a"]] = row["name_a"]
        name_by_accession[row["uniprot_b"]] = row["name_b"]
    for row in compound_source:
        for part in row["pair_key"].split("__"):
            fields = part.split("-", 1)
            if len(fields) == 2:
                name_by_accession.setdefault(fields[0], fields[1])

    pair_status_by_id = {row["original_pair_id"]: row for row in rows(SOURCES["pair_status"])}
    direction_status = rows(SOURCES["direction_status"])
    directions_by_pair: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in direction_status:
        directions_by_pair[row["original_pair_id"]].append(row)

    compounds: list[dict[str, Any]] = []
    compound_counts_by_pair: Counter[str] = Counter()
    compound_counts_by_seed_record: Counter[str] = Counter()
    seed_compound_ids: dict[str, list[str]] = defaultdict(list)
    for source in compound_source:
        compound_id = source["prefilter_id"]
        target, offtarget = directed_pair(source["pair_key"])
        pair_id = unordered_pair(target, offtarget)
        seed_record_id = f'{source["task"]}::{source["seed_zinc_id"]}'
        pose = pose_by_id.get(compound_id, {})
        admet = admet_by_id.get(compound_id, {})
        structure = structure_by_id[compound_id]
        properties = {
            key: number(admet.get(key) or source.get(key))
            for key in ["MW", "cLogP", "TPSA", "rotatable_bonds", "HBD", "HBA", "QED", "SA"]
        }
        admet_public = {name: number(admet.get(field)) for name, field in ADMET_FIELDS.items()}
        record = {
            "compound_id": compound_id,
            "canonical_smiles": source["canonical_smiles"],
            "task": source["task"],
            "branch": source["task"].rsplit("_", 1)[-1],
            "pair_id": pair_id,
            "target_uniprot": target,
            "offtarget_uniprot": offtarget,
            "target_name": name_by_accession.get(target),
            "offtarget_name": name_by_accession.get(offtarget),
            "seed_record_id": seed_record_id,
            "seed_zinc_id": source["seed_zinc_id"],
            "similarity_to_seed": number(source["ecfp4_tanimoto_to_seed"]),
            "structure_download": {
                "availability": "ligand_sdf_and_complex_pdb" if structure["complex_pdbs"] else "ligand_sdf_only",
                "bundle_url": f'downloads/structures/{structure["bundle_filename"]}',
                "bundle_size_bytes": structure["bundle_size_bytes"],
                "bundle_sha256": structure["bundle_sha256"],
                "ligand_sdf_count": 1,
                "complex_pdb_count": len(structure["complex_pdbs"]),
                "complex_pdbs": structure["complex_pdbs"],
                "note": "Computational receptor-ligand complexes from the audited MM/GBSA run; not experimental structures."
                if structure["complex_pdbs"] else
                "The ligand SDF is archived; no computational complex PDB is included in this bundle.",
            },
            "properties": properties,
            "docking": {
                "target_repeats": [number(source[f"target_seed_{seed}"]) for seed in ("17171", "29292", "43434")],
                "offtarget_repeats": [number(source[f"offtarget_seed_{seed}"]) for seed in ("17171", "29292", "43434")],
                "dd_repeats": [number(source[f"dd_seed_{seed}"]) for seed in ("17171", "29292", "43434")],
                "target_median": number(source["detail_target_median"]),
                "offtarget_median": number(source["detail_offtarget_median"]),
                "dd_median": number(source["detail_dd"]),
                "dd_worst": number(source["detail_dd_seed_max"]),
                "dd_sd": number(source["detail_dd_seed_sd"]),
                "seed_dd_median": number(source["seed_detail_dd"]),
                "dd_change_vs_seed": number(source["dd_change_vs_seed"]),
                "sampling_robust": boolean(source["dd_sampling_robust"]),
                "target_pose_stable": boolean(pose.get("target_pose_stable")),
                "offtarget_pose_stable": boolean(pose.get("offtarget_pose_stable")),
                "both_pose_stable": boolean(pose.get("both_pose_stable")),
            },
            "admet": admet_public,
            "evidence": {
                "final_candidate_904": True,
            },
        }
        compounds.append(record)
        compound_counts_by_pair[pair_id] += 1
        compound_counts_by_seed_record[seed_record_id] += 1
        seed_compound_ids[seed_record_id].append(compound_id)

    roster = rows(SOURCES["seed_roster"])
    roster_by_zinc: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in roster:
        roster_by_zinc[row["zinc_id"]].append(row)
    seed_groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in compound_source:
        seed_groups[f'{row["task"]}::{row["seed_zinc_id"]}'].append(row)

    seeds: list[dict[str, Any]] = []
    seed_counts_by_pair: Counter[str] = Counter()
    for seed_record_id, generated in sorted(seed_groups.items()):
        sample = generated[0]
        target, offtarget = directed_pair(sample["pair_key"])
        pair_id = unordered_pair(target, offtarget)
        candidates = roster_by_zinc[sample["seed_zinc_id"]]
        match = next(
            (r for r in candidates if r["target_uniprot"] == target and r["offtarget_uniprot"] == offtarget),
            candidates[0],
        )
        seeds.append(
            {
                "seed_record_id": seed_record_id,
                "seed_zinc_id": sample["seed_zinc_id"],
                "bidirectional_seed_id": match["bidirectional_seed_id"],
                "task": sample["task"],
                "pair_id": pair_id,
                "target_uniprot": target,
                "offtarget_uniprot": offtarget,
                "target_name": name_by_accession.get(target),
                "offtarget_name": name_by_accession.get(offtarget),
                "branch": sample["task"].rsplit("_", 1)[-1],
                "hotspot_bw": match["bw_label"],
                "fast": {
                    "target": number(match["fast_target"]),
                    "offtarget": number(match["fast_offtarget"]),
                    "dd": number(match["fast_dd"]),
                    "rank": number(match["rank"]),
                },
                "detail": {
                    "target_median": number(match["detail_target_median"]),
                    "offtarget_median": number(match["detail_offtarget_median"]),
                    "dd_median": number(match["detail_dd"]),
                    "dd_worst": number(match["detail_dd_seed_max"]),
                    "dd_sd": number(match["detail_dd_seed_sd"]),
                    "target_pose_stable": boolean(match["target_pose_stable"]),
                    "offtarget_pose_stable": boolean(match["offtarget_pose_stable"]),
                },
                "generated_compound_count": compound_counts_by_seed_record[seed_record_id],
                "generated_compound_ids": sorted(seed_compound_ids[seed_record_id]),
            }
        )
        seed_counts_by_pair[pair_id] += 1

    pairs: list[dict[str, Any]] = []
    for row in top_rows:
        pair_id = unordered_pair(row["uniprot_a"], row["uniprot_b"])
        bws = row["top3_bw"].split("|")
        residues = row["top3_residues"].split("|")
        diffs = row["top3_fp_diff"].split("|")
        hotspots = [
            {"hotspot_id": f"{pair_id}::H{i + 1}", "rank": i + 1, "bw": bws[i], "residues": residues[i], "fingerprint_difference": number(diffs[i])}
            for i in range(3)
        ]
        status = pair_status_by_id.get(pair_id, {})
        directions = [
            {
                "direction": d["direction"],
                "task_id": d["directed_task_id"],
                "seed_id": d["bidirectional_seed_id"],
                "seed_zinc_id": d["seed_zinc_id"],
                "terminal_status": d["terminal_status"],
                "redocking_candidates": number(d["redocking_candidates"]),
            }
            for d in directions_by_pair.get(pair_id, [])
        ]
        pairs.append(
            {
                "pair_id": pair_id,
                "rank": number(row["rank"]),
                "receptor_a": {"uniprot": row["uniprot_a"], "name": row["name_a"]},
                "receptor_b": {"uniprot": row["uniprot_b"], "name": row["name_b"]},
                "surface_distance": number(row["distance"]),
                "hotspots": hotspots,
                "directions": directions,
                "pair_terminal_status": status.get("pair_terminal_status"),
                "final_candidate_count": compound_counts_by_pair[pair_id],
                "pocketxmol_seed_record_count": seed_counts_by_pair[pair_id],
                "pocketxmol_compound_count": compound_counts_by_pair[pair_id],
                "input_seed_zinc_ids": sorted({
                    seed["seed_zinc_id"] for seed in seeds
                    if seed["pair_id"] == pair_id and seed["seed_zinc_id"]
                }),
            }
        )

    matrix_path = ROOT / SOURCES["surface_distance_matrix"]
    with matrix_path.open(newline="", encoding="utf-8-sig") as handle:
        matrix_rows = list(csv.reader(handle))
    matrix_accessions = matrix_rows[0][2:]
    matrix_values: list[float] = []
    nearest: dict[str, dict[str, Any]] = {}
    subfamily: dict[str, str | None] = {}
    for i, row in enumerate(matrix_rows[1:]):
        receptor = row[0]
        subfamily[receptor] = None if row[1] == "?" else row[1]
        values = [float(value) for value in row[2:]]
        candidates = [(value, matrix_accessions[j]) for j, value in enumerate(values) if j != i]
        best_distance, best_receptor = min(candidates)
        nearest[receptor] = {"uniprot": best_receptor, "distance": round(best_distance, 4)}
        matrix_values.extend(value for j, value in enumerate(values) if j > i)

    mds_rows = rows(SOURCES["surface_mds"])
    mds_by_id = {
        row["uniprot"]: {"x": number(row["MDS1"]), "y": number(row["MDS2"]), "subfamily": None if row["subfamily"] == "?" else row["subfamily"]}
        for row in mds_rows
    }

    pair_participation: Counter[str] = Counter()
    for pair in pairs:
        pair_participation[pair["receptor_a"]["uniprot"]] += 1
        pair_participation[pair["receptor_b"]["uniprot"]] += 1

    receptors: list[dict[str, Any]] = []
    surface_root = ROOT / SOURCES["surface_assets"]
    for folder in sorted(path for path in surface_root.iterdir() if path.is_dir()):
        receptor = folder.name
        assets = []
        for file in sorted(path for path in folder.iterdir() if path.is_file()):
            assets.append({"filename": file.name, "kind": file.suffix.lstrip("."), "size_bytes": file.stat().st_size})
        expected = {
            f"{receptor}_A_predcoords.npy",
            f"{receptor}_A_predfeatures_emb1.npy",
            f"{receptor}_A_predfeatures_emb2.npy",
            f"{receptor}_A_pred_emb1.vtk",
            f"{receptor}_A_pred_emb2.vtk",
        }
        present = {asset["filename"] for asset in assets}
        receptors.append(
            {
                "uniprot": receptor,
                "name": name_by_accession.get(receptor),
                "subfamily": subfamily.get(receptor),
                "surface_asset_count": len(assets),
                "surface_asset_bytes": sum(asset["size_bytes"] for asset in assets),
                "core_assets_complete": expected.issubset(present),
                "assets": assets,
                "nearest_surface_neighbor": nearest.get(receptor),
                "mds": mds_by_id.get(receptor),
                "selected_pair_count": pair_participation[receptor],
                "public_raw_asset_status": "pending_external_archive",
            }
        )

    with (ROOT / SOURCES["surface_stats"]).open(encoding="utf-8") as handle:
        surface_stats = json.load(handle)
    global_surface = {
        "matrix_receptor_count": len(matrix_accessions),
        "matrix_pair_count": len(matrix_values),
        "distance_summary": {
            "minimum": round(min(matrix_values), 4),
            "median": round(statistics.median(matrix_values), 4),
            "mean": round(statistics.fmean(matrix_values), 4),
            "maximum": round(max(matrix_values), 4),
        },
        "distance_histogram": histogram(matrix_values),
        "selected_pair_distance_histogram": histogram([float(pair["surface_distance"]) for pair in pairs]),
        "structure_stats": surface_stats,
        "mds": [{"uniprot": key, **value, "selected_pair_count": pair_participation[key]} for key, value in mds_by_id.items()],
        "mds_coverage": len(mds_by_id),
        "mds_note": "Legacy frozen coordinates cover 275 of 286 matrix receptors; the full matrix remains the canonical global-distance evidence.",
        "download": "downloads/classA286_distance_matrix.csv",
        "interpretation": "dMaSIF fingerprint distance is a learned surface-representation distance, not sequence identity, experimental affinity, or active-state structural RMSD.",
    }

    active_structures = [
        {
            "pdb_id": row["pdb"],
            "protein": row["protein"],
            "ligand": row["ligand"],
            "function": row["function"],
            "bw_site_count": number(row["n_bw"]),
            "bw_sites": row["bws"].split("|") if row["bws"] else [],
        }
        for row in rows(SOURCES["active_structure_refs"])
    ]
    active_payload = {
        "pairwise_similarity_status": "not_computed",
        "message": "The project currently has no frozen, systematic pairwise active-state GPCR similarity matrix. These 36 records are structure references only.",
        "future_metrics": ["TM-score", "C-alpha RMSD", "binding-pocket RMSD", "shared BW positions"],
        "structures": active_structures,
    }

    summary = {
        "title": "GPCR Selectivity Atlas",
        "version": "0.1.0-mvp",
        "build_date": date.today().isoformat(),
        "creators": [
            {
                "name": "Zhu, Jingyi",
                "name_zh": "朱景一",
                "affiliation": "Shandong University",
                "affiliation_zh": "山东大学",
                "orcid": "0009-0003-8404-0455",
            }
        ],
        "counts": {
            "surface_receptors": len(receptors),
            "surface_matrix_receptors": len(matrix_accessions),
            "selected_receptor_pairs": len(pairs),
            "hotspots": sum(len(pair["hotspots"]) for pair in pairs),
            "directed_tasks": len(direction_status),
            "unique_seed_molecules": len({seed["seed_zinc_id"] for seed in seeds}),
            "seed_task_records": len(seeds),
            "generated_compounds": len(compounds),
            "final_candidates": len(compounds),
            "admet_complete": sum(all(value is not None for value in compound["admet"].values()) for compound in compounds),
            "structure_bundles": structure_manifest["counts"]["compounds"],
            "ligand_sdf_files": structure_manifest["counts"]["ligand_sdf_files"],
            "complex_pdb_files": structure_manifest["counts"]["complex_pdb_files"],
            "active_structure_references": len(active_structures),
            "active_structure_pairwise_similarities": 0,
        },
        "evidence_scope": "computational_prioritization_only",
    }

    source_hashes = {}
    for name, rel in SOURCES.items():
        path = ROOT / rel
        if path.is_file():
            source_hashes[name] = {"source": rel.as_posix(), "sha256": sha256(path)}
        else:
            source_hashes[name] = {"source": rel.as_posix(), "record_type": "directory_manifest"}
    provenance = {
        "atlas_version": summary["version"],
        "build_date": summary["build_date"],
        "creators": summary["creators"],
        "sources": {
            **source_hashes,
            "public_structure_manifest": {
                "source": "site/downloads/structures/manifest.json",
                "sha256": sha256(STRUCTURE_MANIFEST),
            },
        },
        "limitations": [
            "All docking, ADMET, pose-stability, dMaSIF and MM/GBSA values are computational predictions.",
            "The active-state structural similarity matrix has not yet been computed.",
            "Large dMaSIF binary assets remain local until an external public archive and redistribution terms are confirmed.",
            "The frozen MDS table covers 275 receptors while the canonical distance matrix covers 286.",
        ],
    }

    dump_json("summary.json", summary)
    dump_json("receptors.json", receptors)
    dump_json("pairs.json", pairs)
    dump_json("seeds.json", seeds)
    dump_json("compounds.json", compounds)
    dump_json("global_surface.json", global_surface)
    dump_json("active_structures.json", active_payload)
    dump_json("provenance.json", provenance)

    receptor_flat = [
        {
            "uniprot": r["uniprot"], "name": r["name"], "subfamily": r["subfamily"],
            "surface_asset_count": r["surface_asset_count"], "surface_asset_bytes": r["surface_asset_bytes"],
            "core_assets_complete": r["core_assets_complete"], "nearest_neighbor": (r["nearest_surface_neighbor"] or {}).get("uniprot"),
            "nearest_distance": (r["nearest_surface_neighbor"] or {}).get("distance"), "selected_pair_count": r["selected_pair_count"],
            "public_raw_asset_status": r["public_raw_asset_status"],
        }
        for r in receptors
    ]
    write_csv("receptors.csv", receptor_flat, list(receptor_flat[0]))
    pair_flat = []
    for pair in pairs:
        out = {
            "pair_id": pair["pair_id"], "rank": pair["rank"],
            "uniprot_a": pair["receptor_a"]["uniprot"], "name_a": pair["receptor_a"]["name"],
            "uniprot_b": pair["receptor_b"]["uniprot"], "name_b": pair["receptor_b"]["name"],
            "surface_distance": pair["surface_distance"], "pair_terminal_status": pair["pair_terminal_status"],
            "seed_task_records": pair["pocketxmol_seed_record_count"], "input_seed_zinc_ids": "|".join(pair["input_seed_zinc_ids"]),
            "generated_compounds": pair["pocketxmol_compound_count"], "final_candidate_count": pair["final_candidate_count"],
        }
        for i, hotspot in enumerate(pair["hotspots"], 1):
            out[f"hotspot{i}_bw"] = hotspot["bw"]
            out[f"hotspot{i}_residues"] = hotspot["residues"]
            out[f"hotspot{i}_fp_diff"] = hotspot["fingerprint_difference"]
        pair_flat.append(out)
    write_csv("receptor_pairs_top3.csv", pair_flat, list(pair_flat[0]))
    seed_flat = [
        {
            "seed_record_id": s["seed_record_id"], "seed_zinc_id": s["seed_zinc_id"], "bidirectional_seed_id": s["bidirectional_seed_id"],
            "task": s["task"], "pair_id": s["pair_id"], "target_uniprot": s["target_uniprot"], "offtarget_uniprot": s["offtarget_uniprot"],
            "branch": s["branch"], "hotspot_bw": s["hotspot_bw"], "fast_target": s["fast"]["target"], "fast_offtarget": s["fast"]["offtarget"],
            "fast_dd": s["fast"]["dd"], "detail_target_median": s["detail"]["target_median"], "detail_offtarget_median": s["detail"]["offtarget_median"],
            "detail_dd_median": s["detail"]["dd_median"], "detail_dd_worst": s["detail"]["dd_worst"], "detail_dd_sd": s["detail"]["dd_sd"],
            "generated_compound_count": s["generated_compound_count"],
        }
        for s in seeds
    ]
    write_csv("seeds.csv", seed_flat, list(seed_flat[0]))
    compound_flat = []
    for c in compounds:
        out = {
            "compound_id": c["compound_id"], "canonical_smiles": c["canonical_smiles"], "task": c["task"], "branch": c["branch"],
            "pair_id": c["pair_id"], "target_uniprot": c["target_uniprot"], "offtarget_uniprot": c["offtarget_uniprot"],
            "seed_record_id": c["seed_record_id"], "seed_zinc_id": c["seed_zinc_id"], "similarity_to_seed": c["similarity_to_seed"],
            **c["properties"],
            "detail_target_median": c["docking"]["target_median"], "detail_offtarget_median": c["docking"]["offtarget_median"],
            "detail_dd_median": c["docking"]["dd_median"], "detail_dd_worst": c["docking"]["dd_worst"], "detail_dd_sd": c["docking"]["dd_sd"],
            "seed_detail_dd_median": c["docking"]["seed_dd_median"], "dd_change_vs_seed": c["docking"]["dd_change_vs_seed"],
            "target_pose_stable": c["docking"]["target_pose_stable"], "offtarget_pose_stable": c["docking"]["offtarget_pose_stable"],
            **c["admet"],
            "final_candidate_904": c["evidence"]["final_candidate_904"],
            "structure_availability": c["structure_download"]["availability"],
            "structure_bundle_url": c["structure_download"]["bundle_url"],
            "complex_pdb_count": c["structure_download"]["complex_pdb_count"],
        }
        compound_flat.append(out)
    write_csv("compounds_904.csv", compound_flat, list(compound_flat[0]))
    write_csv("active_structure_references.csv", active_structures, ["pdb_id", "protein", "ligand", "function", "bw_site_count", "bw_sites"])
    shutil.copyfile(matrix_path, DOWNLOADS / "classA286_distance_matrix.csv")
    return summary


if __name__ == "__main__":
    print(json.dumps(build(), indent=2, ensure_ascii=False))
