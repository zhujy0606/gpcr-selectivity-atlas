"""Validate and join endpoint rescoring evidence without changing the cohort."""

from collections import defaultdict
import math
import statistics


def unique_rows(rows, key, label):
    result = {}
    for row in rows:
        value = row[key]
        if not value or value in result:
            raise ValueError(f"Missing or duplicate {label}: {value!r}")
        result[value] = row
    return result


def finite(value):
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError("Nonfinite MM-GBSA energy")
    return parsed


def flag(value):
    if value not in ("True", "False"):
        raise ValueError(f"Invalid pose-stability flag: {value!r}")
    return value == "True"


def check_equal(actual, expected, label):
    if actual != expected:
        raise ValueError(f"MM-GBSA {label} mismatch: {actual!r} != {expected!r}")


def check_energy(actual, expected, label):
    if not math.isclose(finite(actual), expected, rel_tol=0, abs_tol=1e-8):
        raise ValueError(f"MM-GBSA {label} arithmetic mismatch")


def load_mmgbsa(masters, paired, inventory, endpoints):
    """Return candidate, actual-parent and sanitized pose records after validation.

    Pair direction and task are mandatory joins. A seed reused in another task
    must retain its own baseline. Pose stability is a flag, not an exclusion.
    """
    master_by_id = unique_rows(masters, "prefilter_id", "master ID")
    pairs = unique_rows(paired, "prefilter_id", "paired candidate ID")
    molecules = unique_rows(inventory, "prefilter_id", "molecule ID")
    check_equal(set(pairs), set(master_by_id), "candidate coverage")
    if any(row["kind"] not in ("candidate", "seed") for row in inventory):
        raise ValueError("Unknown MM-GBSA molecule kind")
    check_equal({key for key, row in molecules.items() if row["kind"] == "candidate"},
                set(master_by_id), "inventory coverage")

    energies = defaultdict(lambda: defaultdict(list))
    seen_poses = set()
    public_poses = []
    for pose in endpoints:
        molecule_id, role = pose["prefilter_id"], pose["role"]
        if molecule_id not in molecules or role not in ("target", "offtarget"):
            raise ValueError("Unknown MM-GBSA pose molecule or receptor role")
        molecule = molecules[molecule_id]
        target, offtarget = (part.split("-", 1)[0] for part in molecule["pair_key"].split("__"))
        check_equal(pose["task"], molecule["task"], "pose task")
        check_equal(pose["receptor_uniprot"], target if role == "target" else offtarget,
                    "pose receptor direction")
        rank = int(pose["cluster_rank"])
        key = (molecule_id, role, rank)
        if rank < 1 or key in seen_poses:
            raise ValueError("Invalid or duplicate MM-GBSA pose cluster")
        seen_poses.add(key)
        energy = finite(pose["delta_total_kcal_mol"])
        energies[molecule_id][role].append(energy)
        # Explicit allowlist: reuse_source contains private host paths.
        public_poses.append({
            "molecule_id": molecule_id, "kind": molecule["kind"],
            "task": molecule["task"],
            "seed_record_id": f'{molecule["task"]}::{molecule["seed_zinc_id"]}',
            "seed_zinc_id": molecule["seed_zinc_id"],
            "pair_id": "__".join(sorted((target, offtarget))),
            "target_uniprot": target, "offtarget_uniprot": offtarget,
            "role": role, "cluster_rank": rank,
            "receptor_uniprot": pose["receptor_uniprot"],
            "endpoint_binding_energy_kcal_mol": energy,
        })
    check_equal(set(energies), set(molecules), "endpoint coverage")

    baselines = {}
    for molecule_id, molecule in molecules.items():
        target_scores = energies[molecule_id]["target"]
        off_scores = energies[molecule_id]["offtarget"]
        if len(target_scores) != 1 or not 1 <= len(off_scores) <= 3:
            raise ValueError(f"Incomplete MM-GBSA poses: {molecule_id}")
        target_stable = flag(molecule["target_pose_stable"])
        off_stable = flag(molecule["offtarget_pose_stable"])
        if off_stable and len(off_scores) != 1:
            raise ValueError("Stable off-target must have one representative")
        baselines[molecule_id] = {
            "record_id": molecule_id, "status": "complete",
            "target_energy_kcal_mol": target_scores[0],
            "offtarget_energy_median_kcal_mol": statistics.median(off_scores),
            "differential_kcal_mol": target_scores[0] - statistics.median(off_scores),
            "target_pose_count": len(target_scores), "offtarget_pose_count": len(off_scores),
            "target_pose_stable": target_stable, "offtarget_pose_stable": off_stable,
        }

    candidates, seeds, used_parents = {}, {}, set()
    for candidate_id, source in master_by_id.items():
        row, molecule = pairs[candidate_id], molecules[candidate_id]
        parent_id = row["parent_prefilter_id"]
        if parent_id not in molecules or molecules[parent_id]["kind"] != "seed":
            raise ValueError(f"Missing actual-parent baseline: {candidate_id}")
        parent = molecules[parent_id]
        for key in ("task", "pair_key", "seed_zinc_id"):
            for linked in (row, molecule, parent):
                check_equal(linked[key], source[key], f"{candidate_id} {key}")
        check_equal(molecule["parent_prefilter_id"], parent_id, "parent ID")
        check_equal(row["kind"], "candidate", "paired kind")
        check_equal(row["status"], "paired_complete", "paired status")
        for status_key in ("candidate_status", "parent_status"):
            check_equal(row[status_key], "complete", status_key)
        for key in ("target_pose_stable", "offtarget_pose_stable"):
            check_equal(flag(row[key]), flag(molecule[key]), f"candidate {key}")
            check_equal(flag(row[f"parent_{key}"]), flag(parent[key]), f"parent {key}")
        candidate_base, parent_base = baselines[candidate_id], baselines[parent_id]
        check_energy(row["candidate_D_median"], candidate_base["differential_kcal_mol"], "candidate D")
        check_energy(row["parent_D_median"], parent_base["differential_kcal_mol"], "parent D")
        delta = candidate_base["differential_kcal_mol"] - parent_base["differential_kcal_mol"]
        check_energy(row["delta_median_kcal_mol"], delta, "candidate-minus-parent delta")
        # Retain the audited source precision after checking reconstructed energies.
        candidates[candidate_id] = {
            **candidate_base, "status": "paired_complete", "parent_id": parent_id,
            "parent_differential_kcal_mol": finite(row["parent_D_median"]),
            "delta_vs_parent_kcal_mol": finite(row["delta_median_kcal_mol"]),
            "improved_vs_parent": delta < 0,
            "parent_target_pose_stable": parent_base["target_pose_stable"],
            "parent_offtarget_pose_stable": parent_base["offtarget_pose_stable"],
        }
        seed_id = f'{source["task"]}::{source["seed_zinc_id"]}'
        if seed_id in seeds:
            check_equal(seeds[seed_id]["record_id"], parent_id, "unique task-specific parent")
        seeds[seed_id] = parent_base
        used_parents.add(parent_id)
    check_equal(used_parents, {key for key, row in molecules.items() if row["kind"] == "seed"},
                "parent coverage")
    public_poses.sort(key=lambda row: (row["molecule_id"], row["role"], row["cluster_rank"]))
    return candidates, seeds, public_poses
