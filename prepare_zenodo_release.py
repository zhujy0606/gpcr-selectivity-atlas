#!/usr/bin/env python3
"""Create Zenodo-ready archives without publishing them."""

from __future__ import annotations

import argparse
import hashlib
import os
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
SURFACES = ROOT / "data/masif/single_gpcr_outputs"
SITE = HERE / "site"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def add_tree(archive: zipfile.ZipFile, source: Path, prefix: str) -> None:
    for path in sorted(
        item for item in source.rglob("*")
        if item.is_file() and item.name != ".DS_Store" and "__MACOSX" not in item.parts
    ):
        archive.write(path, Path(prefix) / path.relative_to(source))


def build(output: Path, version: str) -> list[Path]:
    output.mkdir(parents=True, exist_ok=True)
    table_zip = output / f"gpcr_selectivity_atlas_v{version}_public_tables.zip"
    surface_zip = output / f"gpcr_selectivity_atlas_v{version}_dmasif_assets.zip"
    with zipfile.ZipFile(table_zip, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as archive:
        add_tree(archive, SITE / "data", "data")
        add_tree(archive, SITE / "downloads", "downloads")
        for name in ("README.md", "CITATION.cff", "DATA_DICTIONARY.md", "LICENSE_DATA.md", "PUBLISHING.md"):
            archive.write(HERE / name, name)
    with zipfile.ZipFile(surface_zip, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=4, allowZip64=True) as archive:
        add_tree(archive, SURFACES, "single_gpcr_outputs")
    artifacts = [table_zip, surface_zip]
    sums = output / "SHA256SUMS.txt"
    sums.write_text("".join(f"{sha256(path)}  {path.name}\n" for path in artifacts), encoding="utf-8")
    artifacts.append(sums)
    return artifacts


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", default="0.1.0")
    parser.add_argument("--output", type=Path, default=HERE / "release")
    args = parser.parse_args()
    for artifact in build(args.output, args.version):
        print(f"{artifact.name}\t{os.path.getsize(artifact)} bytes")
