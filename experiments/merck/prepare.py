#!/usr/bin/env python3
"""Download and prepare the pinned BatchIE Merck retrospective screen."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import urllib.request
import venv

DATA_URL = (
    "https://zenodo.org/records/12764821/files/"
    "merck_2016.screen.h5?download=1"
)
DATA_MD5 = "5d7cb090002a01ebcc0c06f0d5f72f67"
SOURCE_URL = (
    "https://zenodo.org/api/records/12765294/files/"
    "tansey-lab/batchie-v0.0.1-zenodo.zip/content"
)
SOURCE_MD5 = "7deaf0ab4d9fb4c935e857d1937762bb"
SOURCE_COMMIT = "6baa258cb5d77430dccab38b73b0908f91503c91"


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, path: Path, checksum: str) -> None:
    if path.exists() and md5(path) == checksum:
        print(f"Reusing verified {path}")
        return
    temporary = path.with_suffix(path.suffix + ".part")
    temporary.unlink(missing_ok=True)
    print(f"Downloading {url}")
    urllib.request.urlretrieve(url, temporary)
    actual = md5(temporary)
    if actual != checksum:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"Checksum mismatch for {path}: {actual} != {checksum}")
    temporary.replace(path)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, default=Path(__file__).parent / "work")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--force-prepare", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if sys.version_info < (3, 11):
        raise RuntimeError("Preparation requires Python 3.11 or newer")

    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    data = work / "merck_2016.screen.h5"
    source = work / "batchie-v0.0.1-zenodo.zip"
    training = work / "merck_training.screen.h5"
    holdout = work / "merck_holdout.screen.h5"
    manifest_path = work / "preparation.json"
    output_rds = work / "merck_batchieR.rds"
    download(DATA_URL, data, DATA_MD5)
    download(SOURCE_URL, source, SOURCE_MD5)

    environment = work / ".venv"
    python = environment / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    executable = environment / (
        "Scripts/prepare_retrospective_simulation.exe"
        if os.name == "nt"
        else "bin/prepare_retrospective_simulation"
    )
    if not executable.exists():
        print(f"Creating upstream environment in {environment}")
        venv.EnvBuilder(with_pip=True).create(environment)
        subprocess.run(
            [str(python), "-m", "pip", "install", str(source)],
            check=True,
        )

    if args.force_prepare or not (training.exists() and holdout.exists()):
        command = [
            str(executable),
            "--data", str(data),
            "--training-output", str(training),
            "--test-output", str(holdout),
            "--plate-generator", "PairwisePlateGenerator",
            "--plate-generator-param", "subset_size=20",
            "--plate-generator-param", "anchor_size=0",
            "--initial-plate-generator", "SparseCoverPlateGenerator",
            "--initial-plate-generator-param", "reveal_single_treatment_experiments=False",
            "--plate-smoother", "BatchieEnsemblePlateSmoother",
            "--plate-smoother-param", "min_size=50",
            "--plate-smoother-param", "n_iterations=1",
            "--plate-smoother-param", "min_n_cell_line_plates=3",
            "--holdout-fraction", "0.1",
            "--seed", "0",
        ]
        print("Running pinned upstream retrospective preparation")
        subprocess.run(command, check=True)

    manifest = {
        "created_at": __import__("datetime").datetime.now(
            __import__("datetime").timezone.utc
        ).isoformat(),
        "data_url": DATA_URL,
        "data_md5": DATA_MD5,
        "upstream_archive_url": SOURCE_URL,
        "upstream_archive_md5": SOURCE_MD5,
        "upstream_commit": SOURCE_COMMIT,
        "preparation_seed": 0,
        "holdout_fraction": 0.1,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    converter = Path(__file__).parent / "convert.R"
    subprocess.run(
        [
            args.rscript,
            "--vanilla",
            str(converter),
            "--training", str(training),
            "--holdout", str(holdout),
            "--manifest", str(manifest_path),
            "--output", str(output_rds),
        ],
        check=True,
    )
    print(f"Prepared batchieR data: {output_rds}")


if __name__ == "__main__":
    main()
