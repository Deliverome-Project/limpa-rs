# SPDX-License-Identifier: GPL-3.0-or-later
"""Bundle an exact git commit and all locked Rust sources for binary recipients."""

import argparse
import subprocess
import tarfile
import tempfile
import tomllib
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref", required=True, help="tested commit or tag")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    commit = subprocess.check_output(
        ["git", "rev-parse", "--verify", args.ref + "^{commit}"], cwd=root, text=True
    ).strip()
    manifest = subprocess.check_output(
        ["git", "show", f"{commit}:Cargo.toml"], cwd=root, text=True
    )
    version = tomllib.loads(manifest)["package"]["version"]
    if args.output.exists():
        raise FileExistsError(args.output)
    with tempfile.TemporaryDirectory(prefix="limpa-source-") as tmp:
        stage = Path(tmp) / f"limpa-rs-{version}-source"
        stage.mkdir()
        archive = Path(tmp) / "project.tar"
        subprocess.run(
            ["git", "archive", "--format=tar", "--output", str(archive), commit],
            cwd=root,
            check=True,
        )
        subprocess.run(["tar", "-xf", str(archive), "-C", str(stage)], check=True)
        vendor = subprocess.run(
            ["cargo", "vendor", "--locked", "vendor"],
            cwd=stage,
            check=True,
            text=True,
            capture_output=True,
        )
        (stage / ".cargo").mkdir(exist_ok=True)
        (stage / ".cargo/config.toml").write_text(vendor.stdout)
        (stage / "SOURCE_COMMIT").write_text(commit + "\n")
        subprocess.run(
            ["cargo", "metadata", "--offline", "--locked", "--format-version", "1"],
            cwd=stage,
            check=True,
            stdout=subprocess.DEVNULL,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(args.output, "x:gz") as bundle:
            bundle.add(stage, arcname=stage.name)
    print(args.output)


if __name__ == "__main__":
    main()
