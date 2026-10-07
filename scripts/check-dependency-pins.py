"""Validate declared Python/R/Rust pins and reject common lock-bypassing installs.

Standard-library-only so CI can run this before installing dependencies.
This is a regression guard, not a sandbox against deliberately obfuscated code.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import sys
import tomllib

ROOT = Path(__file__).resolve().parents[1]
EXACT_PYTHON = re.compile(r"([A-Za-z0-9_.-]+)(?:\[[A-Za-z0-9_,.-]+\])?==([0-9]+(?:\.[0-9]+)*(?:(?:a|b|rc|post|dev)[0-9]+)?)")
EXACT_R = re.compile(r"[0-9]+(?:[.-][0-9]+)+")
BASE_R = {"R", "base", "compiler", "datasets", "graphics", "grDevices", "grid", "methods", "parallel", "splines", "stats", "stats4", "tcltk", "tools", "utils"}
# Vendored upstream bootstrap is the one reviewed installer exception.
BOOTSTRAP_SHA256 = "10cb02fe000afc9814ca287f18c9e525dd8a031712a440e48eef6e1893eed94b"


def normalized(name):
    return re.sub(r"[-_.]+", "-", name).lower()


def check_python(project, lock):
    errors = []
    declared = list(project["project"].get("dependencies", []))
    declared += project["build-system"]["requires"]
    for items in project["project"].get("optional-dependencies", {}).values():
        declared += items
    for items in project.get("dependency-groups", {}).values():
        declared += items
    packages = lock["package"]
    pinned = {(normalized(p["name"]), p["version"]) for p in packages}
    for requirement in declared:
        match = EXACT_PYTHON.fullmatch(requirement) if isinstance(requirement, str) else None
        if not match:
            errors.append(f"Python requirement must use an exact == version: {requirement}")
        elif (normalized(match[1]), match[2]) not in pinned:
            errors.append(f"Python requirement missing from uv.lock: {requirement}")
    for package in packages:
        if package["name"] == "limpa-rs" and package["source"] == {"editable": "."}:
            continue
        if package["source"] != {"registry": "https://pypi.org/simple"}:
            errors.append(f"Unapproved Python source: {package['name']}")
        if not re.fullmatch(r"[0-9][A-Za-z0-9.!+_-]*", package.get("version", "")):
            errors.append(f"Missing/exact version required in uv.lock: {package['name']}")
        artifacts = package.get("wheels", []) + ([package["sdist"]] if "sdist" in package else [])
        if not artifacts or any(not re.fullmatch(r"sha256:[a-f0-9]{64}", a.get("hash", "")) for a in artifacts):
            errors.append(f"Missing Python artifact SHA-256: {package['name']}")
    return errors


def check_r(lock, bootstrap):
    errors = []
    if not EXACT_R.fullmatch(lock["R"].get("Version", "")):
        errors.append("R runtime must have an exact version")
    if not EXACT_R.fullmatch(lock.get("Bioconductor", {}).get("Version", "")):
        errors.append("Bioconductor must have an exact release")
    packages = lock["Packages"]
    for name, package in packages.items():
        if package.get("Package") != name or not EXACT_R.fullmatch(package.get("Version", "")):
            errors.append(f"R package must have an exact version: {name}")
        if package.get("Source") not in {"Repository", "Bioconductor", "GitHub", "Git"}:
            errors.append(f"Unapproved R package source: {name}")
        remote = any(k in package for k in ("RemoteRef", "RemoteSha", "RemoteUrl")) or package.get("Source") in {"GitHub", "Git"}
        if remote and not re.fullmatch(r"[a-f0-9]{40}", package.get("RemoteSha", "")):
            errors.append(f"R remote requires immutable 40-character SHA: {name}")
        if package.get("Source") in {"Repository", "Bioconductor"} and not package.get("Repository"):
            errors.append(f"R package repository missing: {name}")
        for field in ("Depends", "Imports", "LinkingTo"):
            deps = package.get(field, [])
            if isinstance(deps, str):
                deps = deps.split(",")
            for dep in deps:
                dependency = re.split(r"[\s(]", dep.strip())[0]
                if dependency not in BASE_R and dependency not in packages:
                    errors.append(f"Unpinned R transitive dependency: {name} -> {dependency}")
    for name in ("renv", "limpa", "limma", "data.table", "nanoparquet", "iq", "proDA"):
        if name not in packages:
            errors.append(f"Required R dependency absent: {name}")
    version = re.search(r'version <- "([^"]+)"', bootstrap)
    if not version or version[1] != packages.get("renv", {}).get("Version"):
        errors.append("renv bootstrap version differs from renv.lock")
    if hashlib.sha256(bootstrap.encode()).hexdigest() != BOOTSTRAP_SHA256:
        errors.append("renv bootstrap changed: review installer and update its checksum")
    return errors


def check_rust(manifest, lock):
    errors = []
    packages = lock["package"]
    locked = {(p["name"], p["version"]): p for p in packages}
    sections = [manifest, *manifest.get("target", {}).values()]
    if "workspace" in manifest or "patch" in manifest or "replace" in manifest:
        errors.append("Rust workspace/patch/replace requires an explicit pin-policy update")
    for section in sections:
        for kind in ("dependencies", "dev-dependencies", "build-dependencies"):
            for name, spec in section.get(kind, {}).items():
                spec = {"version": spec} if isinstance(spec, str) else spec
                package_name = spec.get("package", name)
                if any(k in spec for k in ("path", "workspace", "registry")):
                    errors.append(f"Unapproved Rust dependency source: {name}")
                elif "git" in spec:
                    rev = spec.get("rev", "")
                    if not re.fullmatch(r"[a-f0-9]{40}", rev) or any(k in spec for k in ("branch", "tag")):
                        errors.append(f"Rust git dependency requires immutable rev: {name}")
                    elif not any(p["name"] == package_name and p.get("source", "").endswith("#" + rev) for p in packages):
                        errors.append(f"Rust git revision absent from lock: {name}")
                else:
                    version = spec.get("version", "")
                    if not re.fullmatch(r"=[0-9]+\.[0-9]+\.[0-9]+(?:-[A-Za-z0-9.-]+)?", version):
                        errors.append(f"Rust dependency needs an exact =version: {name}")
                    elif (package_name, version[1:]) not in locked:
                        errors.append(f"Rust dependency absent from Cargo.lock: {name}")
    for package in packages:
        if package["name"] == manifest["package"]["name"] and "source" not in package:
            continue
        source = package.get("source", "")
        if source == "registry+https://github.com/rust-lang/crates.io-index":
            if not re.fullmatch(r"[a-f0-9]{64}", package.get("checksum", "")):
                errors.append(f"Rust registry checksum missing: {package['name']}")
        elif not (source.startswith("git+https://") and re.search(r"#[a-f0-9]{40}$", source)):
            errors.append(f"Unpinned Rust source in Cargo.lock: {package['name']}")
    return errors


def check_commands(text, name):
    errors = []
    # Join shell continuations so a missing --locked cannot hide across lines.
    text = text.replace("\\\n", " ")
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("#"):
            continue
        if re.search(r"\b(?:pip[0-9.]*|uv\s+pip)\s+install\b", line):
            command = re.sub(r"^-?\s*run:\s*", "", line)
            if command != "uv pip install --no-deps dist/*.whl":
                errors.append(f"{name}: Python installs must use locked sync (only the local wheel may use --no-deps)")
        if re.search(r"\buv\s+(?:sync|run|tool|add)\b", line):
            command = re.sub(r"^-?\s*run:\s*", "", line)
            if command != "uv sync --locked --group ci --no-install-project --no-build":
                errors.append(f"{name}: unsupported Python installation/resolution command")
        if re.search(r"(?:install\.packages|update\.packages|install_(?:github|git|url|version)|(?:renv|pak|BiocManager)::(?:install|pkg_install))\s*\(", line):
            errors.append(f"{name}: R dependencies must use renv::restore from renv.lock")
        if re.search(r"\bcargo\s+(?:build|test|check|clippy|fetch|metadata|vendor|install)\b", line):
            if "--locked" not in line and "--frozen" not in line:
                errors.append(f"{name}: Cargo operations must use --locked or --frozen")
        if re.search(r"python(?:[0-9.]*)?\s+-m\s+build\b", line) and "--no-isolation" not in line:
            errors.append(f"{name}: builds must use the locked environment with --no-isolation")
    return errors


def check_r_imports(text, packages, name):
    # Literal imports and namespace calls; dynamic/obfuscated installation still
    # requires review and is outside this static regression guard's guarantee.
    code = "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("#"))
    names = set(re.findall(r"\b([A-Za-z][A-Za-z0-9.]*)::", code))
    names.update(re.findall(r"\b(?:library|require|requireNamespace)\s*\(\s*['\"]?([A-Za-z][A-Za-z0-9.]*)", code))
    return [f"{name}: R import is absent from renv.lock: {dep}" for dep in sorted(names - set(packages) - BASE_R)]


def main():
    try:
        errors = check_python(tomllib.loads((ROOT / "pyproject.toml").read_text()), tomllib.loads((ROOT / "uv.lock").read_text()))
        errors += check_rust(tomllib.loads((ROOT / "Cargo.toml").read_text()), tomllib.loads((ROOT / "Cargo.lock").read_text()))
        errors += check_r(json.loads((ROOT / "renv.lock").read_text()), (ROOT / "renv/activate.R").read_text())
        for folder in (".github/workflows", "scripts", "python", "R", "validation"):
            for path in (ROOT / folder).rglob("*"):
                if path.is_file() and path.suffix in {".yml", ".yaml", ".sh", ".py", ".R"} and path.resolve() != Path(__file__).resolve():
                    errors += check_commands(path.read_text(), str(path.relative_to(ROOT)))
                    if path.suffix == ".R":
                        errors += check_r_imports(path.read_text(), json.loads((ROOT / "renv.lock").read_text())["Packages"], str(path.relative_to(ROOT)))
    except (OSError, KeyError, TypeError, ValueError) as exc:
        errors = [f"Cannot verify dependency pins: {exc}"]
    for error in errors:
        print(error, file=sys.stderr)
    if not errors:
        print("Python/R/Rust declarations, locks, bootstrap, and install commands are pinned")
    return bool(errors)


if __name__ == "__main__":
    sys.exit(main())
