"""Enforce a seven-day age for every registry version in the CI Python lock.

Uses only the standard library, before dependency installation. Missing metadata,
missing locks, non-PyPI dependencies, and HTTP failures fail closed. This delay is
one supply-chain defense; it does not detect all malicious packages.
"""
from __future__ import annotations

import json
import sys
import tomllib
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MIN_AGE = timedelta(days=7)
# Exact releases allowed through the seven-day rule because they fix a security issue that
# cannot wait. Remove each entry once it is a week old. certifi 2026.10.19 and pip 26.2.2
# drop the root behind the eMudhra intermediates compromised on 2026-10-19
# (https://blog.cloudflare.com/emudhra-certificate-compromise/).
SECURITY_EXCEPTIONS = {("certifi", "2026.10.19"), ("pip", "26.2.2")}


def registry_packages(text):
    packages = tomllib.loads(text)["package"]
    result = []
    for package in packages:
        source = package["source"]
        if package["name"] == "limpa-rs" and source == {"editable": "."}:
            continue
        if source != {"registry": "https://pypi.org/simple"}:
            raise ValueError(f"Unapproved dependency source: {package['name']}")
        result.append((package["name"], package["version"]))
    if not result:
        raise ValueError("No registry packages in uv.lock")
    return sorted(set(result))


def release_time(metadata):
    # The resolver additionally applies the cooldown to each selected artifact.
    times = [datetime.fromisoformat(f["upload_time_iso_8601"].replace("Z", "+00:00"))
             for f in metadata["urls"]]
    if not times or any(t.tzinfo is None for t in times):
        raise ValueError("Missing or invalid release timestamps")
    return min(times)


def check_package(package, now):
    name, version = package
    request = urllib.request.Request(
        f"https://pypi.org/pypi/{name}/{version}/json",
        headers={"User-Agent": "limpa-rs-dependency-age/1.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            published = release_time(json.load(response))
        if now - published < MIN_AGE and package not in SECURITY_EXCEPTIONS:
            return f"{name}=={version}: release is less than seven days old"
    except Exception as exc:
        return f"{name}=={version}: cannot verify release age ({type(exc).__name__})"
    return None


def main():
    try:
        packages = registry_packages((ROOT / "uv.lock").read_text())
    except (OSError, KeyError, ValueError) as exc:
        print(f"Dependency age check failed: {exc}", file=sys.stderr)
        return 1
    now = datetime.now(timezone.utc)
    with ThreadPoolExecutor(max_workers=4) as pool:
        errors = [e for e in pool.map(lambda p: check_package(p, now), packages) if e]
    for error in errors:
        print(error, file=sys.stderr)
    if errors:
        return 1
    print(f"Verified {len(packages)} Python versions are at least seven days old")
    return 0


if __name__ == "__main__":
    sys.exit(main())
