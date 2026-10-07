# SPDX-License-Identifier: GPL-3.0-or-later
"""Copy canonical R/runtime sources into wheels; never maintain duplicate drivers."""

import os
import subprocess
import sys
from pathlib import Path
from shutil import copy2, copytree

from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.command.bdist_wheel import bdist_wheel


# Use one explicit baseline for both the Rust linker and wheel compatibility tag.
if sys.platform == "darwin":
    deployment = tuple(
        int(x)
        for x in os.environ.get("MACOSX_DEPLOYMENT_TARGET", "11.0").split(".")[:2]
    )
    deployment = max(deployment, (11, 0))
    os.environ["MACOSX_DEPLOYMENT_TARGET"] = ".".join(str(x) for x in deployment)


class BuildRuntime(build_py):
    def run(self):
        super().run()
        root = Path(__file__).parent
        dest = Path(self.build_lib) / "limpa_rs" / "_runtime"
        dest.mkdir(parents=True, exist_ok=True)
        for name in ("R", "validation"):
            copytree(root / name, dest / name, dirs_exist_ok=True)
        (dest / "scripts").mkdir(exist_ok=True)
        copy2(
            root / "scripts/limpa_spectronaut.R", dest / "scripts/limpa_spectronaut.R"
        )
        (dest / "renv").mkdir(exist_ok=True)
        copy2(root / "renv/activate.R", dest / "renv/activate.R")
        copy2(root / "renv.lock", dest / "renv.lock")


class BinaryWheel(bdist_wheel):
    def get_tag(self):
        # The bundled executable is platform-specific, but uses no Python ABI.
        _, _, platform = super().get_tag()
        if sys.platform == "darwin":
            # setup-python's universal Python does not make a single Rust binary
            # universal. Tag the actual compiler target, not Python's architecture.
            target = os.environ.get("CARGO_BUILD_TARGET")
            if not target:
                version = subprocess.check_output(["rustc", "-vV"], text=True)
                target = next(
                    line.split(": ", 1)[1]
                    for line in version.splitlines()
                    if line.startswith("host: ")
                )
            arches = {"aarch64-apple-darwin": "arm64", "x86_64-apple-darwin": "x86_64"}
            if target not in arches:
                raise RuntimeError(f"Unsupported macOS Rust target: {target}")
            baseline = os.environ["MACOSX_DEPLOYMENT_TARGET"].replace(".", "_")
            platform = f"macosx_{baseline}_{arches[target]}"
        return "py3", "none", platform


setup(cmdclass={"build_py": BuildRuntime, "bdist_wheel": BinaryWheel})
