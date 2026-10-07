# SPDX-License-Identifier: GPL-3.0-or-later
"""Copy canonical R/runtime sources into wheels; never maintain duplicate drivers."""

from pathlib import Path
from shutil import copy2, copytree

from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.command.bdist_wheel import bdist_wheel


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
        return "py3", "none", platform


setup(cmdclass={"build_py": BuildRuntime, "bdist_wheel": BinaryWheel})
