# SPDX-License-Identifier: GPL-3.0-or-later
"""Explicit, user-writable R runtime setup; importing performs no installation."""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
from pathlib import Path

ASSETS = Path(__file__).resolve().parent / "_runtime"
INSTALL_HINT = "Run `python -m limpa_rs setup-r` to restore the pinned R environment."


def _rscript() -> str:
    candidate = os.environ.get("LIMPA_RS_RSCRIPT") or os.environ.get(
        "DELIVEROME_RSCRIPT"
    )
    found = shutil.which(candidate or "Rscript")
    if not found and not candidate and Path("/usr/local/bin/Rscript").is_file():
        found = "/usr/local/bin/Rscript"
    if not found:
        raise RuntimeError("Rscript not found; install R 4.6.x. " + INSTALL_HINT)
    return found


def _binary() -> Path:
    override = os.environ.get("LIMPA_RS_BIN")
    binary = (
        Path(override).expanduser().resolve()
        if override
        else Path(__file__).parent / "_native"
    )
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise RuntimeError(
            "Packaged Rust executable missing or not executable; reinstall limpa-rs"
        )
    return binary


def _project(project=None) -> Path:
    configured = project or os.environ.get("LIMPA_RS_R_PROJECT")
    if configured:
        return Path(configured).expanduser().resolve()
    digest = hashlib.sha256((ASSETS / "renv.lock").read_bytes()).hexdigest()[:16]
    cache = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
    return (cache / "limpa-rs" / digest).resolve()


def _check_project(project: Path) -> None:
    if (
        not (project / "renv.lock").is_file()
        or not (project / "renv/activate.R").is_file()
    ):
        raise RuntimeError("Pinned R project not configured. " + INSTALL_HINT)
    if (project / "renv.lock").read_bytes() != (ASSETS / "renv.lock").read_bytes():
        raise RuntimeError(
            "R project lock differs from the assessed environment. " + INSTALL_HINT
        )


def _environment(project: Path) -> dict[str, str]:
    env = os.environ.copy()
    env.update(
        RENV_PROJECT=str(project),
        RENV_CONFIG_AUTOLOADER_ENABLED="TRUE",
        RENV_CONFIG_PROMPT_ENABLED="FALSE",
        LIMPA_RS_ASSETS=str(ASSETS),
        R_PROFILE_USER=str(Path(__file__).parent / "runtime-profile.R"),
        LIMPA_RS_BIN=str(_binary()),
    )
    # --no-environ ignores caller .Renviron; use an explicit profile instead of
    # executing whichever .Rprofile happens to be in an analysis directory.
    return env


def _command(project=None):
    project = _project(project)
    _check_project(project)
    return (
        [_rscript(), "--no-environ", "--no-site-file"],
        project,
        _environment(project),
    )


def _run(args, *, project=None, verbose=False):
    command, cwd, env = _command(project)
    proc = subprocess.run(
        command + args, cwd=cwd, env=env, capture_output=True, text=True
    )
    if verbose:
        print(proc.stdout, end="")
        print(proc.stderr, end="")
    if proc.returncode:
        raise RuntimeError(
            f"LIMPA failed (exit {proc.returncode}):\n{proc.stderr.strip()}\n{INSTALL_HINT}"
        )
    return proc


def setup_r(project=None) -> Path:
    """Restore bundled renv.lock. Explicit network/install operation; safe to repeat.

    A custom project must be selected again through LIMPA_RS_R_PROJECT for analyses.
    A pre-existing project is accepted only when its entire lock matches.
    """
    rs = _rscript()
    version = subprocess.run(
        [
            rs,
            "--vanilla",
            "-e",
            'cat(paste(R.version$major, strsplit(R.version$minor,".",fixed=TRUE)[[1]][1],sep="."))',
        ],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    if version != "4.6":
        raise RuntimeError(
            f"Expected R 4.6.x, found {version}; install the pinned R minor release"
        )
    project = _project(project)
    project.mkdir(parents=True, exist_ok=True)
    lock = project / "renv.lock"
    if lock.exists() and lock.read_bytes() != (ASSETS / "renv.lock").read_bytes():
        raise RuntimeError("Refusing to overwrite a different R project lock")
    if not lock.exists():
        shutil.copy2(ASSETS / "renv.lock", lock)
    activate = project / "renv/activate.R"
    activate.parent.mkdir(exist_ok=True)
    if not activate.exists():
        shutil.copy2(ASSETS / "renv/activate.R", activate)
    # Installation output is streamed; no silent package installation during import/run.
    command, cwd, env = _command(project)
    proc = subprocess.run(
        command
        + [
            "-e",
            'renv::restore(prompt=FALSE); source(file.path(Sys.getenv("LIMPA_RS_ASSETS"),"R/limpa_rs.R")); limpa_rs_assert_reference()',
        ],
        cwd=cwd,
        env=env,
    )
    if proc.returncode:
        raise RuntimeError("Pinned R restore failed; see installation output above")
    return project


def runtime_info() -> dict:
    """Verify the installed binary, R and LIMPA pin; raise with setup diagnostics."""
    proc = _run(
        [
            "-e",
            'source(file.path(Sys.getenv("LIMPA_RS_ASSETS"),"R/limpa_rs.R")); limpa_rs_assert_reference(); cat(as.character(packageVersion("limpa")))',
        ]
    )
    return {
        "limpa_version": proc.stdout.strip().splitlines()[-1],
        "rscript": _rscript(),
        "r_project": str(_project()),
        "binary": str(_binary()),
    }


def limpa_available() -> bool:
    """True only when the bundled Rust executable and pinned R reference are usable."""
    try:
        runtime_info()
        return True
    except (RuntimeError, OSError, subprocess.SubprocessError):
        return False


def validate() -> None:
    """Check the installed binary against the frozen synthetic R oracle; raises on drift."""
    for script in ("test_modules.R", "reference.R"):
        _run([str(ASSETS / "validation" / script)], verbose=True)
