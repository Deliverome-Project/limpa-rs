# SPDX-License-Identifier: GPL-3.0-or-later
"""Rust-accelerated LIMPA, with the Deliverome analysis Python interface."""

from .api import COMPARATORS, LimpaResult, detection_probability, run_limpa_spectronaut
from .runtime import limpa_available, runtime_info, setup_r, validate

__version__ = "0.1.0rc1"
__all__ = [
    "COMPARATORS",
    "LimpaResult",
    "detection_probability",
    "run_limpa_spectronaut",
    "limpa_available",
    "runtime_info",
    "setup_r",
    "validate",
]
