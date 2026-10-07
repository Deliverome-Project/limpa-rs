#!/usr/bin/env bash
# Restore the repo's pinned R packages (renv.lock) into the renv project library.
#
# R packages here are managed by renv -- R's equivalent of uv.lock. The lockfile pins every
# package version (CRAN + Bioconductor) and the R / Bioconductor release they belong to.
# Packages install into renv/library/ inside the repo, not your personal R library, so
# they can't clash with other R work on your machine. Safe to re-run.
#
# Only needed for R-backed helpers (currently deliverome_analysis.limpa, for Spectronaut /
# DIA proteomics). Python-only work never needs this.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOCK="$ROOT/renv.lock"
[[ -f "$LOCK" ]] || { echo "No renv.lock at $LOCK" >&2; exit 1; }

# Bioconductor releases are tied to an R minor version (e.g. Bioc 3.23 <-> R 4.6), so the
# installed R must match the lockfile's major.minor or Bioconductor packages won't resolve.
want="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["R"]["Version"])' "$LOCK")"

find_rscript() {
  RSCRIPT="${DELIVEROME_RSCRIPT:-$(command -v Rscript || true)}"
  if [[ -n "$RSCRIPT" && ! -x "$RSCRIPT" ]]; then
    echo "DELIVEROME_RSCRIPT=$RSCRIPT is not an executable; ignoring it." >&2
    RSCRIPT="$(command -v Rscript || true)"
  fi
  # The CRAN macOS installer links /usr/local/bin/Rscript, which may not be on PATH yet.
  [[ -z "$RSCRIPT" && -x /usr/local/bin/Rscript ]] && RSCRIPT=/usr/local/bin/Rscript
  have=""
  [[ -n "$RSCRIPT" ]] && have="$("$RSCRIPT" --vanilla -e 'cat(paste(R.version$major, R.version$minor, sep="."))')"
}

offer_r_install() {
  # $1: why. Installing R is a system change (needs your admin password), so always ask.
  echo "$1 This repo's R packages need R ${want%.*}.x."
  if [[ "$(uname)" == "Darwin" ]] && command -v brew >/dev/null && [[ -t 0 ]]; then
    read -r -p "Install the current R from CRAN via Homebrew now (asks for your password)? [y/N] " ans
    if [[ "$ans" =~ ^[Yy]$ ]]; then
      brew install --cask r || brew upgrade --cask r
      find_rscript
      return
    fi
  fi
  echo "Install R ${want%.*} yourself, open a new terminal, and re-run this script:" >&2
  echo "  macOS:  brew install --cask r     (or the .pkg from https://cloud.r-project.org)" >&2
  echo "  Linux:  https://cloud.r-project.org/bin/linux/" >&2
  exit 1
}

find_rscript
if [[ -z "$RSCRIPT" ]]; then
  offer_r_install "R is not installed."
elif [[ "${want%.*}" != "${have%.*}" ]]; then
  offer_r_install "Found R $have at $RSCRIPT."
fi
if [[ -z "$RSCRIPT" || "${want%.*}" != "${have%.*}" ]]; then
  echo "Still no R ${want%.*}.x (found: ${have:-none}). CRAN may have moved past ${want%.*};" >&2
  echo "install R ${want%.*} specifically, e.g. with rig: https://github.com/r-lib/rig" >&2
  exit 1
fi
echo "R $have at $RSCRIPT; restoring renv.lock (first run builds packages, ~2-5 min)"

# Running from the repo root lets .Rprofile -> renv/activate.R bootstrap renv itself.
cd "$ROOT"
"$RSCRIPT" -e 'renv::restore(prompt = FALSE)'
"$RSCRIPT" -e 'suppressPackageStartupMessages(library(limpa));
  cat("OK: limpa", format(packageVersion("limpa")), "| limma", format(packageVersion("limma")),
      "| library", .libPaths()[1], "\n")'
