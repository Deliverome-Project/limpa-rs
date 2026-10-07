# SPDX-License-Identifier: GPL-3.0-or-later
# check: pinned R and Rust against the committed oracle (default)
# capture VERSION OUTPUT.rds: evaluate installed candidate R against SAME inputs
# compare BASELINE.rds CANDIDATE.rds: compare saved runs, no LIMPA required
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
root <- dirname(dirname(normalizePath(script)))
source(file.path(root, "R/limpa_rs.R"))
source(file.path(root, "validation/modules/fixtures.R"))
source(file.path(root, "validation/modules/compare.R"))
args <- commandArgs(trailingOnly=TRUE)
mode <- if (length(args)) args[1] else "check"
read_snapshot <- function(path) if (grepl("[.]rds$",path)) readRDS(path) else dget(path)
oracle <- file.path(root,"validation/fixtures/limpa-1.4.2.R")
provenance <- function(engine) list(limpa=as.character(packageVersion("limpa")),
  limpa_commit=utils::packageDescription("limpa")$RemoteSha,
  R=as.character(getRversion()), limma=as.character(packageVersion("limma")),
  statmod=as.character(packageVersion("statmod")),
  data.table=as.character(packageVersion("data.table")),
  binary_md5=if (engine == "rust") unname(tools::md5sum(Sys.getenv("LIMPA_RS_BIN"))) else NULL,
  engine=engine, seed=141L)
evaluate <- function(cases, engine, threads=1L) {
  list(schema=1L, provenance=provenance(engine), cases=cases,
    outputs=lapply(cases, function(case) reference_output(run_reference_case(case,engine,threads))))
}
report <- function(x) {
  print(x, row.names=FALSE)
  if (!all(x$pass)) stop("Equivalence failed; do not update the production pin")
}
if (mode == "compare") {
  if (length(args)!=3L) stop("Usage: reference.R compare BASELINE CANDIDATE")
  report(compare_reference_snapshots(read_snapshot(args[2]),read_snapshot(args[3])))
} else {
  suppressPackageStartupMessages(library(limpa))
  if (mode == "capture") {
    if (length(args)!=3L) stop("Usage: reference.R capture EXPECTED_VERSION OUTPUT.rds")
    if (!identical(as.character(packageVersion("limpa")),args[2])) stop("Installed candidate version differs")
    if (file.exists(args[3])) stop("Refusing to overwrite snapshot")
    baseline <- read_snapshot(oracle)
    validate_reference_snapshot(baseline)
    snapshot <- evaluate(baseline$cases,"reference")
    saveRDS(snapshot,args[3],version=3)
    cat("Saved candidate reference; production pin unchanged.\n")
    report(compare_reference_snapshots(baseline,snapshot))
  } else if (mode == "check") {
    if (length(args)>1L) stop("Usage: reference.R check")
    limpa_rs_assert_reference()
    if (!nzchar(Sys.getenv("LIMPA_RS_BIN")))
      Sys.setenv(LIMPA_RS_BIN=file.path(root,"target/release/limpa-rs"))
    baseline <- read_snapshot(oracle)
    if (!identical(baseline$provenance$limpa,limpa_rs_reference$version) ||
        !identical(baseline$provenance$limpa_commit,limpa_rs_reference$commit))
      stop("Oracle provenance differs from production pin")
    report(compare_reference_snapshots(baseline,evaluate(baseline$cases,"reference")))
    for (threads in c(1L,4L)) {
      cat("Rust threads:",threads,"\n")
      report(compare_reference_snapshots(baseline,evaluate(baseline$cases,"rust",threads)))
    }
  } else stop("Unknown mode; use check, capture or compare")
}
