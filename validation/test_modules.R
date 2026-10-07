# SPDX-License-Identifier: GPL-3.0-or-later
# Base-R tests: do not require LIMPA, Rust, or private data.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
root <- dirname(dirname(normalizePath(script)))
source(file.path(root,"R/limpa_rs.R"))
source(file.path(root,"validation/modules/compare.R"))
expect_error <- function(expr, pattern) {
  error <- tryCatch({force(expr); NULL},error=identity)
  stopifnot(inherits(error,"error"), grepl(pattern,conditionMessage(error)))
}
valid <- list(Version=limpa_rs_reference$version, RemoteSha=limpa_rs_reference$commit)
stopifnot(.limpa_rs_validate_reference(valid))
wrong <- valid; wrong$Version <- "99.0.0"
expect_error(.limpa_rs_validate_reference(wrong), "Unsupported LIMPA")
wrong <- valid; wrong$RemoteSha <- paste(rep("0",40),collapse="")
expect_error(.limpa_rs_validate_reference(wrong), "source commit differs")
# Both entry points must check the pin before touching the input or binary.
limpa_rs_assert_reference <- function() stop("pin enforced")
expect_error(limpa_rs_fit(), "pin enforced")
expect_error(dpc_quant_rust(), "pin enforced")
a <- dget(file.path(root,"validation/fixtures/limpa-1.4.2.R"))
stopifnot(all(compare_reference_snapshots(a,a)$pass))
stopifnot(all(rowSums(!is.na(a$cases$pipeline$y)) > 0))
reject <- function(change) {
  b <- change(a)
  stopifnot(!all(compare_reference_snapshots(a,b)$pass))
}
reject(function(b) {b$outputs[[1]]$E[1,1] <- b$outputs[[1]]$E[1,1]+.002; b})
reject(function(b) {b$outputs[[1]]$se[1,1] <- b$outputs[[1]]$se[1,1]+.002; b})
reject(function(b) {b$outputs[[1]]$E[1,1] <- NA_real_; b})
reject(function(b) {b$outputs[[1]]$se[1,1] <- Inf; b})
reject(function(b) {rownames(b$outputs[[1]]$E)[1] <- "wrong-id"; b})
reject(function(b) {b$outputs[[1]]$E <- b$outputs[[1]]$E[,3:1]; b})
reject(function(b) {b$outputs[[1]]$metadata$n.observations[1,1] <- -1; b})
reject(function(b) {b$outputs[[1]]$metadata$genes$extra <- "changed"; b})
b <- a; b$cases[[1]]$y[1,1] <- 0
expect_error(compare_reference_snapshots(a,b), "Fixture inputs differ")
b <- a; b$outputs <- b$outputs[-1]
expect_error(compare_reference_snapshots(a,b), "Invalid reference snapshot")
b <- a; b$schema <- 2L
expect_error(compare_reference_snapshots(a,b), "Invalid reference snapshot")
# A future version is allowed in the comparison tool, never in production.
b <- a; b$provenance$limpa <- "99.0.0"
stopifnot(all(compare_reference_snapshots(a,b)$pass))
stopifnot(reference_numeric_error(c(NA_real_,Inf),c(NA_real_,-Inf))==Inf)
stopifnot(reference_numeric_error(c(NA_real_),c(NaN))==Inf)
cat("Version guards and equivalence failure-detection tests passed.\n")
