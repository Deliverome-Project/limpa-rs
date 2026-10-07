# SPDX-License-Identifier: GPL-3.0-or-later
# Synthetic only. Save the actual inputs alongside the oracle: RNG changes in
# future R versions must not silently change what is being compared.
reference_cases <- function() {
  set.seed(141)
  cases <- list()
  for (n in c(3L, 6L, 32L, 384L)) {
    counts <- c(1L, 3L, 10L, 27L, 5L, 3L)
    pid <- rep(sprintf("P%04d", seq_along(counts)), counts)
    y <- matrix(rnorm(length(pid)*n, 0, .4), length(pid), n) + rnorm(length(pid), 18, 2)
    y[runif(length(y)) > plogis(-11 + .75*y)] <- NA_real_
    y[pid == "P0002", 1:2] <- NA_real_
    keep <- rowSums(!is.na(y)) > 0
    y <- y[keep,,drop=FALSE]; pid <- pid[keep]
    dimnames(y) <- list(paste0("prec", seq_len(nrow(y))), paste0("run", seq_len(n)))
    cases[[paste0("core-", n)]] <- list(kind="core", y=y, pid=pid,
      dpc=c(-11, .75), sigma=.4, prior.mean=18, prior.sd=3, prior.logFC=2)
  }
  z <- cases[["core-6"]]
  z$y[1,] <- NA_real_
  cases[["all-missing-row"]] <- z
  # Exercise hyperparameters, annotations, targets and singleton dispatch.
  set.seed(142)
  pid <- rep(sprintf("G%03d", 1:60), each=3)
  y <- matrix(rnorm(length(pid)*6, 0, .4), length(pid), 6) + rnorm(length(pid), 18, 2)
  y[runif(length(y)) > plogis(-11+.75*y)] <- NA_real_
  keep <- rowSums(!is.na(y)) > 0
  y <- y[keep,,drop=FALSE]; pid <- pid[keep]
  dimnames(y) <- list(paste0("pep", seq_len(nrow(y))), paste0("run", 1:6))
  cases[["pipeline"]] <- list(kind="pipeline", y=y, pid=pid, dpc=c(-11,.75))
  cases[["singletons"]] <- list(kind="pipeline", y=y[1:6,], pid=paste0("S",1:6), dpc=c(-11,.75))
  cases
}

reference_elist <- function(case) {
  new("EList", list(E=case$y,
    genes=data.frame(protein=case$pid, annotation=paste0("annotation-",case$pid)),
    targets=data.frame(sample=colnames(case$y), group=rep(c("A","B"),length.out=ncol(case$y)))))
}

run_reference_case <- function(case, engine=c("reference", "rust"), threads=1L) {
  engine <- match.arg(engine)
  if (case$kind == "core") {
    args <- list(y=case$y, protein.id=case$pid, dpc=case$dpc, sigma=case$sigma,
      prior.mean=case$prior.mean, prior.sd=case$prior.sd, prior.logFC=case$prior.logFC)
    if (engine == "reference") return(do.call(limpa::peptides2Proteins,
      c(args, list(standard.errors=TRUE))))
    return(do.call(limpa_rs_fit, c(args, list(threads=threads))))
  }
  y <- reference_elist(case)
  if (engine == "reference") return(limpa::dpcQuant(y, "protein", dpc=case$dpc, verbose=FALSE))
  dpc_quant_rust(y, "protein", case$dpc, cores=threads)
}

# Plain lists avoid dependence on a future EList class definition when reading
# an older snapshot. Everything besides E and SE is compared exactly.
reference_output <- function(x) {
  list(E=x$E, se=x$other$standard.error, metadata=list(
    genes=x$genes, targets=x$targets, n.observations=x$other$n.observations,
    dpc=x$dpc, prior.mean=x$prior.mean, prior.sd=x$prior.sd, prior.logFC=x$prior.logFC))
}
