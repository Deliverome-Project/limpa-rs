# SPDX-License-Identifier: GPL-3.0-or-later
# Bridge for limpa 1.4.2. Retains its global DPC/hyperparameter estimates and
# initial values. Only independent per-protein posterior fits move to Rust.
# Algorithm reference: SmythLab/limpa (GPL >= 2); see NOTICE.
# Keep this contract synchronized with renv.lock only after upgrade review.
limpa_rs_reference <- list(
  version = "1.4.2",
  commit = "d9ad5218207434cfcbe080328db5e79d39a1623a"
)

.limpa_rs_validate_reference <- function(description) {
  if (!identical(description$Version, limpa_rs_reference$version))
    stop("Unsupported LIMPA version: expected ", limpa_rs_reference$version,
         "; restore renv.lock. Test upgrades with validation/reference.R first.")
  if (!is.null(description$RemoteSha) &&
      !identical(description$RemoteSha, limpa_rs_reference$commit))
    stop("LIMPA source commit differs from the assessed reference; restore renv.lock")
  invisible(TRUE)
}

limpa_rs_assert_reference <- function() {
  .limpa_rs_validate_reference(utils::packageDescription("limpa"))
}

limpa_rs_fit <- function(y, protein.id, dpc, sigma, prior.mean, prior.sd,
                         prior.logFC, binary = Sys.getenv("LIMPA_RS_BIN"), threads = 1L, initialization.groups = NULL) {
  limpa_rs_assert_reference()
  stopifnot(nzchar(binary), file.exists(binary), threads >= 1L)
  y <- as.matrix(y)
  stopifnot(nrow(y) == length(protein.id), ncol(y) >= 2L, !anyNA(protein.id))
  o <- order(protein.id); protein.id <- as.character(protein.id[o]); y <- y[o,,drop=FALSE]
  ids <- unique(protein.id); np <- length(ids); n <- ncol(y)
  starts <- which(!duplicated(protein.id)); ends <- c(starts[-1]-1L,nrow(y))
  if (length(sigma)==1L) sigma <- rep(sigma,np)
  stopifnot(length(sigma)==np, all(is.finite(sigma)), all(sigma>0))
  # Completely unobserved precursors have weakly identified peptide effects.
  # Preserve the reference's weighted-row collapse and optimizer behavior exactly.
  if(any(rowSums(!is.na(y))==0L)) {
    warning("All-missing precursor rows: using the reference LIMPA fitter")
    return(limpa::peptides2Proteins(y,protein.id,sigma=sigma,dpc=dpc,
      prior.mean=prior.mean,prior.sd=prior.sd,prior.logFC=prior.logFC,standard.errors=TRUE))
  }
  if(is.null(initialization.groups)) {
    yi <- limpa::imputeByExpTilt(y,dpc.slope=dpc[2],prior.logfc=prior.logFC)
  } else {
    flat <- unlist(initialization.groups,use.names=FALSE)
    stopifnot(identical(sort(as.character(flat)),sort(ids)))
    yi <- y
    for(g in initialization.groups) {
      k <- protein.id %in% g
      yi[k,] <- limpa::imputeByExpTilt(y[k,,drop=FALSE],dpc.slope=dpc[2],prior.logfc=prior.logFC)
    }
  }
  infile <- tempfile(fileext=".bin"); outfile <- tempfile(fileext=".bin")
  on.exit(unlink(c(infile,outfile)),add=TRUE)
  con <- file(infile,"wb")
  tryCatch({
    writeBin(charToRaw("LIMPAR01"),con)
    writeBin(as.integer(c(np,n)),con,size=4,endian="little")
    writeBin(as.double(c(dpc,prior.mean,prior.sd,prior.logFC)),con,size=8,endian="little")
    for (i in seq_len(np)) {
      rows <- starts[i]:ends[i]; z <- yi[rows,,drop=FALSE]
      b <- rowMeans(z); b <- b-mean(b)
      beta <- colMeans(z)
      missing <- colSums(!is.na(y[rows,,drop=FALSE]))==0
      if (sum(missing)>1L) beta[missing] <- mean(beta[missing])
      if(length(rows)>1L) beta <- c(beta,b[-length(b)])
      writeBin(as.integer(length(rows)),con,size=4,endian="little")
      writeBin(as.double(sigma[i]),con,size=8,endian="little")
      writeBin(as.double(t(y[rows,,drop=FALSE])),con,size=8,endian="little")
      writeBin(as.double(beta),con,size=8,endian="little")
    }
  },finally=close(con))
  status <- system2(binary,c(shQuote(infile),shQuote(outfile),as.integer(threads)))
  if(status!=0L) stop("Rust quantification failed; no result accepted")
  con <- file(outfile,"rb")
  z <- tryCatch({
    stopifnot(rawToChar(readBin(con,"raw",8))=="LIMPAO01")
    stopifnot(identical(readBin(con,"integer",2,size=4,endian="little"),as.integer(c(np,n))))
    values <- readBin(con,"double",np*n*2,size=8,endian="little")
    stopifnot(length(values)==np*n*2,all(is.finite(values)),length(readBin(con,"raw",1))==0L)
    matrix(values,nrow=np,byrow=TRUE)
  },finally=close(con))
  E <- z[,seq_len(n),drop=FALSE]; se <- z[,n+seq_len(n),drop=FALSE]
  dimnames(E)<-dimnames(se)<-list(ids,colnames(y))
  nobs <- t(vapply(seq_len(np),function(i) colSums(!is.na(y[starts[i]:ends[i],,drop=FALSE])),numeric(n)))
  dimnames(nobs)<-dimnames(E)
  new("EList",list(E=E,genes=data.frame(NPrec=as.double(ends-starts+1L),PropObs=rowMeans(nobs)/(ends-starts+1L),row.names=ids),
                  other=list(n.observations=nobs,standard.error=se)))
}

dpc_quant_rust <- function(y, protein_id, dpcfit, cores=1L, reference.cores=1L) {
  limpa_rs_assert_reference()
  dpcv <- if(is.list(dpcfit)) dpcfit$dpc else dpcfit
  pid <- as.character(y$genes[[protein_id]])
  if(!length(pid)||anyNA(pid)) stop("invalid protein IDs")
  # limpa has a distinct all-singleton estimator; retain it exactly.
  if(!anyDuplicated(pid)) return(limpa::dpcQuant(y,protein_id,dpc=dpcfit,verbose=FALSE))
  y$genes[[protein_id]] <- NULL
  o <- order(pid);pid<-pid[o];y<-y[o,]
  h <- limpa::dpcQuantHyperparam(y,protein.id=pid,dpc.slope=dpcv[2])
  h$prior.sd <- max(h$prior.sd,1);h$prior.logFC<-max(h$prior.logFC,1)
  stopifnot(length(reference.cores)==1L,is.finite(reference.cores),reference.cores>=1L)
  ids <- unique(pid)
  groups <- if(reference.cores>1L && length(ids)>1L)
    split(ids,cut(seq_along(ids),min(reference.cores,length(ids)),labels=FALSE)) else NULL
  out <- limpa_rs_fit(y,pid,dpcv,h$sigma,h$prior.mean,h$prior.sd,h$prior.logFC,
                       threads=cores,initialization.groups=groups)
  if(!is.null(y$genes)&&ncol(y$genes)) {
    genes <- y$genes[,!colnames(y$genes)%in%colnames(out$genes),drop=FALSE]
    dup <- duplicated(pid)
    keep <- vapply(seq_len(ncol(genes)),function(i) all(duplicated(genes[,i])[dup]),logical(1))
    genes <- genes[!dup,keep,drop=FALSE];rownames(genes)<-rownames(out$E)
    out$genes<-data.frame(genes,out$genes)
  }
  out$targets<-y$targets;out$dpc<-dpcv
  out$prior.mean<-h$prior.mean;out$prior.sd<-h$prior.sd;out$prior.logFC<-h$prior.logFC
  out
}
