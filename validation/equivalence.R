# Run from a restored PR #141 R environment:
# LIMPA_RS_ROOT=/path/to/limpa-rs Rscript /path/to/limpa-rs/validation/equivalence.R
suppressPackageStartupMessages(library(limpa))
root <- Sys.getenv("LIMPA_RS_ROOT",normalizePath("."))
source(file.path(root,"R/limpa_rs.R"))
Sys.setenv(LIMPA_RS_BIN=file.path(root,"target/release/limpa-rs"))
stopifnot(as.character(packageVersion("limpa"))=="1.4.2")
set.seed(141)
dir.create(file.path(root,"results"),showWarnings=FALSE)
results <- list()
for(n in c(3L,6L,32L,384L)) {
  np<-if(n==384L) 6L else 30L
  counts<-rep(c(1L,3L,10L,27L,50L,5L),length.out=np)
  pid<-rep(sprintf("P%04d",seq_len(np)),counts)
  nr<-length(pid)
  y<-matrix(rnorm(nr*n,0,0.4),nr,n)+rnorm(nr,18,2)
  y[runif(length(y))>plogis(-11+0.75*y)]<-NA
  # Entirely missing samples within a protein; disconnected observations.
  y[pid=="P0002",1:2]<-NA
  y[pid=="P0003",seq_len(max(1,n%/%3))]<-NA
  keep <- rowSums(!is.na(y))>0
  y<-y[keep,,drop=FALSE];pid<-pid[keep];np<-length(unique(pid));sigma<-rep(0.4,np);nr<-nrow(y)
  colnames(y)<-paste0("run",seq_len(n));rownames(y)<-paste0("prec",seq_len(nr))
  d<-c(-11,0.75); sigma<-rep(0.4,np)
  tR<-system.time(ref<-peptides2Proteins(y,pid,sigma=sigma,dpc=d,prior.mean=18,prior.sd=3,prior.logFC=2,standard.errors=TRUE))[["elapsed"]]
  tRust<-system.time(rust<-limpa_rs_fit(y,pid,d,sigma,18,3,2,threads=4L))[["elapsed"]]
  saveRDS(list(y=y,pid=pid,sigma=sigma,d=d,ref=ref,rust=rust),file.path(root,paste0("results/debug-",n,".rds")))
  errE<-max(abs(ref$E-rust$E));errSE<-max(abs(ref$other$standard.error-rust$other$standard.error))
  stopifnot(identical(dimnames(ref$E),dimnames(rust$E)))
  stopifnot(identical(ref$other$n.observations,rust$other$n.observations))
  results[[length(results)+1L]]<-data.frame(samples=n,proteins=np,r_seconds=tR,rust_bridge_seconds=tRust,speedup=tR/tRust,max_log2_error=errE,max_se_error=errSE,pass=errE<=1e-3&&errSE<=1e-3)
  print(results[[length(results)]])
}
result<-do.call(rbind,results)
dir.create(file.path(root,"results"),showWarnings=FALSE)
write.table(result,file.path(root,"results/equivalence.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
if(!all(result$pass)) stop("Reference equivalence tolerance exceeded (1e-3 absolute for log2 and SE)")
