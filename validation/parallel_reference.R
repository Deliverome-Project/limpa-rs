# Repeat scale.R's deterministic workload and compare PR #141's four-worker
# strategy with four Rust threads. Keep global hyperparameters shared.
suppressPackageStartupMessages(library(limpa))
root<-Sys.getenv("LIMPA_RS_ROOT");source(file.path(root,"R/limpa_rs.R"))
Sys.setenv(LIMPA_RS_BIN=file.path(root,"target/release/limpa-rs"))
z<-readRDS(file.path(root,"results/astral.rds"));set.seed(384141)
y<-z$y$E[,rep(seq_len(ncol(z$y$E)),length.out=384),drop=FALSE]
for(j in seq_len(ncol(y))) y[,j]<-y[,j]+rnorm(nrow(y),0,0.1)
colnames(y)<-paste0("expanded_",seq_len(384))
pid<-as.character(z$y$genes$PG.ProteinGroups);o<-order(pid);pid<-pid[o];y<-y[o,]
d<-dpc(y);h<-dpcQuantHyperparam(y,protein.id=pid,dpc.slope=d$dpc[2]);h$prior.sd<-max(h$prior.sd,1);h$prior.logFC<-max(h$prior.logFC,1)
ids<-unique(pid);sel<-unique(round(seq(1,length(ids),length.out=48)));k<-pid%in%ids[sel]
y<-y[k,,drop=FALSE];pid<-pid[k];sigma<-h$sigma[sel];names(sigma)<-ids[sel]
groups<-split(ids[sel],cut(seq_along(sel),4,labels=FALSE))
tref<-system.time(parts<-parallel::mclapply(groups,function(g){
 k<-pid%in%g
 peptides2Proteins(y[k,,drop=FALSE],pid[k],dpc=d$dpc,sigma=sigma[g],prior.mean=h$prior.mean,prior.sd=h$prior.sd,prior.logFC=h$prior.logFC,standard.errors=TRUE)
},mc.cores=4L))[["elapsed"]]
stopifnot(!any(vapply(parts,inherits,logical(1),"try-error")))
refE<-do.call(rbind,lapply(parts,function(p)p$E));refSE<-do.call(rbind,lapply(parts,function(p)p$other$standard.error))
trust<-system.time(out<-limpa_rs_fit(y,pid,d$dpc,sigma,h$prior.mean,h$prior.sd,h$prior.logFC,threads=4L,initialization.groups=groups))[["elapsed"]]
saveRDS(list(y=y,pid=pid,d=d,h=h,sigma=sigma,groups=groups,refE=refE,refSE=refSE),file.path(root,"results/parallel-oracle.rds"))
errE<-max(abs(refE-out$E));errSE<-max(abs(refSE-out$other$standard.error))
metrics<-data.frame(proteins=length(sel),samples=384,r_workers=4,rust_threads=4,reference_parallel_seconds=tref,rust_bridge_seconds=trust,speedup=tref/trust,max_log2_error=errE,max_se_error=errSE,pass=errE<=1e-3&&errSE<=1e-3&&tref/trust>=20)
print(metrics);write.table(metrics,file.path(root,"results/parallel-reference.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
stopifnot(metrics$pass)
