# Expanded REAL precursor workload, not a real 384-sample biological experiment.
# Requires results/astral.rds from real_data.R. Original three columns are repeated
# and perturbed, preserving precursor/protein sizes and empirical missingness.
# Reference timing uses the SAME deterministic 48-protein subset as Rust; no full
# reference runtime is fabricated from an extrapolation.
suppressPackageStartupMessages(library(limpa))
root<-Sys.getenv("LIMPA_RS_ROOT");source(file.path(root,"R/limpa_rs.R"))
Sys.setenv(LIMPA_RS_BIN=file.path(root,"target/release/limpa-rs"))
z<-readRDS(file.path(root,"results/astral.rds"));set.seed(384141)
y<-z$y; y$E<-y$E[,rep(seq_len(ncol(y$E)),length.out=384),drop=FALSE]
# Small independent technical noise, applied one column at a time to bound memory.
for(j in seq_len(ncol(y$E))) y$E[,j]<-y$E[,j]+rnorm(nrow(y$E),0,0.1)
colnames(y$E)<-paste0("expanded_",seq_len(384));y$targets<-NULL
pid<-as.character(y$genes$PG.ProteinGroups);o<-order(pid);pid<-pid[o];y<-y[o,]
ids<-unique(pid)
tprep<-system.time({d<-dpc(y);h<-dpcQuantHyperparam(y,protein.id=pid,dpc.slope=d$dpc[2]);h$prior.sd<-max(h$prior.sd,1);h$prior.logFC<-max(h$prior.logFC,1)})[["elapsed"]]
cat("preparation_seconds=",tprep,"\n");flush.console()
# Even coverage of sorted IDs, identical parameters and starting values on both paths.
sel<-unique(round(seq(1,length(ids),length.out=48)))
k<-pid%in%ids[sel];ys<-y$E[k,,drop=FALSE];ps<-pid[k]
tref<-system.time(ref<-peptides2Proteins(ys,ps,dpc=d$dpc,sigma=h$sigma[sel],prior.mean=h$prior.mean,prior.sd=h$prior.sd,prior.logFC=h$prior.logFC,standard.errors=TRUE))[["elapsed"]]
trustsub<-system.time(rs<-limpa_rs_fit(ys,ps,d$dpc,h$sigma[sel],h$prior.mean,h$prior.sd,h$prior.logFC,threads=4L))[["elapsed"]]
errE<-max(abs(ref$E-rs$E));errSE<-max(abs(ref$other$standard.error-rs$other$standard.error))
cat("subset_reference_seconds=",tref,"subset_rust_seconds=",trustsub,"error_log2=",errE,"error_se=",errSE,"\n");flush.console()
# Full measurement includes global starting-value calculation, IPC, Rust, and readback.
trust<-system.time(full<-limpa_rs_fit(y$E,pid,d$dpc,h$sigma,h$prior.mean,h$prior.sd,h$prior.logFC,threads=4L))[["elapsed"]]
metrics<-data.frame(workload="expanded_three_run_astral_not_real_384_samples",proteins=length(ids),precursors=nrow(y$E),samples=384,threads=4,preparation_seconds=tprep,full_quant_bridge_seconds=trust,total_matrix_to_protein_seconds=tprep+trust,subset_proteins=length(sel),subset_reference_seconds=tref,subset_rust_seconds=trustsub,subset_speedup=tref/trustsub,max_log2_error=errE,max_se_error=errSE,equivalence_pass=errE<=1e-3&&errSE<=1e-3,time_target_pass=tprep+trust<=900,speed_target_subset_pass=tref/trustsub>=20)
print(metrics);write.table(metrics,file.path(root,"results/scale.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
stopifnot(metrics$equivalence_pass,metrics$time_target_pass,metrics$speed_target_subset_pass)
