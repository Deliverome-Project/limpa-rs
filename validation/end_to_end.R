suppressPackageStartupMessages({library(limpa);library(data.table)})
root<-Sys.getenv("LIMPA_RS_ROOT");set.seed(141)
dir.create(file.path(root,"results/e2e"),recursive=TRUE,showWarnings=FALSE)
rows<-list();at<-0L
for(p in seq_len(300)) { mu<-rnorm(1,18,2)
 for(k in seq_len(3)) {offset<-rnorm(1)
  for(j in seq_len(6)) {
   x<-mu+offset+ifelse(p<=40&&j>3,2,0)+rnorm(1,0,.3)
   if(runif(1)>plogis(-12+.75*x)) next
   at<-at+1L;rows[[at]]<-data.table(R.FileName=paste0("run",j),R.Condition=ifelse(j<=3,"A","B"),R.Replicate=(j-1)%%3+1,PG.ProteinGroups=sprintf("P%04d",p),PG.ProteinAccessions=sprintf("P%04d",p),EG.ModifiedSequence=paste0("pep",p,"_",k),FG.Charge=2,`EG.TotalQuantity (Settings)`=2^x,EG.Qvalue=.001,PG.Qvalue=.001,EG.IsImputed=FALSE)
  }
 }
}
a<-rbindlist(rows)
a[PG.ProteinGroups=="P0300"&R.FileName=="run1",c("EG.TotalQuantity (Settings)","EG.IsImputed"):=list(2^40,TRUE)]
a[PG.ProteinGroups=="P0299"&R.FileName=="run4",c("EG.TotalQuantity (Settings)","EG.Qvalue"):=list(2^40,.2)]
report<-file.path(root,"results/e2e/report.tsv");fwrite(a,report,sep="\t")
for(engine in c("reference","rust")) {
 args<-c(shQuote(file.path(root,"scripts/limpa_spectronaut.R")),shQuote(paste0("--report=",report)),shQuote(paste0("--outdir=",file.path(root,"results/e2e",engine))),paste0("--engine=",engine),"--seed=141","--cores=1",shQuote("--formula=~ R.Condition"))
 stopifnot(system2(file.path(R.home("bin"),"Rscript"),args)==0)
}
metrics<-list()
for(name in c("protein_log2.tsv","protein_se.tsv","de_R.ConditionB.tsv")) {
 a<-fread(file.path(root,"results/e2e/reference",name));b<-fread(file.path(root,"results/e2e/rust",name))
 setkeyv(a,"protein_id");setkeyv(b,"protein_id");stopifnot(identical(a$protein_id,b$protein_id),identical(names(a),names(b)))
 cols<-names(a)[vapply(a,is.numeric,logical(1))]
 for(k in setdiff(names(a),cols)) stopifnot(identical(a[[k]],b[[k]]))
 err<-vapply(cols,function(k) max(abs(a[[k]]-b[[k]])),numeric(1))
 # B-statistics and t-statistics may amplify very small input changes; retain
 # all column errors in the report, gate scientific outputs separately.
 for(k in cols) metrics[[length(metrics)+1L]]<-data.frame(file=name,column=k,max_absolute_error=err[k])
 if(startsWith(name,"protein")) stopifnot(max(err)<=1e-3)
 else {stopifnot(err["logFC"]<=1e-3,err["P.Value"]<=1e-3,err["adj.P.Val"]<=1e-3,identical(a$adj.P.Val<.05,b$adj.P.Val<.05))}
}
for(name in c("samples.tsv","protein_annotation.tsv","dpc_points.tsv")) stopifnot(identical(fread(file.path(root,"results/e2e/reference",name)),fread(file.path(root,"results/e2e/rust",name))))
b<-fread(file.path(root,"results/e2e/rust/protein_log2.tsv"));stopifnot(b[protein_id=="P0300",run1]<30,b[protein_id=="P0299",run4]<30)
fwrite(rbindlist(metrics),file.path(root,"results/end-to-end.tsv"),sep="\t")
cat("End-to-end outputs, DE calls, q-value and imputation guards agree.\n")
