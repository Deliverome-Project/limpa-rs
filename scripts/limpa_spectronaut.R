# SPDX-License-Identifier: GPL-3.0-or-later
# Adapted from deliverome-analysis PR #141, commit 32c563c95410637e0431faf27eab1e7a9786464b.
# limpa quantification (+ optional differential expression) of a Spectronaut report.
#
# Called by deliverome_analysis.limpa.run_limpa_spectronaut(); not meant to be run by hand,
# but it works standalone:
#
#   Rscript limpa_spectronaut.R --report=Report.tsv --outdir=out \
#       (or --matrix-dir=<deliverome_analysis.spectronaut_trim output> instead of --report)
#       [--cores=N]  per-protein fits in parallel (forked; matches serial to optimizer tolerance) \
#       [--protein-id=PG.ProteinGroups] [--q-cutoff=0.01] \
#       [--samples=samples.tsv] [--formula="~ R.Condition"] [--contrasts="B-A,C-A"] \
#       [--compare=maxlfq_proda]
#
# Input: a Spectronaut *Normal Report* (long format, one row per precursor per run).
# limpa::readSpectronaut keeps EG.Qvalue and PG.Qvalue <= q-cutoff and drops any value with
# EG.IsImputed == TRUE, so Spectronaut imputation can never leak into the quantification.
#
# Outputs (all in --outdir):
#   protein_log2.tsv        protein x run, DPC-Quant log2 expression (complete, no NAs)
#   protein_se.tsv          protein x run, standard error of each value
#   protein_annotation.tsv  protein id, annotation columns, NPrec, PropObs
#   samples.tsv             run-level table (Spectronaut run info + any --samples columns)
#   dpc_points.tsv          precursors used to fit the DPC: limpa's mean estimate, n detected
#   summary.json            DPC coefficients, counts, missingness, package versions
#   de_<name>.tsv           one limma topTable per coefficient / contrast (only with --formula)
#   maxlfq_log2.tsv         protein x run MaxLFQ (iq::fast_MaxLFQ), NAs kept (with --compare)
#   cmp_<method>__<name>.tsv  comparator DE tables, columns as limma topTable (with --compare)
#
# The comparator follows the limpa paper (Li, Cobbold & Smyth 2025), run on the same filtered
# precursor matrix limpa uses:
#   maxlfq_proda    MaxLFQ, then proDA (probabilistic dropout model)

suppressPackageStartupMessages({
  library(limpa)
  library(data.table)
})

parse_args <- function(argv) {
  out <- list()
  for (a in argv) {
    m <- regmatches(a, regexec("^--([^=]+)=(.*)$", a))[[1]]
    if (length(m) != 3) stop("Bad argument (expected --key=value): ", a)
    out[[m[2]]] <- m[3]
  }
  out
}

json_escape <- function(x) gsub('"', '\\\\"', gsub("\\\\", "\\\\\\\\", x))

to_json <- function(x) {
  # Minimal JSON writer for a flat named list of scalars / short vectors (avoids a jsonlite dep).
  fmt <- function(v) {
    if (is.null(v) || length(v) == 0) return("null")
    if (is.character(v)) {
      s <- paste0('"', json_escape(v), '"')
    } else if (is.logical(v)) {
      s <- tolower(as.character(v))
    } else {
      s <- ifelse(is.finite(v), format(v, digits = 10, scientific = FALSE, trim = TRUE), "null")
    }
    if (length(v) == 1 && is.null(names(v))) s else {
      if (!is.null(names(v))) paste0("{", paste0('"', names(v), '": ', s, collapse = ", "), "}")
      else paste0("[", paste(s, collapse = ", "), "]")
    }
  }
  paste0("{\n", paste0('  "', names(x), '": ', vapply(x, fmt, ""), collapse = ",\n"), "\n}\n")
}

write_matrix <- function(m, id_name, path) {
  d <- data.table(rownames(m), m)
  setnames(d, 1, id_name)
  fwrite(d, path, sep = "\t", na = "NA")
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
report <- args[["report"]]
matrix_dir <- args[["matrix-dir"]]
outdir <- args[["outdir"]]
if (is.null(outdir) || (is.null(report) == is.null(matrix_dir))) {
  stop("--outdir and exactly one of --report / --matrix-dir are required")
}
cores <- if (is.null(args[["cores"]])) 1L else as.integer(args[["cores"]])
protein_id <- if (is.null(args[["protein-id"]])) "PG.ProteinGroups" else args[["protein-id"]]
q_cutoff <- if (is.null(args[["q-cutoff"]])) 0.01 else as.numeric(args[["q-cutoff"]])
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- read ---------------------------------------------------------------------------------
read_report <- function() {
  # Spectronaut reports differ in which PG annotation columns they carry (the BGS Factory
  # report has no PG.Genes), so ask only for the ones that exist. readSpectronaut errors on a
  # missing column rather than skipping it.
  header <- names(fread(report, nrows = 0, sep = "\t", check.names = FALSE))
  required <- c("R.FileName", "EG.ModifiedSequence", "FG.Charge",
                "EG.TotalQuantity (Settings)", "EG.Qvalue", "PG.Qvalue", protein_id)
  missing_cols <- setdiff(required, header)
  if (length(missing_cols)) {
    stop("Report is missing required column(s): ", paste(missing_cols, collapse = ", "),
         ". Export a Spectronaut *Normal Report* that includes them.")
  }
  has_flag <- "EG.IsImputed" %in% header
  if (!has_flag) {
    message("WARNING: report has no EG.IsImputed column, so Spectronaut-imputed values cannot be ",
            "removed. Re-export with that column, or set imputation to 'No Imputing' in Spectronaut.")
  }
  annotation <- unique(intersect(c(protein_id, "PG.ProteinAccessions", "PG.Genes"), header))
  y <- readSpectronaut(
    report,
    annotation.columns = annotation,
    q.cutoffs = q_cutoff,
    filter.columns = if (has_flag) "EG.IsImputed" else NULL,
    filter.values = if (has_flag) TRUE else NULL,
    verbose = TRUE
  )
  list(y = y, has_flag = has_flag, source = normalizePath(report))
}

read_matrix_dir <- function() {
  # Output of deliverome_analysis.spectronaut_trim: already filtered exactly as readSpectronaut
  # filters, log2, NaN = not detected.
  man <- jsonlite_free_read(file.path(matrix_dir, "manifest.json"))
  wide <- as.data.frame(nanoparquet::read_parquet(file.path(matrix_dir, "precursor_log2.parquet")))
  ann <- as.data.frame(nanoparquet::read_parquet(file.path(matrix_dir, "precursors.parquet")))
  E <- as.matrix(wide[, -1, drop = FALSE])
  E[is.nan(E)] <- NA
  rownames(E) <- wide$precursor
  ann <- ann[match(wide$precursor, ann$precursor), , drop = FALSE]
  names(ann)[names(ann) == "protein_id"] <- protein_id
  genes <- ann[, setdiff(names(ann), "precursor"), drop = FALSE]
  rownames(genes) <- wide$precursor
  smp <- as.data.frame(fread(file.path(matrix_dir, "samples.tsv"), sep = "\t"))
  rownames(smp) <- smp$run
  tg <- smp[colnames(E), setdiff(names(smp), "run"), drop = FALSE]
  y <- new("EList", list(E = E, genes = genes, targets = tg))
  list(y = y, has_flag = isTRUE(man$imputed_flag_present),
       source = paste0(normalizePath(matrix_dir), " (from ", man$source_report, ")"))
}

# Minimal reader for the trimmer's flat manifest (avoids a jsonlite dependency).
jsonlite_free_read <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "")
  out <- list()
  for (m in regmatches(txt, gregexpr('"[^"]+"\\s*:\\s*("[^"]*"|true|false|[-0-9.eE]+)', txt))[[1]]) {
    kv <- regmatches(m, regexec('"([^"]+)"\\s*:\\s*(.*)$', m))[[1]]
    v <- gsub('^"|"$', "", kv[3])
    out[[kv[2]]] <- if (v %in% c("true", "false")) v == "true" else v
  }
  out
}

input <- if (!is.null(report)) read_report() else read_matrix_dir()
y <- input$y
has_imputed_flag <- input$has_flag
if (ncol(y$E) < 2) stop("limpa needs at least 2 runs; this report has ", ncol(y$E))

# ---- detection probability curve + protein quantification ----------------------------------
# Parallel version of limpa::dpcQuant.EList: the Bayes hyperparameters are estimated once from
# all proteins, then each protein's fit is independent given them, so proteins are split across
# forked workers (same order, same annotation). Not bit-identical to serial: peptides2Proteins
# derives optimizer starting values from a simple imputation over its whole input, which differs
# per chunk, so estimates agree to optimizer tolerance (max 4e-4 log2 on the liver export).
dpc_quant_parallel <- function(y, protein_id, dpcfit, cores) {
  if (cores <= 1L) return(dpcQuant(y, protein_id, dpc = dpcfit, verbose = FALSE))
  dpcv <- dpcfit$dpc
  pid <- as.character(y$genes[[protein_id]])
  y$genes[[protein_id]] <- NULL
  if (!anyDuplicated(pid)) return(dpcQuant(y, pid, dpc = dpcfit, verbose = FALSE))
  o <- order(pid); pid <- pid[o]; y <- y[o, ]
  h <- dpcQuantHyperparam(y, protein.id = pid, dpc.slope = dpcv[2])
  h$prior.sd <- max(h$prior.sd, 1); h$prior.logFC <- max(h$prior.logFC, 1)
  prots <- unique(pid)
  grp <- split(prots, cut(seq_along(prots), min(cores, length(prots)), labels = FALSE))
  parts <- parallel::mclapply(grp, function(g) {
    k <- pid %in% g
    peptides2Proteins(y[k, ], protein.id = pid[k], dpc = dpcv, sigma = h$sigma[g],
                      prior.mean = h$prior.mean, prior.sd = h$prior.sd,
                      prior.logFC = h$prior.logFC, standard.errors = TRUE)
  }, mc.cores = cores)
  if (any(vapply(parts, inherits, logical(1), "try-error"))) stop("a parallel worker failed")
  yp <- new("EList", list(
    E = do.call(rbind, lapply(parts, function(p) p$E)),
    genes = do.call(rbind, lapply(parts, function(p) p$genes)),
    other = list(
      n.observations = do.call(rbind, lapply(parts, function(p) p$other$n.observations)),
      standard.error = do.call(rbind, lapply(parts, function(p) p$other$standard.error)))))
  d <- !duplicated(pid)
  if (!is.null(y$genes) && ncol(y$genes)) {
    genes <- y$genes[, !(colnames(y$genes) %in% colnames(yp$genes)), drop = FALSE]
    dup <- duplicated(pid)
    keep <- vapply(seq_len(ncol(genes)), function(i) all(duplicated(genes[, i])[dup]), logical(1))
    genes <- genes[d, keep, drop = FALSE]
    rownames(genes) <- rownames(yp$E)
    yp$genes <- data.frame(genes, yp$genes)
  }
  yp$targets <- y$targets
  yp$dpc <- dpcv
  yp$prior.mean <- h$prior.mean; yp$prior.sd <- h$prior.sd; yp$prior.logFC <- h$prior.logFC
  yp
}

engine <- if(is.null(args[["engine"]])) "rust" else args[["engine"]]
if(!engine %in% c("rust","reference")) stop("--engine must be rust or reference")
if(engine=="rust") {
  script <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1])
  root <- dirname(dirname(normalizePath(script)))
  source(file.path(root,"R/limpa_rs.R"))
  if(!nzchar(Sys.getenv("LIMPA_RS_BIN"))) Sys.setenv(LIMPA_RS_BIN=file.path(root,"target/release/limpa-rs"))
  reference_cores <- if(is.null(args[["reference-cores"]])) 1L else as.integer(args[["reference-cores"]])
  dpc_quant_parallel <- function(y,protein_id,dpcfit,cores)
    dpc_quant_rust(y,protein_id,dpcfit,cores,reference.cores=reference_cores)
}
set.seed(if(is.null(args[["seed"]])) 141L else as.integer(args[["seed"]]))
dpcfit <- dpc(y)
t_quant <- system.time(yp <- dpc_quant_parallel(y, protein_id, dpcfit, cores))[["elapsed"]]
# The precursors dpc() fitted on (a random subset), with limpa's estimate of each one's
# complete-data mean intensity: the x-axis of limpa::plotDPC, used for the QC fit plot.
fwrite(data.table(precursor = names(dpcfit$mu), mu = unname(dpcfit$mu),
                  n_detected = unname(dpcfit$n.detected), nsamples = dpcfit$nsamples),
       file.path(outdir, "dpc_points.tsv"), sep = "\t")

write_matrix(yp$E, "protein_id", file.path(outdir, "protein_log2.tsv"))
write_matrix(yp$other$standard.error, "protein_id", file.path(outdir, "protein_se.tsv"))
ann <- data.table(protein_id = rownames(yp$E), as.data.table(yp$genes))
fwrite(ann, file.path(outdir, "protein_annotation.tsv"), sep = "\t", na = "NA")

# ---- sample table ---------------------------------------------------------------------------
targets <- data.frame(run = colnames(y$E), stringsAsFactors = FALSE)
if (!is.null(y$targets)) targets <- cbind(targets, y$targets[targets$run, , drop = FALSE])
if (!is.null(args[["samples"]])) {
  extra <- as.data.frame(fread(args[["samples"]], sep = "\t", check.names = FALSE))
  if (!"run" %in% names(extra)) stop("--samples TSV must have a 'run' column matching R.FileName")
  unknown <- setdiff(extra$run, targets$run)
  absent <- setdiff(targets$run, extra$run)
  if (length(unknown)) stop("--samples has runs not in the report: ", paste(unknown, collapse = ", "))
  if (length(absent)) stop("--samples is missing runs from the report: ", paste(absent, collapse = ", "))
  extra <- extra[, c("run", setdiff(names(extra), names(targets))), drop = FALSE]
  targets <- merge(targets, extra, by = "run", sort = FALSE)
  targets <- targets[match(colnames(y$E), targets$run), , drop = FALSE]
}
fwrite(targets, file.path(outdir, "samples.tsv"), sep = "\t", na = "NA")

# ---- optional differential expression --------------------------------------------------------
de_tables <- character(0)
if (!is.null(args[["formula"]])) {
  design <- model.matrix(as.formula(args[["formula"]]), data = targets)
  colnames(design) <- make.names(colnames(design))
  resid_df <- nrow(design) - qr(design)$rank
  if (resid_df < 1) {
    stop("Design '", args[["formula"]], "' leaves ", resid_df, " residual degrees of freedom for ",
         nrow(design), " runs: there are no replicates to estimate variance from. ",
         "Differential expression needs >= 2 runs per group (>= 3 recommended).")
  }
  fit <- dpcDE(yp, design, plot = FALSE)
  if (!is.null(args[["contrasts"]])) {
    cons <- trimws(strsplit(args[["contrasts"]], ",")[[1]])
    cm <- makeContrasts(contrasts = cons, levels = design)
    fit <- contrasts.fit(fit, cm)
    coefs <- colnames(cm)
  } else {
    coefs <- setdiff(colnames(design), "X.Intercept.")
  }
  fit <- eBayes(fit)
  for (cf in coefs) {
    tt <- topTable(fit, coef = cf, number = Inf, sort.by = "P")
    tt <- data.table(protein_id = rownames(tt), tt)
    fname <- paste0("de_", gsub("[^A-Za-z0-9._-]+", "_", cf), ".tsv")
    fwrite(tt, file.path(outdir, fname), sep = "\t", na = "NA")
    de_tables <- c(de_tables, setNames(fname, cf))
  }
}

# ---- optional comparator pipelines -----------------------------------------------------------
compare <- if (is.null(args[["compare"]])) character(0) else trimws(strsplit(args[["compare"]], ",")[[1]])
known <- c("maxlfq_proda")
if (length(setdiff(compare, known))) {
  stop("Unknown --compare method(s): ", paste(setdiff(compare, known), collapse = ", "),
       ". Known: ", paste(known, collapse = ", "))
}
cmp_tables <- character(0)

write_de <- function(tt, method, cf) {
  fname <- paste0("cmp_", method, "__", gsub("[^A-Za-z0-9._-]+", "_", cf), ".tsv")
  fwrite(tt, file.path(outdir, fname), sep = "\t", na = "NA")
  cmp_tables[paste0(method, "::", cf)] <<- fname
}

if (length(compare)) {
  suppressPackageStartupMessages(library(iq))
  # MaxLFQ on the same precursor log2 matrix and the same protein grouping limpa used.
  # No iq median normalization: Spectronaut already normalized EG.TotalQuantity (Settings).
  obs <- which(!is.na(y$E), arr.ind = TRUE)
  long <- list(
    protein_list = as.character(y$genes[[protein_id]][obs[, 1]]),
    sample_list = colnames(y$E)[obs[, 2]],
    id = rownames(y$E)[obs[, 1]],
    quant = y$E[obs]
  )
  mq <- iq::fast_MaxLFQ(long)$estimate
  mq <- mq[, colnames(y$E), drop = FALSE]
  write_matrix(mq, "protein_id", file.path(outdir, "maxlfq_log2.tsv"))

  if (is.null(args[["formula"]])) {
    message("--compare without --formula: wrote maxlfq_log2.tsv only (no DE to compare)")
  } else {
    if ("maxlfq_proda" %in% compare) {
      suppressPackageStartupMessages(library(proDA))
      keep <- rowSums(!is.na(mq)) > 0
      pfit <- proDA::proDA(mq[keep, , drop = FALSE], design = design)
      for (cf in coefs) {
        td <- proDA::test_diff(pfit, contrast = cf, sort_by = "pval")
        tt <- data.table(protein_id = td$name, logFC = td$diff, AveExpr = td$avg_abundance,
                         t = td$t_statistic, P.Value = td$pval, adj.P.Val = td$adj_pval,
                         se = td$se, df = td$df, n_obs = td$n_obs)
        write_de(tt, "maxlfq_proda", cf)
      }
    }
  }
}

# ---- summary ---------------------------------------------------------------------------------
summary <- list(
  engine = engine,
  reference_initialization_cores = if(engine=="rust") reference_cores else cores,
  seed = if(is.null(args[["seed"]])) 141L else as.integer(args[["seed"]]),
  input = input$source,
  cores = cores,
  quant_seconds = t_quant,
  protein_id_column = protein_id,
  q_cutoff = q_cutoff,
  imputed_values_filtered = has_imputed_flag,
  n_runs = ncol(y$E),
  n_precursors = nrow(y$E),
  n_proteins = nrow(yp$E),
  precursor_missing_fraction = mean(is.na(y$E)),
  observed_precursors_per_run = colSums(!is.na(y$E)),
  dpc_beta0 = unname(dpcfit$dpc[1]),
  dpc_beta1 = unname(dpcfit$dpc[2]),
  formula = args[["formula"]],
  de_tables = if (length(de_tables)) de_tables else NULL,
  compare_tables = if (length(cmp_tables)) cmp_tables else NULL,
  iq_version = if (length(compare)) as.character(packageVersion("iq")) else NULL,
  proda_version = if ("maxlfq_proda" %in% compare) as.character(packageVersion("proDA")) else NULL,
  limpa_version = as.character(packageVersion("limpa")),
  limma_version = as.character(packageVersion("limma")),
  r_version = R.version.string
)
writeLines(to_json(summary), file.path(outdir, "summary.json"))
message("limpa: ", nrow(y$E), " precursors -> ", nrow(yp$E), " proteins across ", ncol(y$E), " runs")
