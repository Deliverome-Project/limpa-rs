# SPDX-License-Identifier: GPL-3.0-or-later
reference_numeric_error <- function(a, b) {
  if (is.null(a) && is.null(b)) return(0)
  if (!is.numeric(a) || !is.numeric(b) || !identical(dim(a), dim(b)) ||
      !identical(dimnames(a), dimnames(b)) || !identical(names(a), names(b)) || length(a) != length(b) ||
      !identical(is.na(a), is.na(b)) || !identical(is.nan(a), is.nan(b)) ||
      !identical(is.infinite(a), is.infinite(b))) return(Inf)
  inf <- is.infinite(a)
  if (!identical(a[inf], b[inf])) return(Inf)
  finite <- is.finite(a)
  if (!any(finite)) return(0)
  max(abs(a[finite] - b[finite]))
}

compare_reference_outputs <- function(reference, candidate) {
  if (!identical(names(reference), names(candidate))) stop("Output fields differ")
  e <- reference_numeric_error(reference$E, candidate$E)
  se <- reference_numeric_error(reference$se, candidate$se)
  metadata <- identical(reference$metadata, candidate$metadata)
  # Fixed acceptance limits; changing these is a separate scientific decision.
  data.frame(max_log2_error=e, max_se_error=se, metadata_equal=metadata,
             pass=is.finite(e) && is.finite(se) && e <= 1e-3 && se <= 1e-3 && metadata)
}

validate_reference_snapshot <- function(x) {
  if (!identical(x$schema, 1L) || !is.list(x$provenance) ||
      !all(c("limpa", "R", "limma", "engine") %in% names(x$provenance)) ||
      !is.list(x$cases) || !length(x$cases) || anyDuplicated(names(x$cases)) ||
      !identical(names(x$cases), names(x$outputs))) stop("Invalid reference snapshot")
  for (out in x$outputs) {
    if (!identical(names(out), c("E","se","metadata")) ||
        !is.matrix(out$E) || !is.numeric(out$E) || !length(out$E) ||
        !is.list(out$metadata) || !identical(names(out$metadata),
          c("genes","targets","n.observations","dpc","prior.mean","prior.sd","prior.logFC")))
      stop("Invalid reference output schema")
  }
  invisible(TRUE)
}

compare_reference_snapshots <- function(reference, candidate) {
  validate_reference_snapshot(reference); validate_reference_snapshot(candidate)
  if (!identical(reference$cases, candidate$cases)) stop("Fixture inputs differ")
  do.call(rbind, lapply(names(reference$cases), function(id) {
    cbind(case=id, reference_limpa=reference$provenance$limpa,
      candidate_limpa=candidate$provenance$limpa,
      compare_reference_outputs(reference$outputs[[id]], candidate$outputs[[id]]))
  }))
}
