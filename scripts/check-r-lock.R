# Fail before numerical tests if the actual R environment drifts from renv.lock.
args <- commandArgs(trailingOnly = TRUE)
lock_path <- if (length(args)) args[[1L]] else "renv.lock"
lock <- renv::lockfile_read(lock_path)
errors <- character()
if (as.character(getRversion()) != lock$R$Version)
  errors <- c(errors, sprintf("R version: expected %s, installed %s", lock$R$Version, getRversion()))
for (name in names(lock$Packages)) {
  expected <- lock$Packages[[name]]
  actual <- tryCatch(utils::packageDescription(name), error = function(e) NULL)
  if (is.null(actual) || !is.list(actual) || !identical(actual$Version, expected$Version)) {
    errors <- c(errors, sprintf("%s: installed version does not match %s", name, expected$Version))
    next
  }
  # Repository-built packages do not always preserve RemoteSha. Compare it when
  # present; git/GitHub installs must retain it. This is not an archive hash check.
  if (!is.null(expected$RemoteSha)) {
    if (!is.null(actual$RemoteSha) && !identical(actual$RemoteSha, expected$RemoteSha))
      errors <- c(errors, paste(name, "installed source SHA differs from lock"))
    if (expected$Source %in% c("Git", "GitHub") && is.null(actual$RemoteSha))
      errors <- c(errors, paste(name, "installed source SHA is missing"))
  }
}
if (length(errors)) stop(paste(errors, collapse = "\n"), call. = FALSE)
cat("Verified", length(lock$Packages), "installed R package versions against renv.lock\n")
