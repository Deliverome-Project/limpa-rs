# SPDX-License-Identifier: GPL-3.0-or-later
if (getRversion() < "4.6" || getRversion() >= "4.7")
  stop("limpa-rs requires the assessed R 4.6.x environment")
source(file.path(Sys.getenv("RENV_PROJECT"), "renv/activate.R"))
