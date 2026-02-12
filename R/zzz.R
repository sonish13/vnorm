.onLoad <- function(libname, pkgname) {
  options(cmdstanr_no_ver_check = TRUE)
}

.onAttach <- function(libname, pkgname) {
  suppressPackageStartupMessages({

    if (requireNamespace("instantiate", quietly = TRUE)) {
      require("instantiate", character.only = TRUE)
    }

    if (requireNamespace("mpoly", quietly = TRUE)) {
      require("mpoly", character.only = TRUE)
    }

    if (requireNamespace("cmdstanr", quietly = TRUE)) {
      ok <- TRUE
      ok <- ok && !identical(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL), NULL)
      ok <- ok && !identical(cmdstanr::cmdstan_version(), "none")

      if (isTRUE(ok)) {
        require("cmdstanr", character.only = TRUE)
      }
    }
  })
}
