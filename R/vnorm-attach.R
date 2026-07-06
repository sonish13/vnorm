.onLoad <- function(libname, pkgname) {
  # silence cmdStanr version-check startup chatter for package users
  if (is.null(getOption("cmdstanr_no_ver_check"))) {
    options(cmdstanr_no_ver_check = TRUE)
  }
}
