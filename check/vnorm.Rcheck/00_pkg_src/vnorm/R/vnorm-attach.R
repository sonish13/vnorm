.onLoad <- function(libname, pkgname) {
  # Silence cmdstanr version-check startup chatter for package users.
  options(cmdstanr_no_ver_check = TRUE)
}
