#!/usr/bin/env Rscript

# Local CRAN-style check helper.
# - Disables bashisms check on non-Linux environments.
# - Leaves test/example execution enabled by default.

devtools::document(".")

devtools::check(
  ".",
  cran = TRUE,
  env_vars = c(
    `_R_CHECK_BASHISMS_` = "FALSE",
    `_R_CHECK_FUTURE_FILE_TIMESTAMPS_` = "FALSE"
  )
)
