args <- commandArgs(trailingOnly = TRUE)

min_cov <- Sys.getenv("VNORM_COVERAGE_MIN", unset = "0")
if (length(args) > 0) {
  min_cov <- args[[1]]
}
min_cov <- suppressWarnings(as.numeric(min_cov))
if (!is.finite(min_cov)) {
  stop("Coverage threshold must be numeric.", call. = FALSE)
}

if (!requireNamespace("covr", quietly = TRUE)) {
  stop("Package 'covr' is required. Install it with install.packages('covr').", call. = FALSE)
}

install_path <- tempfile(pattern = "vnorm-covr-lib-")

cov <- covr::package_coverage(type = "tests", install_path = install_path)
pct <- covr::percent_coverage(cov)
xml_out <- Sys.getenv("VNORM_COVERAGE_XML", unset = "")

cat(sprintf("Total coverage: %.2f%%\n", pct))
cov_df <- as.data.frame(cov)
if (nrow(cov_df) > 0 && all(c("filename", "value") %in% names(cov_df))) {
  by_file <- aggregate(value ~ filename, data = cov_df, FUN = function(x) {
    mean(x > 0) * 100
  })
  by_file <- by_file[order(by_file$filename), , drop = FALSE]
  names(by_file) <- c("file", "line_coverage_pct")
  cat("\nApproximate line coverage by file (%):\n")
  print(by_file, row.names = FALSE)
}

if (pct < min_cov) {
  stop(
    sprintf(
      "Coverage %.2f%% is below required threshold %.2f%%.",
      pct, min_cov
    ),
    call. = FALSE
  )
}

if (nzchar(xml_out)) {
  tryCatch(
    {
      covr::to_cobertura(cov, filename = xml_out)
    },
    error = function(e) {
      xml <- covr::to_cobertura(cov)
      writeLines(as.character(xml), con = xml_out)
    }
  )
  cat(sprintf("\nWrote Cobertura XML: %s\n", xml_out))
}

invisible(cov)
