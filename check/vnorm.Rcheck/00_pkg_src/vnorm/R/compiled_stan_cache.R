.vnorm_state <- new.env(parent = emptyenv())

get_compiled_stan_info <- function() {
  # Read cache metadata (name/path) for user-compiled Stan models.
  info <- .vnorm_state$compiled_stan_info
  if (is.null(info)) {
    return(data.frame(name = character(), path = character(), stringsAsFactors = FALSE))
  }
  info
}

set_compiled_stan_info <- function(info) {
  # Replace the full cache table in package-private state.
  .vnorm_state$compiled_stan_info <- info
  invisible(info)
}

add_compiled_stan_info <- function(name, path) {
  # Append a model entry and keep unique rows.
  info <- get_compiled_stan_info()
  new_row <- data.frame(name = name, path = path, stringsAsFactors = FALSE)
  info <- unique(rbind(info, new_row))
  set_compiled_stan_info(info)
}

clear_compiled_stan_info <- function() {
  # Drop cached model metadata.
  .vnorm_state$compiled_stan_info <- NULL
  invisible(NULL)
}
