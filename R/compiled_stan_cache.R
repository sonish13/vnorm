.vnorm_state <- new.env(parent = emptyenv())

get_compiled_stan_info <- function() {
  info <- .vnorm_state$compiled_stan_info
  if (is.null(info)) {
    return(data.frame(name = character(), path = character(), stringsAsFactors = FALSE))
  }
  info
}

set_compiled_stan_info <- function(info) {
  .vnorm_state$compiled_stan_info <- info
  invisible(info)
}

add_compiled_stan_info <- function(name, path) {
  info <- get_compiled_stan_info()
  new_row <- data.frame(name = name, path = path, stringsAsFactors = FALSE)
  info <- unique(rbind(info, new_row))
  set_compiled_stan_info(info)
}

clear_compiled_stan_info <- function() {
  .vnorm_state$compiled_stan_info <- NULL
  invisible(NULL)
}
