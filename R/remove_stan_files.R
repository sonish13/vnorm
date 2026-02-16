#' Remove files created by [compile_stan_code()]
#'
#' Deletes all the variables and Stan files created by [compile_stan_code()].
#'
#' @param path Path where [compile_stan_code()] was used and files need deletion.
#' Defaults to the current working directory.
#'
#'
#' @export
remove_stan_files <- function(path = getwd()) {
  compiled_stan_info <- get_compiled_stan_info()

  if (nrow(compiled_stan_info) > 0) {
    stan_file_names <- compiled_stan_info$path
    exe_file_names <- sub("\\.stan$", "", stan_file_names)

    file.remove(stan_file_names[file.exists(stan_file_names)])
    file.remove(exe_file_names[file.exists(exe_file_names)])
    clear_compiled_stan_info()
  }
  message("Files deleted")
}
