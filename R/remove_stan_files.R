#' Remove files created by[compile_stan_code()]
#'
#' Deletes all the variables and stan files created by [compile_stan_code()]
#'
#'
#' @param path Path where [compile_stan_code()] was used and files need deletion.
#' Defaults to curerent working directory.
#' @usage remove_stan_file()
#'
#'
#' @export
remove_stan_files <- function(path = getwd()) {
  if(exists("compiled_stan_info", envir = .GlobalEnv)){
  stan_file_names <-  compiled_stan_info$path
  exe_file_names <- sub("\\.stan$", "", stan_file_names)
  file.remove(stan_file_names[file.exists(stan_file_names)])
  file.remove(exe_file_names[file.exists(exe_file_names)])
  remove(compiled_stan_info, envir = .GlobalEnv)
  }
  return("Files deleted")
}
