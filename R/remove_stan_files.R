#' Remove files created by [compile_stan_code()]
#'
#' Deletes all the variables and Stan files created by [compile_stan_code()].
#'
#' @param path Path where [compile_stan_code()] was used and files need deletion.
#' Defaults to the current working directory.
#' @usage remove_stan_files()
#'
#'
#' @export
remove_stan_files <- function(path = getwd()) {
  #check to see if the files to be deleted are present
  if (exists("compiled_stan_info", envir = .GlobalEnv)) {
    stan_file_names <- compiled_stan_info$path
    exe_file_names <- sub("\\.stan$", "", stan_file_names)
    #remove the files
    file.remove(stan_file_names[file.exists(stan_file_names)])
    file.remove(exe_file_names[file.exists(exe_file_names)])

    remove(compiled_stan_info, envir = .GlobalEnv)
  }
  message("Files deleted")
}
