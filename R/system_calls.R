#' System calls and stderr stdout reporting
#'
#' Robustly pass commands to system() and capture stderr and stdout data
#'
#' @param command [Required] String that is passed directly to the base::system() command.
#'
#' 
#'
#' @return A list of length 3. list[[1]] is the exit code for the command. list[[2]] is the "stdout" text. list[[3]] is the "stderr" text.
#'
#'
#'
#' This function was mostly copied from:
#' [https://stackoverflow.com/a/17090601/2367748](https://stackoverflow.com/a/17090601/2367748)
#'
#' @examples
#' robust_list <- robust_system("date; cal")
#' 
#' @seealso [robust_system_check()]
#'
#' @export
robust_system <- function(command) {
    stdoutFile <- tempfile(pattern = "R_robust_system_stdout.", fileext = as.character(Sys.getpid()))
    stderrFile <- tempfile(pattern = "R_robust_system_stderr.", fileext = as.character(Sys.getpid()))
    
    cat("robust_system temp files (below) are deleted when command completes.", sep = "\n")
    cat(command, sep = "\n")
    cat(stdoutFile, sep = "\n")
    cat(stderrFile, sep = "\n")
    
    retval <- list()
    retval$command <- command
    retval$exit_status <- system(paste0(command, " 2> ", shQuote(stderrFile), " > ", shQuote(stdoutFile)))
    retval$stdout <- readLines(stdoutFile)
    retval$stderr <- readLines(stderrFile)


    unlink(c(stdoutFile, stderrFile))
    return(retval)
}






#' Print the data captured by robust_system() and stop if nonzero exit code
#'
#' @param robust_list [Required] List object that is returned by the [robust_system()] function.
#' @param name [Optional. Default: name of `robust_list` object] String that provides a name for the system call captured in the `robust_list`.
#' 
#'
#' @return This function prints the contents of the `robust_list` in a nice way, directly to stderr and stdout streams. It also checks the exit code and will `stop()` if not zero. 
#'
#'
#' @examples
#' robust_system_check(robust_list)
#' 
#' @seealso [robust_system()]
#'
#' @export
robust_system_check <- function(robust_list = NULL, name = deparse(substitute(robust_list))) {
    # Write info to stdout
    cat_stdout <- function(x) {
        if (is.null(x)) {
            cat("NULL", sep = "\n", append = TRUE, file = stdout())
        } else {
            cat(x, sep = "\n", append = TRUE, file = stdout())
        }
    }
    cat_stderr <- function(x) {
        if (is.null(x)) {
            cat("NULL", sep = "\n", append = TRUE, file = stderr())
        } else {
            cat(x, sep = "\n", append = TRUE, file = stderr())
        }
    }
    
    # Write info to stdout
    cat_stdout("#######################################################################")
    cat_stdout(name)
    cat_stdout("#######################################################################")
    cat_stdout("---------------------------- command ----------------------------------")
    cat_stdout(robust_list$command)
    cat_stdout("---------------------------- exit_status ------------------------------")
    cat_stdout(robust_list$exit_status)
    cat_stdout("---------------------------- stdout -----------------------------------")
    cat_stdout(robust_list$stdout)


    # Write info to stderr
    cat_stderr("#######################################################################")
    cat_stderr(name)
    cat_stderr("#######################################################################")
    cat_stderr("---------------------------- command ----------------------------------")
    cat_stderr(robust_list$command)
    cat_stderr("---------------------------- exit_status ------------------------------")
    cat_stderr(robust_list$exit_status)
    cat_stderr("---------------------------- stderr -----------------------------------")
    cat_stderr(robust_list$stderr)




    if (robust_list$exit_status != "0") {
        stop(paste0("ERROR: ", name, " section failed with non-zero exit code"))
    }
}


