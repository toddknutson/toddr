#' Write out current R session information
#'
#' The function writes our R session info, with a time and date stamped into the filename. The function will overwrite previous files with updated timestamp.
#'
#' @param out_file_prefix File path, including a filename prefix if desired.
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' out_file_prefix <- "code_out/dir/050_"
#' write_session_info(out_file_prefix)
#' }
#'
#' @export
write_session_info <- function(out_file_prefix) {
	out_dir <- dirname(as.character(out_file_prefix))

	# If previous session info is present, delete.
	if (length(grep("session_info", list.files(out_dir))) > 0) { 
	invisible(file.remove(grep("session_info", list.files(out_dir, full.names = TRUE), value = TRUE)))
	}
	# Write out session info
	writeLines(capture.output(devtools::session_info()), paste0(out_file_prefix, "session_info_", gsub(" |:", "_", Sys.time()), ".txt"))
	print(Sys.time())
}






#' Release memory
#'
#' Run after deleting large objects to initiate garbage collection multiple times.
#' [Reference](https://stackoverflow.com/a/1467334/2367748)

#' @param n Integer, Number of times to run garbage collection. Default is 10.
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' clean_mem()
#' }
#'
#' @export
clean_mem <- function(n = 10) {
    for (i in 1:n) gc()
}



