% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{write_session_info}
\alias{write_session_info}
\title{Write out current R session information}
\usage{
write_session_info(out_file_prefix)
}
\arguments{
\item{out_file_prefix}{File path, including a filename prefix if desired.}
}
\value{
None
}
\description{
The function writes our R session info, with a time and date stamped into the filename. The function will overwrite previous files with updated timestamp.
}
\examples{
\dontrun{
out_file_prefix <- "code_out/dir/050_"
write_session_info(out_file_prefix)
}

}
