#' read in single chromosome data
#'
#' Function to read in data for single chromosome output from
#' \code{diffsplice_bam} function from Yin Hu.
#' Last updated on 09/08/2014.
#'
#' @param file character string to data directory
#' @param n integer value specifying number of samples
#' 
#' @return a \code{data.frame} object containing splicing information for data
#' 
#' @export
#' 
#' @author Patrick

readchr <- function(file, n) {
    ##for faster loading, specify column types
    col_format <- c("character", "numeric", "numeric",
                    "character", "numeric", "numeric",
                    rep("numeric", n))
    ##for consistency, provide column names
    col_names <- c("gIdx", "gStart", "gStop", "kind", "start", "stop",
                   paste0("s", 1:n))
    
    x <- read.table(file, colClasses=col_format)
    names(x) <- col_names
    
    x
}

