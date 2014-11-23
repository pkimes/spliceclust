#' read in data for single chromosome
#'
#' Function to read in data for single chromosome output from
#' \code{diffsplice_bam} function provided by Yin Hu.
#'
#' @param file character string to data directory
#' @param n integer value specifying number of samples
#' 
#' @return
#' a \code{data.frame} object containing splicing information in
#' \code{6 + n} columns. The columns are:
#' \code{gIdx}, \code{gStart}, \code{gStop}, \code{kind}, \code{start},
#' \code{stop}, \code{s1}, \code{s2}, ...
#' 
#' @export
#' @author Patrick Kimes

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

