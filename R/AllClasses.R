
## construct union class to allow for empty eCoverage jCoverage slots
setOldClass("matrix")
setClassUnion("MatrixOrNULL", c("matrix", "NULL"))


#' connected component object
#'
#' S4 class encapsulating exon and junction
#' information for single connected component
#'
#' @slot eRanges a \code{GRanges} object containing the exonic information
#' @slot jRanges a \code{GRanges} object containing the junction information
#' @slot eCoverage a matrix of exon coverage values, e.g. RPKM, with
#'       number of rows equal to the length of \code{exons}, or \code{NULL}
#'       if not available
#' @slot jCoverage a matrix of junction coverage values, e.g. number of
#'       spanning reads, with number of rows equal to the length of
#'       \code{juncs}, or \code{NULL} if not available
#'
#' @param obj a \code{concomp} object
#' @param value new value for \code{concomp} slot
#' 
#' @details
#' The columns of \code{eCoverage} and \code{jCoverage} correspond
#' to individual samples or experiments, and should be of equal length.
#'
#' @exportClass concomp
#' @name concomp-class
#' @aliases concomp-class
#' @rdname concomp-class
#' @author Patrick Kimes

setClass("concomp",
         slots=list(
             eRanges = "GRanges",
             jRanges = "GRanges",
             eCoverage = "MatrixOrNULL",
             jCoverage = "MatrixOrNULL"))

