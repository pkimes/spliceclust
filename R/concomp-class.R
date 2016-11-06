#' Connected Component
#'
#' S4 class encapsulating exon and junction information for single
#' connected component
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
#' @param obj a \code{GRangesList} with two \code{GRanges} named "e" and "j"
#'        corresponding to exon and junction ranges and coverage for a
#'        collection of samples to be used to construct a \code{concomp} object.
#'        Alternative, can be a \code{data.frame} with gene model columns:
#'        'chr', 'seqlengths', 'gStart', 'gStop', and 'kind' (e.g. 'e' or 'j'), 
#'        as well as per-sample coverage columns: 's1', 's2', etc.
#' @param ... other parameters to be passed to constructor
#' 
#' @param object a \code{concomp} object for show method
#'
#' @param obj a \code{concomp} object for getters and setters
#' @param value new value for \code{concomp} slot
#' 
#' @details
#' The columns of \code{eCoverage} and \code{jCoverage} correspond
#' to individual samples or experiments, and should be of equal length.
#'
#' @name concomp-class
#' @import GenomicRanges IRanges
#' @exportClass concomp
#' @author Patrick Kimes
NULL

