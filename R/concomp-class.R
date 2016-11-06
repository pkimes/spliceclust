#' Connected Component
#'
#' S4 class encapsulating exon and junction information for a single
#' "connected component" object. 
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
#' @param obj primary input to \code{concomp} constructor (see Details for
#'        more information on supported signatures)
#' @param ... other parameters to be passed to constructor
#' 
#' @param object a \code{concomp} object for slot setters and getters, and show method
#' @param value new value for \code{concomp} slot
#' 
#' @details
#' The constructor method \code{obj} parameter can be either:
#' \itemize{
#' \item a \code{GRangesList} with two \code{GRanges} named "e" and "j"
#'       corresponding to exon and junction ranges and coverage for a
#'       collection of samples to be used to construct a \code{concomp} object.
#'       Both \code{GRanges} must have one metadata for each sample
#'       (named sequentially: 's1', 's2', ...), giving the corresponding exon or
#'       splice junction coverages. 
#' \item a \code{data.frame} with gene model columns: 'chr', 'seqlengths',
#'       'gStart', 'gStop', and 'kind' (e.g. 'e' or 'j'), 
#'       as well as per-sample coverage columns: 's1', 's2', ...
#' }
#'
#' @name concomp-class
#' @import GenomicRanges IRanges
#' @exportClass concomp
#' @author Patrick Kimes
NULL

