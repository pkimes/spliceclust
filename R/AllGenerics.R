## ###################################################################
## ###################################################################
## constructor methods

#' @name concomp
#' @export
#' @docType methods
#' @rdname concomp-constructor
setGeneric("concomp",
           valueClass="concomp",
           function(obj, ...) standardGeneric("concomp"))



## ###################################################################
## ###################################################################
## plot methods

#' Splice Graph HeatMap (SpliceGraHM)
#'
#' Method for plotting both the exon and junction values of a single
#' connected component
#'
#' @export
#' @docType methods
#' @name splicegrahm-generic
#' @keywords internal
setGeneric("splicegrahm",
           function(obj, ...)  standardGeneric("splicegrahm"))


#' Two-class Splice Graph HeatMap (SpliceGraHM2)
#'
#' Method for plotting both the exon and junction values of a single
#' connected component for two separate clusters or populations
#'
#' @export
#' @docType methods
#' @name splicegrahm2-generic
#' @keywords internal
setGeneric("splicegrahm2",
           function(obj1, obj2, ...)  standardGeneric("splicegrahm2"))


#' Splice Graph Loadings Plot (SpliceGraLP)
#'
#' Method for plotting different loading vectors for a \code{concomp}
#' object, e.g. principal component directions
#'
#' @export
#' @docType methods
#' @name splicegralp-generic
#' @keywords internal
setGeneric("splicegralp",
           function(obj, ...) standardGeneric("splicegralp"))


#' Splice Graph PCA Plot
#'
#' Method for plotting PCA loading vectors for a \code{concomp}
#' object (wrapper to splicegralp)
#'
#' @export
#' @docType methods
#' @name splicepca-generic
#' @keywords internal
setGeneric("splicepca",
           function(obj, ...) standardGeneric("splicepca"))


#' Splice Graph Parallel Coordinates Plot
#'
#' Method for plotting exon expression for a \code{concomp}
#' object with each exon given equal width with heights corresponding
#' to expression level
#'
#' @export
#' @docType methods
#' @name splicepcp-generic
#' @keywords internal
setGeneric("splicepcp",
           function(obj, ...) standardGeneric("splicepcp"))



## ###################################################################
## ###################################################################
## getter methods

#' eRanges
#' return \code{eRanges} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{eRanges}
#' @keywords internal
#' @rdname eRanges
setGeneric("eRanges", function(obj) standardGeneric("eRanges"))

#' jRanges
#' return \code{jRanges} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{jRanges}
#' @keywords internal
#' @rdname jRanges
setGeneric("jRanges", function(obj) standardGeneric("jRanges"))

#' eCoverage
#' return \code{eCoverage} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{eCoverage}
#' @keywords internal
#' @rdname eCoverage
setGeneric("eCoverage", function(obj) standardGeneric("eCoverage"))

#' jCoverage
#' return \code{jCoverage} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{jCoverage}
#' @keywords internal
#' @rdname jCoverage
setGeneric("jCoverage", function(obj) standardGeneric("jCoverage"))



## ###################################################################
## ###################################################################
## setter methods

#' eRanges replace
#' replace \code{eRanges} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{eRanges}
#' @keywords internal
#' @rdname eRanges-replace
setGeneric("eRanges<-", function(obj, ..., value) standardGeneric("eRanges<-"))

#' jRanges replace
#' replace \code{jRanges} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{jRanges}
#' @keywords internal
#' @rdname jRanges-replace
setGeneric("jRanges<-", function(obj, ..., value) standardGeneric("jRanges<-"))

#' eCoverage replace
#' replace \code{eCoverage} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{eCoverage}
#' @keywords internal
#' @rdname eCoverage-replace
setGeneric("eCoverage<-", function(obj, ..., value) standardGeneric("eCoverage<-"))

#' jCoverage replace
#' replace \code{jCoverage} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{jCoverage}
#' @keywords internal
#' @rdname jCoverage-replace
setGeneric("jCoverage<-", function(obj, ..., value) standardGeneric("jCoverage<-"))
