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

#' exons
#' return \code{exons} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{exons}
#' @keywords internal
#' @rdname exons
setGeneric("exons", function(obj) standardGeneric("exons"))

#' juncs
#' return \code{juncs} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{juncs}
#' @keywords internal
#' @rdname juncs
setGeneric("juncs", function(obj) standardGeneric("juncs"))

#' exonValues
#' return \code{exonValues} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{exonValues}
#' @keywords internal
#' @rdname exonValues
setGeneric("exonValues", function(obj) standardGeneric("exonValues"))

#' juncValues
#' return \code{juncValues} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{juncValues}
#' @keywords internal
#' @rdname juncValues
setGeneric("juncValues", function(obj) standardGeneric("juncValues"))



## ###################################################################
## ###################################################################
## setter methods

#' exons replace
#' replace \code{exons} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{exons}
#' @keywords internal
#' @rdname exons-replace
setGeneric("exons<-", function(obj, ..., value) standardGeneric("exons<-"))

#' juncs replace
#' replace \code{juncs} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{juncs}
#' @keywords internal
#' @rdname juncs-replace
setGeneric("juncs<-", function(obj, ..., value) standardGeneric("juncs<-"))

#' exonValues replace
#' replace \code{exonValues} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{exonValues}
#' @keywords internal
#' @rdname exonValues-replace
setGeneric("exonValues<-", function(obj, ..., value) standardGeneric("exonValues<-"))

#' juncValues replace
#' replace \code{juncValues} slot
#' @export
#' @docType methods
#' @param obj object with slot \code{juncValues}
#' @keywords internal
#' @rdname juncValues-replace
setGeneric("juncValues<-", function(obj, ..., value) standardGeneric("juncValues<-"))
