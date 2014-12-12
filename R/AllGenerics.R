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
#' @name splicePCA-generic
#' @keywords internal
setGeneric("splicePCA",
           function(obj, ...) standardGeneric("splicePCA"))



## ###################################################################
## ###################################################################
## analysis methods

#' cluster
#' cluster object and return class labels
#' @export
#' @docType methods
#' @param obj object which can be clustered
#' @keywords internal
setGeneric("cluster",
           function(obj, ...) standardGeneric("cluster"))



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
