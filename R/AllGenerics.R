## ###################################################################
## ###################################################################
## constructor methods

#' @export
#' @rdname concomp-class
setGeneric("concomp",
           valueClass="concomp",
           function(obj, ...) standardGeneric("concomp"))



## ###################################################################
## ###################################################################
## plot methods

#' @export
#' @rdname splicegrahm
#' @docType methods
setGeneric("splicegrahm",
           function(obj, ...)  standardGeneric("splicegrahm"))

#' @export
#' @rdname splicegrahm2
#' @docType methods
setGeneric("splicegrahm2",
           function(obj1, obj2, ...)  standardGeneric("splicegrahm2"))

#' @export
#' @docType methods
#' @rdname splicegralp
setGeneric("splicegralp",
           function(obj, ...) standardGeneric("splicegralp"))

#' @export
#' @docType methods
#' @rdname splicepca
setGeneric("splicepca",
           function(obj, ...) standardGeneric("splicepca"))

#' @export
#' @docType methods
#' @rdname splicepcp
setGeneric("splicepcp",
           function(obj, ...) standardGeneric("splicepcp"))



## ###################################################################
## ###################################################################
## getter methods

#' @export
#' @rdname slots-generic
setGeneric("eRanges", function(obj) standardGeneric("eRanges"))

#' @export
#' @rdname slots-generic
setGeneric("jRanges", function(obj) standardGeneric("jRanges"))

#' @export
#' @rdname slots-generic
setGeneric("eCoverage", function(obj) standardGeneric("eCoverage"))

#' @export
#' @rdname slots-generic
setGeneric("jCoverage", function(obj) standardGeneric("jCoverage"))



## ###################################################################
## ###################################################################
## setter methods

#' @export
#' @rdname slots-generic
setGeneric("eRanges<-", function(obj, ..., value) standardGeneric("eRanges<-"))

#' @export
#' @rdname slots-generic
setGeneric("jRanges<-", function(obj, ..., value) standardGeneric("jRanges<-"))

#' @export
#' @rdname slots-generic
setGeneric("eCoverage<-", function(obj, ..., value) standardGeneric("eCoverage<-"))

#' @export
#' @rdname slots-generic
setGeneric("jCoverage<-", function(obj, ..., value) standardGeneric("jCoverage<-"))
