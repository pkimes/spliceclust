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
           function(obj, sort_sep = FALSE, sort_idx = 1,
                    log_base = 10, log_shift = 1, bin = TRUE,
                    genomic = TRUE, ex_use = 2/3, flip_neg = TRUE, 
                    j_incl = FALSE, highlight = NULL,
                    use_blk = FALSE, eps = 1e4, txlist = NULL,
                    txdb = NULL, orgdb = NULL, title="", ...)
               standardGeneric("splicegrahm"))

#' @export
#' @rdname splicegrahm2
#' @docType methods
setGeneric("splicegrahm2",
           function(obj1, obj2, sort_sep = FALSE,
                    sort_idx1 = 1, sort_idx2 = 1,
                    log_base = 10, log_shift = 1, bin = TRUE,
                    genomic = TRUE, ex_use = 2/3, flip_neg = TRUE, 
                    j_incl = FALSE, use_blk = FALSE, eps = 1e4,
                    txlist = NULL, txdb = NULL, orgdb = NULL, title="",
                    mirror = TRUE, same_scale = TRUE, ...)
               standardGeneric("splicegrahm2"))

#' @export
#' @docType methods
#' @rdname splicegralp
setGeneric("splicegralp",
           function(obj, e_loads, j_loads = NULL, load_lims = NULL, 
                    genomic = TRUE, ex_use = 2/3,
                    flip_neg = TRUE, use_blk = FALSE,
                    txlist = NULL, txdb = NULL,
                    orgdb = NULL,...)
               standardGeneric("splicegralp"))

#' @export
#' @docType methods
#' @rdname splicepca
setGeneric("splicepca",
           function(obj, npc = 3, pc_sep = TRUE, ej_w = c(1, 1),
                    log_base = 10, log_shift = 1,
                    genomic = TRUE, ex_use = 2/3,
                    flip_neg = TRUE, use_blk = FALSE,
                    txlist = NULL, txdb = NULL,
                    orgdb = NULL, scores = FALSE,
                    plot = TRUE, ...)
               standardGeneric("splicepca"))

#' @export
#' @docType methods
#' @rdname splicepcp
setGeneric("splicepcp",
           function(obj, log_base = 10, log_shift = 1, genomic = TRUE,
                    ex_use = 2/3, flip_neg = TRUE, imodel = TRUE, 
                    highlight = NULL, eps = 1e4,
                    txlist = NULL, txdb = NULL, orgdb = NULL, ...)
               standardGeneric("splicepcp"))



## ###################################################################
## ###################################################################
## getter methods

#' @export
#' @rdname slots-generic
setGeneric("eRanges", function(object) standardGeneric("eRanges"))

#' @export
#' @rdname slots-generic
setGeneric("jRanges", function(object) standardGeneric("jRanges"))

#' @export
#' @rdname slots-generic
setGeneric("eCoverage", function(object) standardGeneric("eCoverage"))

#' @export
#' @rdname slots-generic
setGeneric("jCoverage", function(object) standardGeneric("jCoverage"))



## ###################################################################
## ###################################################################
## setter methods

#' @export
#' @rdname slots-generic
setGeneric("eRanges<-", function(object, ..., value) standardGeneric("eRanges<-"))

#' @export
#' @rdname slots-generic
setGeneric("jRanges<-", function(object, ..., value) standardGeneric("jRanges<-"))

#' @export
#' @rdname slots-generic
setGeneric("eCoverage<-", function(object, ..., value) standardGeneric("eCoverage<-"))

#' @export
#' @rdname slots-generic
setGeneric("jCoverage<-", function(object, ..., value) standardGeneric("jCoverage<-"))
