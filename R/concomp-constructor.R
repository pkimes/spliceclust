#' Connected Component (concomp) object
#'
#' construct a connected component object from various input objects
#' 
#' @return \code{concomp} object
#'
#' @import GenomicRanges IRanges
#' @name concomp-constructor
#' @author Patrick Kimes
NULL


.concomp.GRangesList <- function(obj) {  
    if (!setequal(names(obj), c("e", "j")))
        stop("GRangesList must only contain two objects: 'e' and 'j'")

    gr_e <- obj[["e"]]
    gr_j <- obj[["j"]]

    vals_e <- mcols(gr_e)
    vals_j <- mcols(gr_j)

    vals_e <- vals_e[grep("^s[0-9]", names(vals_e), value=TRUE)]
    vals_j <- vals_j[grep("^s[0-9]", names(vals_j), value=TRUE)]

    if (any(names(vals_e) != names(vals_j)))
        stop("'e' and 'j' GRanges must have same samples in same order")

    vals_e <- as.data.frame(vals_e)
    vals_j <- as.data.frame(vals_j)

    mcols(gr_e) <- NULL
    mcols(gr_j) <- NULL
    
    new("concomp",
        exons = gr_e,
        juncs = gr_j,
        exonValues = as.matrix(vals_e),
        juncValues = as.matrix(vals_j))
}


.concomp.data.frame <- function(obj) {
    if (!("chr" %in% names(obj)))
        stop("'chr' column must be specified with e.g. 'chr9'.")

    if (!("kind" %in% names(obj)))
        stop(paste0("'kind' column must be specified with 'e' and 'j'",
                    " to differentiate exons and juncs."))

    if (!("seqlengths" %in% names(obj)))
        stop(paste0("'seqlengths' column must be specified, see manual for",
                    " directions on how to obtain appropriate values,",
                    " note that only the first value is used."))
    
    gr <- makeGRangesFromDataFrame(obj, seqnames.field="chr")
    seqlengths(gr) <- obj$seqlengths[1]
    gr_e <- gr[obj$kind == "e"]
    gr_j <- gr[obj$kind == "j"]

    val_cols <- grep("^s[0-9]", names(obj), value=TRUE)
    if (length(val_cols) > 0) {
        vals <- as.matrix(obj[val_cols])
        vals_e <- vals[obj$kind == "e", , drop=FALSE]
        vals_j <- vals[obj$kind == "j", , drop=FALSE]
    } else {
        vals_e <- NULL
        vals_j <- NULL
    }
        
    new("concomp",
        exons = gr_e,
        juncs = gr_j,
        exonValues = vals_e,
        juncValues = vals_j)
}


#' @param obj a \code{GRangesList} with two elements named "e" and "j"
#' @rdname concomp-constructor
setMethod("concomp",
          signature(obj = "GRangesList"),
          .concomp.GRangesList)


#' @param obj a \code{data.frame} containing all of the information
#' @rdname concomp-constructor
setMethod("concomp",
          signature(obj = "data.frame"),
          .concomp.data.frame)
