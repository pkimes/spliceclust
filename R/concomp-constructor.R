.concomp.GRangesList <- function(obj) {  
    if (!setequal(names(obj), c("e", "j")))
        stop("GRangesList must only contain two objects: 'e' and 'j'")

    gr_e <- sort(obj[["e"]])
    gr_j <- sort(obj[["j"]])

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
        eRanges = gr_e,
        jRanges = gr_j,
        eCoverage = as.matrix(vals_e),
        jCoverage = as.matrix(vals_j))
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
    gr_e <- sort(gr[obj$kind == "e"])
    gr_j <- sort(gr[obj$kind == "j"])

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
        eRanges = gr_e,
        jRanges = gr_j,
        eCoverage = vals_e,
        jCoverage = vals_j)
}


#' @rdname concomp-class
setMethod("concomp",
          signature(obj = "GRangesList"),
          .concomp.GRangesList)


#' @rdname concomp-class
setMethod("concomp",
          signature(obj = "data.frame"),
          .concomp.data.frame)
