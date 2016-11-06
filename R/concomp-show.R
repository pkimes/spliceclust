.show.concomp <- function(object) {
    cat("\n concomp object summary: \n")
    cat(paste0("    n sampl: ", ncol(eCoverage(object)), "\n"))
    cat(paste0("    p exons: ", length(eRanges(object)), "\n"))
    cat(paste0("    p juncs: ", length(jRanges(object)), "\n"))
}

#' @rdname concomp-class
setMethod("show",
          signature(object = "concomp"),
          function(object) .show.concomp(object)) 
