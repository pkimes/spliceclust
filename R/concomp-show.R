
.show.concomp <- function(object) {
    cat("\n concomp object summary: \n")
    cat(paste0("    n sampl: ", ncol(eCoverage(object)), "\n"))
    cat(paste0("    p exons: ", length(eRanges(object)), "\n"))
    cat(paste0("    p juncs: ", length(jRanges(object)), "\n"))
}


#' show summary of concomp object
#'
#' show quick summary of \code{concomp} object
#' 
#' @param object a \code{concomp} object
#' 
#' @return Descrption of \code{concomp} object
#'
#' @keywords internal
#' @name show-concomp
#' @aliases show,concomp-method
#' @author Patrick Kimes
setMethod("show",
          signature(object = "concomp"),
          function(object) .show.concomp(object)) 
