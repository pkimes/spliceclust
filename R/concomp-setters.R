#' @describeIn concomp set \code{eRanges} slot
setReplaceMethod("eRanges", "concomp",
                 function(object, value) { object@eRanges <- value; object })


#' @describeIn concomp set \code{jRanges} slot
setReplaceMethod("jRanges", "concomp",
                 function(object, value) { object@jRanges <- value; object })


#' @describeIn concomp set \code{eCoverage} slot
setReplaceMethod("eCoverage", "concomp",
                 function(object, value) { object@eCoverage <- value; object })


#' @describeIn concomp set \code{jCoverage} slot
setReplaceMethod("jCoverage", "concomp",
                 function(object, value) { object@jCoverage <- value; object })



