#' @describeIn concomp set \code{eRanges} slot
setReplaceMethod("eRanges", "concomp",
                 function(obj, value) { obj@eRanges <- value; obj })


#' @describeIn concomp set \code{jRanges} slot
setReplaceMethod("jRanges", "concomp",
                 function(obj, value) { obj@jRanges <- value; obj })


#' @describeIn concomp set \code{eCoverage} slot
setReplaceMethod("eCoverage", "concomp",
                 function(obj, value) { obj@eCoverage <- value; obj })


#' @describeIn concomp set \code{jCoverage} slot
setReplaceMethod("jCoverage", "concomp",
                 function(obj, value) { obj@jCoverage <- value; obj })



