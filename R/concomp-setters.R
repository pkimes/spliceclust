
#' @describeIn concomp set \code{eRanges} slot
#' @aliases eRanges<-,concomp-method
setReplaceMethod("eRanges", "concomp",
                 function(obj, value) { obj@eRanges <- value; obj })


#' @describeIn concomp set \code{jRanges} slot
#' @aliases jRanges<-,concomp-method
setReplaceMethod("jRanges", "concomp",
                 function(obj, value) { obj@jRanges <- value; obj })


#' @describeIn concomp set \code{eCoverage} slot
#' @aliases eCoverage<-,concomp-method
setReplaceMethod("eCoverage", "concomp",
                 function(obj, value) { obj@eCoverage <- value; obj })


#' @describeIn concomp set \code{jCoverage} slot
#' @aliases jCoverage<-,concomp-method
setReplaceMethod("jCoverage", "concomp",
                 function(obj, value) { obj@jCoverage <- value; obj })



