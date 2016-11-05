
#' @describeIn concomp get \code{eRanges} slot
#' @aliases eRanges,concomp-method
setMethod("eRanges", "concomp", function(obj) return(obj@eRanges))


#' @describeIn concomp get \code{jRanges} slot
#' @aliases jRanges,concomp-method
setMethod("jRanges", "concomp", function(obj) return(obj@jRanges))


#' @describeIn concomp get \code{eCoverage} slot
#' @aliases eCoverage,concomp-method
setMethod("eCoverage", "concomp", function(obj) return(obj@eCoverage))


#' @describeIn concomp get \code{jCoverage} slot
#' @aliases jCoverage,concomp-method
setMethod("jCoverage", "concomp", function(obj) return(obj@jCoverage))



