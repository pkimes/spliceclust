#' @describeIn concomp get \code{eRanges} slot
setMethod("eRanges", "concomp", function(obj) return(obj@eRanges))


#' @describeIn concomp get \code{jRanges} slot
setMethod("jRanges", "concomp", function(obj) return(obj@jRanges))


#' @describeIn concomp get \code{eCoverage} slot
setMethod("eCoverage", "concomp", function(obj) return(obj@eCoverage))


#' @describeIn concomp get \code{jCoverage} slot
setMethod("jCoverage", "concomp", function(obj) return(obj@jCoverage))



