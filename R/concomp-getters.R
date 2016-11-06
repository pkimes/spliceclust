#' @describeIn concomp get \code{eRanges} slot
setMethod("eRanges", "concomp", function(object) return(object@eRanges))


#' @describeIn concomp get \code{jRanges} slot
setMethod("jRanges", "concomp", function(object) return(object@jRanges))


#' @describeIn concomp get \code{eCoverage} slot
setMethod("eCoverage", "concomp", function(object) return(object@eCoverage))


#' @describeIn concomp get \code{jCoverage} slot
setMethod("jCoverage", "concomp", function(object) return(object@jCoverage))



