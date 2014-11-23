
#' @describeIn concomp get \code{exons} slot
#' @aliases exons,concomp-method
setMethod("exons", "concomp", function(obj) return(obj@exons))


#' @describeIn concomp get \code{juncs} slot
#' @aliases juncs,concomp-method
setMethod("juncs", "concomp", function(obj) return(obj@juncs))


#' @describeIn concomp get \code{exonValues} slot
#' @aliases exonValues,concomp-method
setMethod("exonValues", "concomp", function(obj) return(obj@exonValues))


#' @describeIn concomp get \code{juncValues} slot
#' @aliases juncValues,concomp-method
setMethod("juncValues", "concomp", function(obj) return(obj@juncValues))



