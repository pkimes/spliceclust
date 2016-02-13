
#' @describeIn concomp set \code{exons} slot
#' @aliases exons<-,concomp-method
setReplaceMethod("exons", "concomp",
                 function(obj, value) { obj@exons <- value; obj })


#' @describeIn concomp set \code{juncs} slot
#' @aliases juncs<-,concomp-method
setReplaceMethod("juncs", "concomp",
                 function(obj, value) { obj@juncs <- value; obj })


#' @describeIn concomp set \code{exonValues} slot
#' @aliases exonValues<-,concomp-method
setReplaceMethod("exonValues", "concomp",
                 function(obj, value) { obj@exonValues <- value; obj })


#' @describeIn concomp set \code{juncValues} slot
#' @aliases juncValues<-,concomp-method
setReplaceMethod("juncValues", "concomp",
                 function(obj, value) { obj@juncValues <- value; obj })



