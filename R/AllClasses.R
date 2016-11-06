## ###################################################################
## ###################################################################
## construct union class to allow for empty eCoverage jCoverage slots

setOldClass("matrix")
setClassUnion("MatrixOrNULL", c("matrix", "NULL"))


## ###################################################################
## ###################################################################
## main connected component class

#' @name concomp-class
setClass("concomp",
         slots=list(
             eRanges = "GRanges",
             jRanges = "GRanges",
             eCoverage = "MatrixOrNULL",
             jCoverage = "MatrixOrNULL"))

