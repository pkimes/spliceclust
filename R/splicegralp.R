#' SpliceGraLP: splice graph loadings plot
#'
#' Splice graphs with loadings shown along the exons and junctions to illustrate
#' expression patterns across a connected component. The function can be used to
#' illustrate, e.g. differences across groups/clusters or principal component
#' loadings.
#'
#' @param obj a \code{concomp} object with exon and junction information
#' @param e_loads a matrix of exon loadings in columns
#' @param j_loads a matrix of junction loadings in columns (default = NULL)
#' @param load_lims a vector of length 2 that should be used as the bounds on
#'        the color scale, e.g. \code{c(-1, 1)} if plotting PCA loadings 
#' @param genomic a logical whether genomic coordinates should be used to
#'        plot the heatmap (default = TRUE)
#' @param ex_use a numeric specifying the proportion of the plot exons should occupy if
#'        non-genomic coordinate plotting is desired (default = 2/3)
#' @param flip_neg a logical whether to flip plotting of genes on negative strand
#'        to run left to right (default = TRUE)
#' @param use_blk a logical whether to use a black background (default = FALSE)
#' @param txlist a GRangesList of transcripts or genes which should be queried and
#'        added to the plot if falling within the region of the connected component
#'        (default = NULL)
#' @param txdb a transcript database which can be used to query the transcript IDs
#'        identified from txlist (default = NULL)
#' @param orgdb a database that can be queried using keys obtained from \code{txdb}
#'        to determine corresponding gene symbols (default = NULL)
#' @param ... other parameters to be passed
#' 
#' @return
#' a ggplot2 plot showing the specified number of loadings
#'
#' @name splicegralp
#' @import ggplot2
#' @importFrom ggbio geom_alignment
#' @importFrom grid arrow unit
#' @importFrom reshape2 melt
#' @author Patrick Kimes
NULL
