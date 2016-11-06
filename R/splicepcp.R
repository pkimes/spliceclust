#' SplicePCP: splice graph parallel coordinates plot
#'
#' Splice graph expression plot with exons along horizontal axis with vertical
#' axis showing expression level. The plot provides an alternative to the
#' heatmaps used to visualize expression using the \code{splicegrahm} and
#' \code{splicegralp} functions.
#'
#' @param obj a \code{concomp} object
#' @param log_base a numeric specifying the scaling of expression values at each exon,
#'        which 0 resulting in no long scaling being applied (default = 10)
#' @param log_shift a numeric specifying the shift to be used in the log transformation
#'        for handling 0 values (default = 1)
#' @param genomic a logical whether genomic coordinates should be used to
#'        plot the heatmap (default = TRUE)
#' @param ex_use a numeric specifying the proportion of the plot exons should occupy if
#'        non-genomic coordinate plotting is desired (default = 2/3)
#' @param flip_neg a logical whether to flip plotting of genes on negative strand
#'        to run left to right (default = TRUE)
#' @param highlight a vector of labels to highlight samples in groups or clusters
#'        (default = NULL)
#' @param imodel a logical whether to include the connected component plot generated
#'        by splicegralp or splicegrahm to the plot (default = TRUE)
#' @param eps a numeric value specifying the number of base pairs around \code{obj} to look
#'        for overlapping gene models, if eps = NULL, then all overlapping gene models are
#'        included (default = 1e4) 
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
#' a ggplot2 plot showing the splice graph with heatmaps shown at each node
#' and each splice
#'
#' @name splicepcp
#' @import ggplot2 GenomicRanges
#' @importFrom ggbio geom_alignment
#' @importFrom grid arrow unit
#' @importFrom reshape2 melt
#' @author Patrick Kimes
NULL
