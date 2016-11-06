#' SpliceGraHM2: two-class splice graph heatmap
#'
#' Splice graphs as disjoint heatmaps for the comparative visualization
#' of splice graphs across multiple samples in a single figure. Differs from
#' \code{splicegrahm} by allowing 2 graphs to be drawn in parallel along the
#' horizontal axis flipped along the vertical axis.
#'
#' @param obj1 a \code{concomp} object
#' @param obj2 a \code{concomp} object
#' @param sort_sep a logical whether to sort each exon, junction separately
#'        (default = FALSE)
#' @param sort_idx1 an integer value specifying the order of the obj1 samples in
#'        each exon, see details for more information on all possible
#'        input, if length is n, then this ordering is used (default = 1)
#' @param sort_idx2 an integer value specifying the order of the obj2 samples in
#'        each exon, see details for more information on all possible
#'        input, if length is n, then this ordering is used (default = 1)
#' @param log_base a numeric specifying the scale of the binning for the
#'        plotting of expression values at each exon, which 0 resulting in no long
#'        scaling being applied (default = 10)
#' @param log_shift a numeric specifying the shift to be used in the log transformation
#'        for handling 0 values (default = 1)
#' @param bin a logical whether to bin the values for easier viewing (default = TRUE)
#' @param genomic a logical whether genomic coordinates should be used to
#'        plot the heatmap (default = TRUE)
#' @param ex_use a numeric specifying the proportion of the plot exons should occupy if
#'        non-genomic coordinate plotting is desired (default = 2/3)
#' @param flip_neg a logical whether to flip plotting of genes on negative strand
#'        to run left to right (default = TRUE)
#' @param j_incl a logical whether to include heatmaps for junctions
#'        (default = FALSE)
#' @param use_blk a logical whether to use a black background (default = FALSE)
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
#' @param title a character string title printed at the top of plot (default = "")
#' @param mirror a logical whether bottom plot should be flipped along vertical
#'        axis for 'mirrored' effect (default = TRUE)
#' @param same_scale a logical whether top and bottom plots should be shown on same
#'        vertical scaling, i.e. to highlight group size differences (default = TRUE)
#' @param ... other parameters to be passed
#' 
#' @return
#' a ggplot2 plot showing the splice graph with heatmaps shown at each node
#' and each splice
#'
#' @details
#' sort_idx can take values of either:
#' \itemize{
#' \item{\code{1}}: sort on first exon
#' \item{\code{2}}: sort on PC 2
#' \item{\code{3}}: sort on mean exon coverage
#' \item{\code{4}}: sort on mean exon log-coverage
#' \item{\code{5}}: rev sort on mean exon coverage
#' \item{\code{6}}: rev sort on mean exon log-coverage
#' }
#' 
#' @name splicegrahm2
#' @import ggplot2 GenomicRanges
#' @importFrom ggbio ggbio geom_alignment autoplot tracks
#' @importFrom grid arrow unit
#' @importFrom reshape2 melt
#' @author Patrick Kimes
NULL
