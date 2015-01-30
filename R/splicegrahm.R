#' SpliceGraHM: splice graph heatmap
#'
#' Splice graphs as disjoint heatmaps for the comparative visualization
#' of splice graphs across multiple samples in a single figure.
#'
#' @param obj a \code{concomp} object
#' @param sort_sep a logical whether to sort each exon, junction separately
#'        (default = FALSE)
#' @param sort_idx an integer value specifying the order of the samples in
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
#' @param highlight a vector of labels to highlight samples in groups or clusters
#'        (default = NULL)
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
#' @param ... other parameters to be passed
#' 
#' @return
#' a ggplot2 plot showing the splice graph with heatmaps shown at each node
#' and each splice
#'
#' @details
#' sort_idx can take values of either:
#' \itemize{
#' \item{\code{1}}: sort based on first exon
#' \item{\code{2}}: sort based on PC 2
#' }
#' 
#' @name splicegrahm
#' @export
#' @import ggplot2 GenomicRanges
#' @importFrom ggbio geom_alignment autoplot ggplot
#' @importFrom grid arrow unit
#' @importFrom reshape2 melt
#' @author Patrick Kimes
NULL

.splicegrahm.concomp <- function(obj, sort_sep = FALSE, sort_idx = 1,
                                 log_base = 10, log_shift = 1, bin = TRUE,
                                 genomic = TRUE, ex_use = 2/3, flip_neg = TRUE, 
                                 j_incl = FALSE, highlight = NULL,
                                 use_blk = FALSE, eps = 1e4, txlist = NULL,
                                 txdb = NULL, orgdb = NULL, ...) {

    ##exonValues and juncValues must be specified
    if (is.null(exonValues(obj)) || is.null(juncValues(obj)))
        stop(paste0("exonValues and juncValues cannot be NULL for splicegrahm, \n",
                    "consider using splicegralp instead."))

    
    ##unpack concomp
    gr_e <- exons(obj)
    gr_j <- juncs(obj)
    vals_e <- exonValues(obj)
    vals_j <- juncValues(obj)
    
    ##dataset dimension
    n <- ncol(vals_e)
    p_e <- nrow(vals_e)
    p_j <- nrow(vals_j)
    dna_len <- width(range(gr_e))
    rna_len <- sum(width(gr_e))
    
    n <- ncol(vals_e)
    p_e <- nrow(vals_e)
    p_j <- nrow(vals_j)
    dna_len <- width(range(gr_e))
    rna_len <- sum(width(gr_e))

    ##determine overlapping annotations
    if (is.null(txlist)) {
        tx_plot <- NULL
    } else {
        tx_plot <- find_annotations(obj, txlist, txdb, orgdb, eps)
    }
    

    ##change GRanges coordinates if non-genomic coordinates are desired
    if (genomic) {
        if (is.null(tx_plot)) {
            annot_track <- NULL
        } else { 
            annot_track <- ggplot() +
                geom_alignment(tx_plot, gap.geom="arrow", aes(group=tx)) +
                    theme_bw()
        }
    } else {
        adj_out <- adj_ranges(gr_e, gr_j, tx_plot, ex_use)
        gr_e <- adj_out$gr_e
        gr_j <- adj_out$gr_j
        annot_track <- adj_out$annot_track
    }
    

    ##determine whether plots should be flipped
    if (all(strand(gr_e) == "*") && !is.null(txlist) && !is.null(tx_plot)) {
        iflip <- flip_neg && all(strand(tx_plot) == '-')
    } else {
        iflip <- flip_neg && all(strand(gr_e) == '-')
    }

    ##force highlight to be ordering for quicker plotting
    if (!is.null(highlight)) {
        sort_idx <- order(highlight)
        highlight <- sort(highlight)
    }
    
    ##determine order of samples
    if (sort_sep) {
        vals_e <- t(apply(vals_e, 1, sort))
        vals_j <- t(apply(vals_j, 1, sort))
    } else {
        idx <- sampl_sort(sort_idx, vals_e, vals_j, n)
        vals_e <- vals_e[, idx]
        vals_j <- vals_j[, idx]
    }

    
    ##create dataframe for plotting
    sg_df <- sg_create(gr_e, gr_j, vals_e, vals_j, j_incl,
                      log_base, log_shift, bin, n, p_j)

    
    ##plot on genomic coordinates
    sg_obj <- sg_drawbase(sg_df, use_blk, j_incl, genomic,
                          gr_e, log_base, bin, n, highlight, p_j)
    
    
    ##add arrow information if needed
    sg_obj <- sg_drawjuncs(sg_obj, sg_df, j_incl, use_blk, iflip,
                           gr_e, gr_j, vals_j, n, p_j, highlight)
    
    
    ##plot with horizontal axis oriented on negative strand
    if (iflip) { sg_obj <- sg_obj + scale_x_reverse() }


    ##add annotations if txdb was passed to function
    if (!is.null(txlist) && !is.null(tx_plot)) {
        if (iflip) { annot_track <- annot_track + scale_x_reverse() }
        sg_obj <- tracks(sg_obj, annot_track, heights=c(2, 1))
    }

    sg_obj
}

#' @keywords internal
#' @title splicegrahm method
#' @name splicegrahm-concomp
#' @aliases splicegrahm,concomp-method
setMethod("splicegrahm",
          signature(obj = "concomp"),
          function(obj, ... ) .splicegrahm.concomp(obj, ...))



