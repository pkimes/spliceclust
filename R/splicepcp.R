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
#' @export
#' @import ggplot2 GenomicRanges
#' @importFrom ggbio geom_alignment
#' @importFrom grid arrow unit
#' @importFrom reshape2 melt
#' @author Patrick Kimes
NULL

.splicepcp.concomp <- function(obj, log_base = 10, log_shift = 1, genomic = TRUE,
                               ex_use = 2/3, flip_neg = TRUE, imodel = TRUE, 
                               highlight = NULL, eps = 1e4,
                               txlist = NULL, txdb = NULL, orgdb = NULL, ...) {
    
    ## eCoverage and jCoverage must be specified
    if (is.null(eCoverage(obj)) || is.null(jCoverage(obj)))
        stop(paste0("eCoverage and jCoverage cannot be NULL for splicepcp, \n",
                    "consider using splicegralp instead."))

    ##unpack concomp
    gr_e <- eRanges(obj)
    gr_j <- jRanges(obj)
    vals_e <- eCoverage(obj)
    vals_j <- jCoverage(obj)
    
    ##dataset dimension
    n <- ncol(vals_e)
    p_e <- nrow(vals_e)
    p_j <- nrow(vals_j)
    dna_len <- width(range(gr_e))
    rna_len <- sum(width(gr_e))

    ## highlight must be vector of length n
    if (!is.null(highlight)) {
        if (length(highlight) != n) {
            stop("highlight must be a vector of length n.")
        }
        if (length(unique(highlight)) > 9) {
            stop("highlight currently only support up to 9 unique groups.")
        }
    }
    
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
        genomic <- adj_out$genomic
        if (genomic) {
            warning("Using 'genomic = TRUE' since exons occupy more than specified 'ex_use' ",
                    "proportion of plot in genomic coordinates. No need to squish.")
        }
    }

    
    ##determine whether plots should be flipped
    if (all(strand(gr_e) == "*") && !is.null(txlist) && !is.null(tx_plot)) {
        iflip <- flip_neg && all(strand(tx_plot) == '-')
    } else {
        iflip <- flip_neg && all(strand(gr_e) == '-')
    }


    ##construct data.frame
    sp_df <- sp_create(gr_e, gr_j, vals_e, vals_j,
                       log_base, log_shift, n, p_e, p_j)
    
    ## sp_drawbase -- same as sg_drawbase?
    sp_obj <- sp_drawbase(sp_df, highlight, genomic,
                          log_base, gr_e, p_e)

    
    ##combine gene model from splicegrahm or splicegralp to plot
    if (imodel) {
        sp_mobj <- sp_drawmodel(sp_df, gr_e, gr_j, p_j, genomic)
        if (iflip) { sp_mobj <- sp_mobj + scale_x_reverse() }
    }

    
    ##plot with horizontal axis oriented on negative strand
    if (iflip) { sp_obj <- sp_obj + scale_x_reverse() }

    if (iflip && !is.null(tx_plot)) {
        annot_track <- annot_track + scale_x_reverse()
    }

    
    iannots <- !(is.null(txlist) || is.null(tx_plot))
    
    ##combine generated tracks
    if (imodel && iannots) {
        sp_obj <- tracks(sp_obj, sp_mobj, annot_track, heights=c(3, 1, 1.5))
    } else if (imodel) {
        sp_obj <- tracks(sp_obj, sp_mobj, heights=c(3, 1))
    } else if (iannots) {
        sp_obj <- tracks(sp_obj, annot_track, heights=c(2, 1))
    }
    
    sp_obj
}

#' @keywords internal
#' @title splicepcp method
#' @rdname splicepcp
#' @aliases splicepcp,concomp-method
setMethod("splicepcp",
          signature(obj = "concomp"),
          function(obj, ...) .splicepcp.concomp(obj, ...))
