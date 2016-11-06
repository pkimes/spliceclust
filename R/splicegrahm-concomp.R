.splicegrahm.concomp <- function(obj, sort_sep = FALSE, sort_idx = 1,
                                 log_base = 10, log_shift = 1, bin = TRUE,
                                 genomic = TRUE, ex_use = 2/3, flip_neg = TRUE, 
                                 j_incl = FALSE, highlight = NULL,
                                 use_blk = FALSE, eps = 1e4, txlist = NULL,
                                 txdb = NULL, orgdb = NULL, title="", ...) {

    ## eCoverage and jCoverage must be specified
    if (is.null(eCoverage(obj)) || is.null(jCoverage(obj)))
        stop(paste0("eCoverage and jCoverage cannot be NULL for splicegrahm, \n",
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
    
    ##determine overlapping annotations
    if (is.null(txlist)) {
        tx_plot <- NULL
    } else {
        tx_plot <- find_annotations(obj, txlist, txdb, orgdb, eps)
    }

    if (p_e < 2) { genomic <- TRUE }

    ##change GRanges coordinates if non-genomic coordinates are desired
    if (genomic) {
        if (is.null(tx_plot)) {
            annot_track <- NULL
        } else { 
            annot_track <- ggbio() +
                geom_alignment(tx_plot, gap.geom="arrow", aes(group=tx)) +
                    theme_bw()
        }
    } else if (p_e > 1) {
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
        vals_e <- vals_e[, idx, drop=FALSE]
        vals_j <- vals_j[, idx, drop=FALSE]
    }

    
    ##create dataframe for plotting
    sg_df <- sg_create(gr_e, gr_j, vals_e, vals_j, j_incl,
                      log_base, log_shift, bin, n, p_j)
    if (bin) {
        v_max <- length(levels(sg_df$value)) - 1
        levels(sg_df$value) <- 0:v_max
    }
    
    ##plot on genomic coordinates
    sg_obj <- sg_drawbase(sg_df, use_blk, j_incl, genomic,
                          gr_e, log_base, bin, n, highlight, p_j, iflip)
    
    
    ##add arrow information if needed
    if (p_j > 0) {
        sg_obj <- sg_drawjuncs(sg_obj, sg_df, j_incl, use_blk, iflip,
                               gr_e, gr_j, vals_j, n, p_j, highlight)
    }
    
    ##plot with horizontal axis oriented on negative strand
    if (iflip) { sg_obj <- sg_obj + scale_x_reverse() }


    ##add annotations if txdb was passed to function
    if (!is.null(txlist) && !is.null(tx_plot)) {
        if (iflip) { annot_track <- annot_track + scale_x_reverse() }
        sg_obj <- tracks(sg_obj, annot_track, heights=c(2, 1), title=title)
    } else {
        sg_obj <- sg_obj + ggtitle(title)
    }

    sg_obj
}

#' @rdname splicegrahm
setMethod("splicegrahm",
          signature(obj = "concomp"),
          function(obj, ... ) .splicegrahm.concomp(obj, ...))



