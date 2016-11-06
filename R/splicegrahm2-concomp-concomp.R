.splicegrahm2.concomp.concomp <- function(obj1, obj2, sort_sep, sort_idx1, sort_idx2,
                                          log_base, log_shift, bin, genomic, ex_use,
                                          flip_neg, j_incl, use_blk, eps, txlist, txdb,
                                          orgdb, title, mirror, same_scale, ...) {

    ## eCoverage and jCoverage must be specified
    if (is.null(eCoverage(obj1)) || is.null(jCoverage(obj1)) ||
        is.null(eCoverage(obj2)) || is.null(jCoverage(obj2)))
        stop(paste0("eCoverage and jCoverage cannot be NULL for splicegrahm2, \n",
                    "consider using splicegralp instead."))


    ##unpack and compute dimensions for concomp 1
    gr_e1 <- eRanges(obj1)
    gr_j1 <- jRanges(obj1)
    vals_e1 <- eCoverage(obj1)
    vals_j1 <- jCoverage(obj1)
    n1 <- ncol(vals_e1)
    p_e1 <- nrow(vals_e1)
    p_j1 <- nrow(vals_j1)
    

    ##unpack and compute dimensions for concomp 2
    gr_e2 <- eRanges(obj2)
    gr_j2 <- jRanges(obj2)
    vals_e2 <- eCoverage(obj2)
    vals_j2 <- jCoverage(obj2)
    n2 <- ncol(vals_e2)
    p_e2 <- nrow(vals_e2)
    p_j2 <- nrow(vals_j2)

    ##take maximum of both datasets to get dimensions of plot
    n_max <- max(n1, n2)
    
    ##create single concomp with both ranges
    gr_e <- c(eRanges(obj1), eRanges(obj2))
    gr_j <- c(jRanges(obj1), jRanges(obj2))
    strand(gr_e) <- "*"
    strand(gr_j) <- "*"
    obj <- concomp(GRangesList("e"=gr_e, "j"=gr_j))

    ##determine overlapping annotations for concomps 1 and 2
    if (is.null(txlist)) {
        tx_plot <- NULL
    } else {
        tx_plot <- find_annotations(obj, txlist, txdb, orgdb, eps)
    }

    if (max(p_e1, p_e2) < 2) { genomic <- TRUE }

    ##change GRanges coordinates if non-genomic coordinates are desired
    if (genomic) {
        if (is.null(tx_plot)) {
            annot_track <- NULL
        } else { 
            annot_track <- ggbio() +
                geom_alignment(tx_plot, gap.geom="arrow", aes(group=tx)) +
                    theme_bw()
        }
    } else {
        adj_out <- adj_ranges(gr_e1, gr_j1, tx_plot, ex_use, eRanges(obj))
        gr_e1 <- adj_out$gr_e
        gr_j1 <- adj_out$gr_j

        adj_out <- adj_ranges(gr_e2, gr_j2, tx_plot, ex_use, eRanges(obj))
        gr_e2 <- adj_out$gr_e
        gr_j2 <- adj_out$gr_j
        annot_track <- adj_out$annot_track
        genomic <- adj_out$genomic
        if (genomic) {
            warning("Using 'genomic = TRUE' since exons occupy more than specified 'ex_use' ",
                    "proportion of plot in genomic coordinates. No need to squish.")
        }
    }
    

    ##determine whether plots should be flipped
    if (all(strand(gr_e1) == "*") && all(strand(gr_e2) == "*") &&
        !is.null(txlist) && !is.null(tx_plot)) {
        iflip <- flip_neg && all(strand(tx_plot) == '-')
    } else {
        iflip <- flip_neg && all(strand(gr_e1) == '-') && all(strand(gr_e2) == '-')
    }


    ##determine order of samples
    if (sort_sep) {
        vals_e1 <- t(apply(vals_e1, 1, sort))
        vals_j1 <- t(apply(vals_j1, 1, sort))
    } else {
        idx <- sampl_sort(sort_idx1, vals_e1, vals_j1, n1)
        vals_e1 <- vals_e1[, idx, drop=FALSE]
        vals_j1 <- vals_j1[, idx, drop=FALSE]

        idx <- sampl_sort(sort_idx2, vals_e2, vals_j2, n2)
        vals_e2 <- vals_e2[, idx, drop=FALSE]
        vals_j2 <- vals_j2[, idx, drop=FALSE]
    }
    
    
    ##create dataframe for plotting
    sg_df1 <- sg_create(gr_e1, gr_j1, vals_e1, vals_j1, j_incl,
                        log_base, log_shift, bin, n1, p_j1, FALSE,
                        ifelse(same_scale, n_max, n1))
    sg_df2 <- sg_create(gr_e2, gr_j2, vals_e2, vals_j2, j_incl,
                        log_base, log_shift, bin, n2, p_j2, mirror,
                        ifelse(same_scale, n_max, n2))
    if (bin) {
        v_max <- max(length(levels(sg_df1$value)), length(levels(sg_df2$value))) - 1
        levels(sg_df1$value) <- 0:v_max
        levels(sg_df2$value) <- 0:v_max
    }

    
    ##plot on genomic coordinates
    sg_obj1 <- sg_drawbase(sg_df1, use_blk, j_incl, genomic,
                           gr_e1, log_base, bin, n1, NULL, p_j1, iflip, FALSE,
                           ifelse(same_scale, n_max, n1))
    sg_obj2 <- sg_drawbase(sg_df2, use_blk, j_incl, genomic,
                           gr_e2, log_base, bin, n2, NULL, p_j2, iflip, mirror,
                           ifelse(same_scale, n_max, n2))
    
    
    ##add arrow information if needed
    if (p_j1 > 0) {
        sg_obj1 <- sg_drawjuncs(sg_obj1, sg_df1, j_incl, use_blk, iflip,
                                gr_e1, gr_j1, vals_j1, n1, p_j1, NULL, FALSE,
                                ifelse(same_scale, n_max, n1))
    }
    if (p_j2 > 0) {
        sg_obj2 <- sg_drawjuncs(sg_obj2, sg_df2, j_incl, use_blk, iflip,
                                gr_e2, gr_j2, vals_j2, n2, p_j2, NULL, mirror,
                                ifelse(same_scale, n_max, n2))
    }

    
    ##plot with horizontal axis oriented on negative strand
    if (iflip) {
        sg_obj1 <- sg_obj1 + scale_x_reverse()
        sg_obj2 <- sg_obj2 + scale_x_reverse()
    }
    
    
    ##add annotations if txdb was passed to function
    if (!is.null(txlist) && !is.null(tx_plot)) {
        if (iflip) { annot_track <- annot_track + scale_x_reverse() }
        sg_obj <- tracks(sg_obj1, annot_track, sg_obj2, heights=c(2, 1, 2), title=title)
    } else {
        sg_obj <- tracks(sg_obj1, sg_obj2, heights=c(2, 2), title=title)
    }

    sg_obj
}


#' @rdname splicegrahm2
setMethod("splicegrahm2",
          signature(obj1 = "concomp", obj2 = "concomp"),
          .splicegrahm2.concomp.concomp)
