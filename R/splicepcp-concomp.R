.splicepcp.concomp <- function(obj, log_base, log_shift, genomic, ex_use,
                               flip_neg, imodel, highlight, eps,
                               txlist, txdb, orgdb, ...) {
    
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


#' @rdname splicepcp
setMethod("splicepcp",
          signature(obj = "concomp"),
          .splicepcp.concomp)
