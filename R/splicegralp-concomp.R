.splicegralp.concomp <- function(obj, e_loads, j_loads = NULL, load_lims = NULL, 
                                 genomic = TRUE, ex_use = 2/3,
                                 flip_neg = TRUE, use_blk = FALSE,
                                 txlist = NULL, txdb = NULL,
                                 orgdb = NULL, ...) {

    ##unpack concomp
    gr_e <- eRanges(obj)
    gr_j <- jRanges(obj)

    
    ##dataset dimension
    p_e <- length(gr_e)
    p_j <- length(gr_j)
    dna_len <- width(range(gr_e))
    rna_len <- sum(width(gr_e))

    n_vec <- ncol(e_loads)
    colnames(e_loads) <- paste0("Dir", 1:n_vec)
    colnames(j_loads) <- paste0("Dir", 1:n_vec)


    ##check sizes of loading matrices
    if (nrow(e_loads) != p_e ||
        nrow(j_loads) != p_j) {
        stop("loading dimensions not same as obj dimensions")
    }
    

    ##determine overlapping annotations
    if (is.null(txlist)) {
        tx_plot <- NULL
    } else {
        tx_plot <- find_annotations(obj, txlist, txdb, orgdb, eps = NULL)
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
    if (all(strand(gr_e) == "*") && !is.null(txlist) && length(tx_plot) > 0) {
        iflip <- flip_neg && all(strand(tx_plot) == '-')
    } else {
        iflip <- flip_neg && all(strand(gr_e) == '-')
    }

    
    ##create dataframe for plotting
    sl_df <- sl_create(gr_e, e_loads)


    ##create base plot
    sl_obj <- sl_drawbase(sl_df, gr_e, gr_j, p_e, p_j, e_loads, j_loads, 
                          genomic, use_blk, load_lims, n_vec)

    
    ##plot with horizontal axis oriented on negative strand
    if (iflip) { sl_obj <- sl_obj + scale_x_reverse() }


    ##add annotations if txdb was passed to function
    if (!is.null(txlist) && length(tx_plot) > 0) {
        if (iflip) { annot_track <- annot_track + scale_x_reverse() }
        sl_obj <- tracks(sl_obj, annot_track, heights=c(2, 1))
    }

    sl_obj
}


#' @rdname splicegralp
setMethod("splicegralp",
          signature(obj = "concomp"),
          function(obj, ...) .splicegralp.concomp(obj, ...))
