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
#' @export
#' @import ggplot2
#' @importFrom ggbio geom_alignment
#' @importFrom grid arrow unit
#' @importFrom reshape2 melt
#' @author Patrick Kimes
NULL

.splicegralp.concomp <- function(obj, e_loads, j_loads = NULL, load_lims = NULL, 
                                 genomic = TRUE, ex_use = 2/3,
                                 flip_neg = TRUE, use_blk = FALSE,
                                 txlist = NULL, txdb = NULL,
                                 orgdb = NULL, ...) {

    ##unpack concomp
    gr_e <- eRanges(obj)
    gr_j <- eRanges(obj)

    
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

#' splicegralp method
#' 
#' \code{splicegralp} method for \code{concomp} class object.
#' See \code{splicegralp} documentation for more details.
#'
#' @keywords internal
#' @seealso splicegralp
#' @name splicegralp-concomp
#' @aliases splicegralp,concomp-method
setMethod("splicegralp",
          signature(obj = "concomp"),
          function(obj, ...) .splicegralp.concomp(obj, ...))
