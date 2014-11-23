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
#' @import ggplot2 RColorBrewer reshape2 grid
#' @author Patrick Kimes
NULL

.splicegralp.concomp <- function(obj, e_loads, j_loads = NULL, load_lims = NULL, 
                                 genomic = TRUE, ex_use = 2/3,
                                 flip_neg = TRUE, use_blk = FALSE,
                                 txlist = NULL, txdb = NULL,
                                 orgdb = NULL, ...) {
    
    ##currently can't include gene models if not plotting on genomic scale
    if (!is.null(txlist) && !genomic)
        cat("ignoring txlist input since plotting on non-genomic scale. \n")

    
    ##unpack concomp
    gr_e <- exons(obj)
    gr_j <- juncs(obj)
    
    ##dataset dimension
    p_e <- length(gr_e)
    p_j <- length(gr_j)
    n_vec <- ncol(e_loads)
    colnames(e_loads) <- paste0("Dir", 1:n_vec)
    colnames(j_loads) <- paste0("Dir", 1:n_vec)
    

    ##check sizes of loading matrices
    if (nrow(e_loads) != p_e ||
        nrow(j_loads) != p_j) {
        stop("loading dimensions not same as obj dimensions")
    }

    
    ##color generators
    crp <- colorRampPalette(c("#053061", "#f7f7f7", "#67001f"))
    crp2 <- colorRamp(c("#053061", "#f7f7f7", "#67001f"))


    ##change GRanges coordinates if non-genomic coordinates are desired
    if (!genomic) {
        dna_len <- width(range(gr_e))
        rna_len <- sum(width(gr_e))

        shift_1 <- -start(ranges(gr_e))[1]+1
        ranges(gr_e) <- shift(ranges(gr_e), shift_1)
        ranges(gr_j) <- shift(ranges(gr_j), shift_1)

        shrink_by <- 1 - (1 - ex_use)/ex_use * rna_len/(dna_len - rna_len)
        gps <- distance(ranges(gr_e)[-p_e], ranges(gr_e)[-1]) * shrink_by

        eeb <- start(ranges(gr_e))[-1]
        for (ii in (p_e-1):1) {
            i_adj <- start(ranges(gr_e)) >= eeb[ii]
            start(ranges(gr_e))[i_adj] <- start(ranges(gr_e))[i_adj] - gps[ii]
            
            i_adj <- end(ranges(gr_e)) >= eeb[ii]
            end(ranges(gr_e))[i_adj] <- end(ranges(gr_e))[i_adj] - gps[ii]

            i_adj <- start(ranges(gr_j)) >= eeb[ii]
            start(ranges(gr_j))[i_adj] <- start(ranges(gr_j))[i_adj] - gps[ii]
            
            i_adj <- end(ranges(gr_j)) >= eeb[ii]
            end(ranges(gr_j))[i_adj] <- end(ranges(gr_j))[i_adj] - gps[ii]
        }
    }


    ##determine whether plots should be flipped
    iflip <- flip_neg && all(strand(gr_e) == '-')

    
    ##strand of junctions for arrow heads
    arrowhead <- ifelse(as.character(strand(gr_j)) == "-",
                        "last", "first")

    
    ##construct ggplot2 object
    gg_e <- data.frame(xmin=start(ranges(gr_e)),
                       xmax=end(ranges(gr_e)),
                       e_loads)
    gg_e <- reshape2::melt(gg_e, id.vars=c("xmin", "xmax"))
    gg_e$ymin <- -1
    gg_e$ymax <- 1
    
    
    ##plot on genomic coordinates
    g_obj <- ggplot(gg_e, aes(xmin=xmin, xmax=xmax,
                              ymin=ymin, ymax=ymax,
                              color=value, fill=value))


    ##add horizontal line first
    g_obj <- g_obj + 
        geom_hline(yintercept=0,
                   color=ifelse(use_blk, "#F0F0F0", "#3C3C3C"))

    
    ##add basic plot structure
    g_obj <-
        g_obj + 
        geom_rect() +
        facet_grid(variable~.) +
        scale_y_continuous(breaks=NULL, limits=c(-1, 3)) +
        ylab("") +
        xlab(ifelse(genomic, paste0("Genomic Coordinates, ", seqnames(gr_e[1])),
                    "non-genomic coordinates")) +
        theme_bw()
    

    ##frame exons if not using black background
    if (use_blk) {
        g_obj <- g_obj +
            theme(panel.grid.minor.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.background = element_rect(fill="#3C3C3C"))
    } else {
        g_obj <- g_obj +
            annotate("rect", size=.25,
                     xmin=start(ranges(gr_e)) - 1,
                     xmax=end(ranges(gr_e)) + 1,
                     ymin=-1, ymax=1,
                     alpha=1, color="#3C3C3C", fill=NA)
    }


    ##determine the same coloring bounds for splices and exons
    if (is.null(load_lims)) {
        color_lim <- max(max(abs(j_loads)), max(abs(e_loads)))
        color_lim <- c(-color_lim, color_lim)
    } else {
        color_lim <- load_lims
    }

    
    ##add continuous color palette
    g_obj <- g_obj +
        scale_color_gradient2("expr", limits=color_lim,
                              low="#053061", mid="#f7f7f7", high="#67001f",
                              guide="none") +
        scale_fill_gradient2("expr", limits=color_lim,
                             low="#053061", mid="#f7f7f7", high="#67001f")


    ##add splicing arrows to plot
    width_props <- width(ranges(gr_j)) / width(range(gr_e))
    j_col_scale <- j_loads / diff(color_lim) + .5
    for (j in 1:p_j) {
        w_prop <- width_props[j]
        for (ipc in 1:n_vec) {
            circle1 <- .pseudoArc(xmin=start(ranges(gr_j))[j],
                                  xmax=end(ranges(gr_j))[j],
                                  ymin=1, height=2*sqrt(w_prop))
            circle2 <- cbind(circle1, xmin=gg_e$xmin[1], xmax=gg_e$xmax[1],
                             ymin=gg_e$xmin[1], ymax=gg_e$ymax[1],
                             variable=paste0("Dir", ipc), value=0)
            g_obj <- g_obj +
                geom_path(data=circle2, size=.75, aes(x=x, y=y),
                          color=.rgb2hex(crp2(j_col_scale[j, ipc])),
                          arrow=arrow(length=unit(.015*n_vec, "npc"),
                              ends=arrowhead[j]))
        }
    }


    ##plot with horizontal axis oriented on negative strand
    if (iflip) {
        g_obj <- g_obj + scale_x_reverse()
    }


    ##add annotations if txdb was passed to function
    if (!is.null(txlist) && genomic) {
        cand_idx <- unique(queryHits(findOverlaps(txlist, gr_e)))

        if (is.null(txdb)) {
            cand_names <- as.character(cand_idx)
        } else {
            cand_names <- concomp2name(obj, txlist, txdb, orgdb)
        }
        
        if (length(cand_idx) > 0) {
            tx_match <- txlist[cand_idx]
            names(tx_match) <- make.unique(cand_names)
            annot_track <- ggplot(tx_match) +
                geom_alignment() + theme_bw()            
            if (iflip) { annot_track <- annot_track + scale_x_reverse() }

            g_obj <- tracks(g_obj, annot_track, heights=c(2, 1))
        }
    }
    
    g_obj
}



#' @rdname splicegralp
#' @aliases splicegralp,concomp-method
setMethod("splicegralp",
          signature(obj = "concomp"),
          function(obj, ...) .splicegralp.concomp(obj, ...))
