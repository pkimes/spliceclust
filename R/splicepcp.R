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
#' @import ggbio ggplot2 RColorBrewer reshape2 grid GenomicRanges
#' @author Patrick Kimes
NULL

.splicepcp.concomp <- function(obj, log_base = 10, log_shift = 1, genomic = TRUE,
                               ex_use = 2/3, flip_neg = TRUE, imodel = TRUE, 
                               j_incl = FALSE, highlight = NULL, txlist = NULL,
                               txdb = NULL, orgdb = NULL, ...) {
    
    ##exonValues and juncValues must be specified
    if (is.null(exonValues(obj)) || is.null(juncValues(obj)))
        stop(paste0("exonValues and juncValues cannot be NULL for splicepcp, \n",
                    "consider using splicegralp instead."))


    ##currently can't include gene models if not plotting on genomic scale
    if (!is.null(txlist) && !genomic)
        cat("ignoring txlist input since plotting on non-genomic scale. \n")

    
    ##color generators
    crp <- colorRampPalette(c("#f7fbff", "#08306b"))
    crp2 <- colorRamp(c("#f7fbff", "#047760"))
    hl_cols <- c("#646464", "#CA4942", "#5A3589")

    ##long list of letters
    ll_LETTERS <- LETTERS
    for (k in 2:5) {
        ll_LETTERS <- c(ll_LETTERS,
                        do.call("paste0", rep(list(LETTERS), k)))
    }

    
    ##unpack concomp
    gr_e <- exons(obj)
    gr_j <- juncs(obj)
    vals_e <- exonValues(obj)
    vals_j <- juncValues(obj)
    
    ##dataset dimension
    n <- ncol(vals_e)
    p_e <- nrow(vals_e)
    p_j <- nrow(vals_j)

    
    ##change GRanges coordinates if non-genomic coordinates are desired
    if (!genomic) {
        dna_len <- width(range(gr_e))
        rna_len <- sum(width(gr_e))

        ##don't try to shrink gene model if intronic space is already small
        if (rna_len/dna_len <= ex_use) {
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
            
        } else {
            genomic <- TRUE
        }
        
    }


    ##determine whether plots should be flipped
    iflip <- flip_neg && all(strand(gr_e) == '-')

    ##if negative strand ordering, flip order of obj
    if (iflip) {
        gr_e <- rev(gr_e)
        gr_j <- rev(gr_j)
        vals_e <- vals_e[p_e:1, , drop=FALSE]
        vals_j <- vals_j[p_j:1, , drop=FALSE]
    }

    
    ##strand of junctions for arrow heads
    arrowhead <- ifelse(as.character(strand(gr_j)) == "-",
                        "last", "first")
    

    ##construct ggplot2 object
    gg_e <- data.frame(xmin=start(ranges(gr_e)),
                       xmax=end(ranges(gr_e)),
                       vals_e)
    gg_e <- reshape2::melt(gg_e, id.vars=c("xmin", "xmax"))

    ##change y values 
    gg_e$ymin <- gg_e$value
    gg_e$ymax <- gg_e$value

    
    ##transform if desired
    if (log_base > 0) {
        gg_e$ymin <- log(log_shift + gg_e$ymin, base=log_base)
        gg_e$ymax <- log(log_shift + gg_e$ymax, base=log_base)
        v_max <- max(gg_e$ymax)
    }
    
    
    ##color by groups if specified for highlight
    if (is.null(highlight)) {
        g_obj <- ggplot(gg_e, aes(x=xmin, xend=xmax,
                                  y=ymin, yend=ymax))
    } else {
        gg_e$hl <- rep(as.factor(highlight), each=p_e)
        g_obj <- ggplot(gg_e, aes(x=xmin, xend=xmax,
                                  y=ymin, yend=ymax, color=hl))
    }        

    ##add horizontal line first
    g_obj <- g_obj + 
        geom_hline(yintercept=0, color="#3C3C3C")

    
    ##add main segments
    g_obj <- g_obj + geom_segment(alpha=2/3, size=1/3)

    
    ##add basic plot structure
    g_obj <- g_obj + 
        ylab(ifelse(log_base==0, "Expression", "log Expression")) +
        xlab(ifelse(genomic, paste0("Genomic Coordinates, ", seqnames(gr_e[1])),
                    "non-genomic coordinates")) +
        theme_bw()

    
    ## ##add discrete color palette
    ## if (!is.null(highlight)) {
    ##     g_obj <- g_obj +
    ##         scale_color_manual("expr", values=hl_cols[1:nlevels(gg_e$hl)])
    ## }
    

    ##plot with horizontal axis oriented on negative strand
    if (iflip) {
        g_obj <- g_obj + scale_x_reverse()
    }


    ##warm user that txlist won't be processed unless on genomic coordinates
    if (!is.null(txlist) && !genomic) {
        cat(paste0("!!!  can't plot model transcripts unless ",
                   "genomic=TRUE  !!!\n"))
    }

    
    ##combine gene model from splicegrahm or splicegralp to plot
    if (imodel) {
        ##construct ggplot2 object
        gg_mod <- data.frame(xmin=start(ranges(gr_e)),
                             xmax=end(ranges(gr_e)))
        gg_mod <- reshape2::melt(gg_mod, id.vars=c("xmin", "xmax"))
        gg_mod$ymin <- -1
        gg_mod$ymax <- 1
        
        ##plot on genomic coordinates
        g_mobj <- ggplot(gg_mod, aes(xmin=xmin, xmax=xmax,
                                     ymin=ymin, ymax=ymax))
        
        ##add horizontal line first
        g_mobj <- g_mobj + 
            geom_hline(yintercept=0, color="#3C3C3C")
        
        ##add basic plot structure
        g_mobj <- g_mobj + 
            geom_rect() +
            scale_y_continuous(breaks=NULL, limits=c(-1, 7)) +
            ylab("") +
            xlab(ifelse(genomic, paste0("Genomic Coordinates, ", seqnames(gr_e[1])),
                        "non-genomic coordinates")) +
            theme_bw()
    
        ##frame exons if not using black background
        g_mobj <- g_mobj +
            annotate("rect", size=.25,
                     xmin=start(ranges(gr_e)) - 1,
                     xmax=end(ranges(gr_e)) + 1,
                     ymin=-1, ymax=1,
                     alpha=1, color="#3C3C3C", fill=NA)

        if (iflip) { g_mobj <- g_mobj + scale_x_reverse() }

        ##add splicing arrows to plot
        width_props <- width(ranges(gr_j)) / width(range(gr_e))
        for (j in 1:p_j) {
            w_prop <- width_props[j]
            circle1 <- .pseudoArc(xmin=start(ranges(gr_j))[j],
                                  xmax=end(ranges(gr_j))[j],
                                  ymin=1, height=6*sqrt(w_prop))
            circle2 <- cbind(circle1, xmin=gg_e$xmin[1], xmax=gg_e$xmax[1],
                             ymin=gg_e$xmin[1], ymax=gg_e$ymax[1],
                             value=0)
            g_mobj <- g_mobj +
                geom_path(data=circle2, size=.25, aes(x=x, y=y),
                          arrow=grid::arrow(length=grid::unit(.025, "npc"),
                              ends=arrowhead[j]))
        }
    }

    
    ##add annotations if txdb was passed to function
    if (!is.null(txlist) && genomic) {
        cand_names <- concomp2name(obj, txlist, txdb, orgdb)
        cand_idx <- cand_names[, 1]
        cand_names <- cand_names[, 2]
        
        if (length(cand_idx) > 0) {
            tx_match <- txlist[cand_idx]
            names(tx_match) <- make.unique(cand_names)
            annot_track <- ggplot(tx_match) +
                geom_alignment() + theme_bw()
            if (iflip) { annot_track <- annot_track + scale_x_reverse() }

            itxdb <- TRUE
        }
    } else {
        itxdb <- FALSE
    }

    
    ##combine generated tracks
    if (imodel && itxdb) {
        g_obj <- tracks(g_obj, g_mobj, annot_track, heights=c(2, 1, 1))
    } else if (imodel) {
        g_obj <- tracks(g_obj, g_mobj, heights=c(2, 1))
    } else if (itxdb) {
        g_obj <- tracks(g_obj, annot_track, heights=c(2, 1))
    }

    
    g_obj
}



    

#' @rdname splicepcp
#' @aliases splicepcp,concomp-method
setMethod("splicepcp",
          signature(obj = "concomp"),
          function(obj, ...) .splicepcp.concomp(obj, ...))
