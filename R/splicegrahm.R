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
                                 use_blk = FALSE, txlist = NULL,
                                 txdb = NULL, orgdb = NULL, ...) {

    ##exonValues and juncValues must be specified
    if (is.null(exonValues(obj)) || is.null(juncValues(obj)))
        stop(paste0("exonValues and juncValues cannot be NULL for splicegrahm, \n",
                    "consider using splicegralp instead."))


    ##can't include gene models if not plotting on genomic scale
    if (!is.null(txlist) && !genomic) {
        cat("since txlist provided, plotting on genomic scale. \n")
        genomic <- TRUE
    }
    
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


    ##determine overlapping gene models
    if (!is.null(txlist) && genomic) {
        cand_names <- concomp2name(obj, txlist, txdb, orgdb)
        cand_idx <- cand_names[, 1]
        cand_names <- cand_names[, 2]
        
        if (length(cand_idx) > 0) {
            tx_match <- txlist[cand_idx]
            names(tx_match) <- make.unique(cand_names)
            annot_track <- ggplot(tx_match) + geom_alignment() + theme_bw()
        }
    }

    
    ##determine whether plots should be flipped based on passed data
    ## or matched models -- remove arrows if not passed !!!!!
    if (all(strand(gr_e) == "*") &&
        !is.null(txlist) &&
        length(cand_idx) > 0) {
        iflip <- flip_neg && all(strand(tx_match) == '-')
    } else {
        iflip <- flip_neg && all(strand(gr_e) == '-')
    }
    
    
    ##determine order of samples
    if (sort_sep) {
        vals_e <- t(apply(vals_e, 1, sort))
        vals_j <- t(apply(vals_j, 1, sort))

    } else {
        if (length(sort_idx) == 1) {
            if (sort_idx == 1) {
                idx <- order(vals_e[1, ])

            } else if (sort_idx == 2) {
                pca <- prcomp(t(vals_e))
                idx <- order(pca$x[, 2])

            } else {
                idx <- 1:n
            }
            
        } else if (length(sort_idx) == n) {
            idx <- sort_idx

        } else {
            idx <- 1:n
        }

        vals_e <- vals_e[, idx]
        vals_j <- vals_j[, idx]
    }

    
    ##strand of junctions for arrow heads
    arrowhead <- ifelse(as.character(strand(gr_j)) == "-",
                        "last", "first")
    

    ##construct ggplot2 object
    gg_e <- data.frame(xmin=start(ranges(gr_e)),
                       xmax=end(ranges(gr_e)),
                       vals_e)
    gg_e <- reshape2::melt(gg_e, id.vars=c("xmin", "xmax"))
    gg_e$ymin <- as.numeric(gg_e$variable)
    gg_e$ymax <- gg_e$ymin + 1

    
    ##transform if desired
    if (log_base > 0) {
        gg_e$value <- log(log_shift + gg_e$value, base=log_base)
        v_max <- max(gg_e$value)
    }
    
    
    ##perform binning
    if (bin && log_base > 0) {
        gg_e$value <- factor(floor(gg_e$value))
    }
            
    
    ##include junctions if necessary
    if (j_incl) {
        junc_x <- seq(min(min(ranges(gr_e))),
                      max(max(ranges(gr_e))),
                      length.out=p_j+2)
        junc_x <- junc_x[2:(p_j+1)]
        
        ##flip order of plotting to match negative strand
        ## if (iflip) {
        ##     junc_x <- rev(junc_x)
        ## }
        
        junc_y <- 2.25*n
        s_size <- .5
        junc_w <- width(range(gr_e)) / (2.5*(max(p_j, 5)+1))
        ## junc_w <- min(junc_w, 50)

        ##create data.frame
        gg_j <- data.frame(xmin=junc_x - junc_w, xmax=junc_x + junc_w, vals_j)
        gg_j <- reshape2::melt(gg_j, id.vars=c("xmin", "xmax"))
        gg_j$ymin <- as.numeric(gg_j$variable)*s_size + junc_y
        gg_j$ymax <- gg_j$ymin + s_size
        gg_j$kind <- "j"

        ##transform if desired
        if (log_base > 0) {
            gg_j$value <- log(log_shift + gg_j$value, base=log_base)
        }    

        ##perform binning
        if (bin && log_base > 0) {
            gg_j$value <- factor(floor(gg_j$value))
        }

    }


    ##merge exon and junction data
    if (j_incl) {
        gg_e$kind <- "e"
        gg_e <- rbind(gg_e, gg_j)
    }

    
    ##plot on genomic coordinates
    g_obj <- ggplot(gg_e, aes(xmin=xmin, xmax=xmax,
                              ymin=ymin, ymax=ymax,
                              color=value, fill=value))


    ##add horizontal line first
    g_obj <- g_obj + 
        geom_hline(yintercept=n/2,
                   color=ifelse(use_blk, "#F0F0F0", "#3C3C3C"))


    ##add annotation colors if specified
    if (!is.null(highlight)) {
        wrc <- width(range(gr_e))
        minrc <- min(c(start(gr_e), end(gr_e)))
        maxrc <- max(c(start(gr_e), end(gr_e)))
        for (i in 1:n) {
            if (highlight[i] > 0)
                g_obj <- g_obj +
                    annotate("rect",
                             xmin=c(minrc-wrc*0.05, maxrc),
                             xmax=c(minrc, maxrc+wrc*0.05),
                             ymin=c(i, i), ymax=c(i+1, i+1),
                             fill=hl_cols[highlight[i]], alpha=1)
        }
        
        ##also add for junciton level
        if (j_incl) {
            for (i in 1:n) {
                if (highlight[i] > 0)
                    g_obj <- g_obj +
                        annotate("rect",
                                 xmin=c(min(junc_x)-junc_w*2, max(junc_x)+junc_w),
                                 xmax=c(min(junc_x)-junc_w, max(junc_x)+junc_w*2),
                                 ymin=c(i, i)*s_size+junc_y, ymax=c(i+1, i+1)*s_size+junc_y,
                                 fill=hl_cols[highlight[i]], alpha=1)
            }
        }
    }

    
    ##add basic plot structure
    g_obj <- g_obj + 
        geom_rect() +
        scale_y_continuous(breaks=NULL, limits=c(0, (2.15+j_incl)*n)) +
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
                     ymin=1, ymax=n+1,
                     alpha=1, color="#3C3C3C", fill=NA)
    }


    
    ##add continuous or discrete color palette
    if (log_base > 0 && bin) {
        g_obj <- g_obj + 
            scale_color_manual("expr", breaks=0:v_max, values=crp(v_max+1), guide="none") +
            scale_fill_manual("expr", breaks=0:v_max, values=crp(v_max+1),
                              labels=paste0("<", log_base^(1:(v_max+1))))
    } else {
        g_obj <- g_obj +
            scale_color_continuous("expr", low="#f7fbff", high="#08306b", guide="none") +
            scale_fill_continuous("expr", low="#f7fbff", high="#08306b")
    }

    
    ##add splicing arrows to plot
    e_prop <- rowMeans(vals_j > 0) 
    for (j in 1:p_j) {
        w_prop <- width(ranges(gr_j)[j]) / width(range(gr_e))
        circle1 <- .pseudoArc(xmin=start(ranges(gr_j))[j],
                              xmax=end(ranges(gr_j))[j],
                              ymin=n+1, height=1*n*sqrt(w_prop))

        ##only include arrows if direction is known
        if (all(strand(gr_e) == "*")) {
            g_obj <- g_obj +
                annotate("path", size=.75,
                         x=circle1$x, y=circle1$y,
                         color=.rgb2hex(crp2(e_prop[j])))
        } else {
            g_obj <- g_obj +
                annotate("path", size=.75,
                         x=circle1$x, y=circle1$y,
                         color=.rgb2hex(crp2(e_prop[j])),
                         arrow=grid::arrow(length=grid::unit(.015, "npc"), ends=arrowhead[j]))
        }
    }
    

    ##add text labels for splicing arrows
    if (j_incl) {
        abc <- ll_LETTERS[1:p_j]
        if (iflip) { abc <- rev(abc) }

        w_prop <- width(ranges(gr_j)) / width(range(gr_e))
        g_obj <- g_obj +
            annotate("text", size=3,
                     x=(start(ranges(gr_j)) + end(ranges(gr_j))) / 2,
                     y=(n+1)*(1 + sqrt(w_prop)), vjust=0, 
                     label=abc,
                     color=ifelse(use_blk, "#F0F0F0", "#3C3C3C"))

        g_obj <- g_obj +
            annotate("text", size=3,
                     x=junc_x, y=rep(junc_y + (n+1)*(.025 + s_size), p_j),
                     label=abc, vjust=0,
                     color=ifelse(use_blk, "#F0F0F0", "#3C3C3C"))

        ##add rectangles around 
        if (!use_blk) {
            g_obj <- g_obj +
                annotate("rect", size=.25,
                         xmin=junc_x - junc_w - 1, xmax=junc_x + junc_w + 1,
                         ymin=junc_y + s_size, ymax=junc_y + (n+1)*s_size,
                         alpha=1, color=.rgb2hex(crp2(1)), fill=NA)
        }
        
    }


    ##plot with horizontal axis oriented on negative strand
    if (iflip) {
        g_obj <- g_obj + scale_x_reverse()
    }


    ##add annotations if txdb was passed to function
    if (!is.null(txlist) && genomic && length(cand_idx) > 0) {
        if (iflip) { annot_track <- annot_track + scale_x_reverse() }
        g_obj <- tracks(g_obj, annot_track, heights=c(2, 1))
    }

    g_obj
}

#' @keywords internal
#' @title splicegrahm method
#' @name splicegrahm-concomp
#' @aliases splicegrahm,concomp-method
setMethod("splicegrahm",
          signature(obj = "concomp"),
          function(obj, ... ) .splicegrahm.concomp(obj, ...))



