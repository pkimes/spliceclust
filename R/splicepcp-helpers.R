#' construct ggplot2 object for plotting
#'
#' @param gr_e \code{GenomicRanges} for exons
#' @param gr_j \code{GenomicRanges} for junctions
#' @param vals_e matrix of exon coverages
#' @param vals_j matrix of junction coverages
#' @param n number of samples in \code{vals_e}, \code{vals_j}
#' @param p_j number of junctions in \code{gr_j}
#' @param log_base see \code{splicegrahm} documentation
#' @param log_shift see \code{splicegrahm} documentation
#'
#' @keywords internal
#' @author Patrick Kimes
sp_create <- function(gr_e, gr_j, vals_e, vals_j,
                      log_base, log_shift, n, p_e, p_j) {

    ##exon expression data
    sp_df <- data.frame(xmin=start(ranges(gr_e)),
                       xmax=end(ranges(gr_e)),
                       vals_e)
    sp_df <- reshape2::melt(sp_df, id.vars=c("xmin", "xmax"))

    sp_df$ymin <- sp_df$value
    sp_df$ymax <- sp_df$value
    sp_df$transp <- 2/3
    sp_df <- sp_df[c("xmin", "xmax", "variable", "ymin", "ymax", "transp")]

    
    ##between-exon start data
    sp_gap_df1 <- data.frame(xmin = end(ranges(gr_e))[-p_e],
                             xmax = start(ranges(gr_e))[-1],
                             vals_e[-p_e, ])
    sp_gap_df1 <- reshape2::melt(sp_gap_df1, id.vars=c("xmin", "xmax"))
    names(sp_gap_df1) <- c("xmin", "xmax", "variable", "ymin")
    
    ##between-exon end data
    sp_gap_df2 <- data.frame(xmin = end(ranges(gr_e))[-p_e],
                             xmax = start(ranges(gr_e))[-1],
                             vals_e[-1, ])
    sp_gap_df2 <- reshape2::melt(sp_gap_df2, id.vars=c("xmin", "xmax"))
    names(sp_gap_df2) <- c("xmin", "xmax", "variable", "ymax")

    ##complete between-exon data 
    sp_gap_df <- merge(sp_gap_df1, sp_gap_df2)
    sp_gap_df$transp <- 1/4
    sp_gap_df <- sp_gap_df[order(sp_gap_df$variable), ]
    

    ##combine
    sp_df <- rbind(sp_df, sp_gap_df)

    
    ##transform if desired
    if (log_base > 0) {
        sp_df$ymin <- log(log_shift + sp_df$ymin, base=log_base)
        sp_df$ymax <- log(log_shift + sp_df$ymax, base=log_base)
    }

    sp_df
}


#' construct base ggplot2 object for splicegrahm plot
#'
#' @param sp_df data.frame output from \code{sp_create}
#' @param gr_e \code{GenomicRanges} for exons
#' @param p_e number of exons in \code{gr_e}
#' @param highlight see \code{splicepcp} documentation
#' @param genomic see \code{splicepcp} documentation
#' @param log_base see \code{splicepcp} documentation
#' 
#' @keywords internal
#' @author Patrick Kimes
sp_drawbase <- function(sp_df, highlight, genomic,
                        log_base, gr_e, p_e) {

    pal <- plot_colors()

    ##color by groups if specified for highlight
    if (is.null(highlight)) {
        sp_obj <- ggplot(sp_df, aes(x=xmin, xend=xmax,
                                  y=ymin, yend=ymax))
    } else {
        sp_df$hl <- as.factor(c(rep(highlight, each=p_e),
                               rep(highlight, each=p_e-1)))
        sp_obj <- ggplot(sp_df, aes(x=xmin, xend=xmax,
                                  y=ymin, yend=ymax, color=hl)) +
                                      scale_color_brewer("Cluster", palette="Set1")
    }        

    ##add horizontal line first
    sp_obj <- sp_obj + 
        geom_hline(yintercept=0, color="#3C3C3C")

    
    ##add main segments
    sp_obj <- sp_obj + geom_segment(aes(alpha=transp), size=1/3) +
        scale_alpha_continuous(guide=FALSE)

    
    ##add basic plot structure
    sp_obj <- sp_obj + 
        ylab(ifelse(log_base == 0, "expression",
                    paste0("log", log_base, " expression"))) +
            xlab(ifelse(genomic,
                        paste0("Genomic Coordinates, ", seqnames(gr_e[1])),
                        "Non-Genomic Coordinates")) +
                            theme_bw()

    sp_obj
}



#' construct ggplot object with connected component model for splicepcp
#'
#' @param sp_df data.frame output from \code{sg_create}
#' @param gr_e \code{GenomicRanges} for exons
#' @param gr_j \code{GenomicRanges} for junctions
#' @param p_j number of junctions in \code{gr_j}
#' @param genomic see \code{splicepcp} documentation
#' 
#' @keywords internal
#' @author Patrick Kimes
sp_drawmodel <- function(sp_df, gr_e, gr_j, p_j, genomic) {

    ##strand of junctions for arrow heads
    arrowhead <- ifelse(as.character(strand(gr_j)) == "-", "last", "first")

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
                    xlab(ifelse(genomic,
                                paste0("Genomic Coordinates, ", seqnames(gr_e[1])),
                                "non-genomic coordinates")) +
                                    theme_bw()
    
    ##frame exons if not using black background
    g_mobj <- g_mobj +
        annotate("rect", size=.25,
                 xmin=start(ranges(gr_e)),
                 xmax=end(ranges(gr_e)),
                 ymin=-1, ymax=1,
                 alpha=1, color="#3C3C3C", fill=NA)

    ##add splicing arrows to plot
    width_props <- width(ranges(gr_j)) / width(range(gr_e))
    for (j in 1:p_j) {
        w_prop <- width_props[j]
        circle1 <- .pseudoArc(xmin=start(ranges(gr_j))[j],
                              xmax=end(ranges(gr_j))[j],
                              ymin=1, height=6*sqrt(w_prop))
        circle2 <- cbind(circle1, xmin=sp_df$xmin[1], xmax=sp_df$xmax[1],
                         ymin=sp_df$xmin[1], ymax=sp_df$ymax[1],
                         value=0)

        ##only include arrows if direction is known
        if (all(strand(gr_e) == "*")) {
            g_mobj <- g_mobj +
                geom_path(data=circle2, size=.25, aes(x=x, y=y))
        } else {
            g_mobj <- g_mobj +
                geom_path(data=circle2, size=.25, aes(x=x, y=y),
                          arrow=grid::arrow(length=grid::unit(.025, "npc"),
                              ends=arrowhead[j]))
        }
    }

    g_mobj
}


