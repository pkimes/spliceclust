#' construct ggplot2 object for plotting
#'
#' @param gr_e \code{GenomicRanges} for exons
#' @param gr_j \code{GenomicRanges} for junctions
#' @param vals_e matrix of exon coverages
#' @param vals_j matrix of junction coverages
#' @param j_incl see \code{splicegrahm} documentation
#' @param log_base see \code{splicegrahm} documentation
#' @param log_shift see \code{splicegrahm} documentation
#' @param bin see \code{splicegrahm} documentation
#' @param n number of samples in \code{vals_e}, \code{vals_j}
#' @param p_j number of junctions in \code{gr_j}
#' @param y_flip logical whether model should be flipped on
#'        vertical axis (default = FALSE)
#' @param same_scale_n number of samples that should be used to set vertical scaling of
#'        \code{splicegrahm2} plot (if \code{splicegrahm2} parameter
#'        \code{same_scale = FALSE}, then simply \code{n}) (default = n)
#'
#' @keywords internal
#' @author Patrick Kimes
sg_create <- function(gr_e, gr_j, vals_e, vals_j, j_incl,
                      log_base, log_shift, bin, n, p_j,
                      y_flip = FALSE, same_scale_n = n) {
    
    sg_df <- data.frame(xmin=start(ranges(gr_e)),
                        xmax=end(ranges(gr_e)),
                        vals_e)
    sg_df <- reshape2::melt(sg_df, id.vars=c("xmin", "xmax"))
    if (length(gr_e) == 1) { sg_df$variable <- 1:n }
    sg_df$ymin <- as.numeric(sg_df$variable) + .25
    sg_df$ymax <- sg_df$ymin + .50

    
    ##transform if desired
    if (log_base > 0) {
        sg_df$value <- log(log_shift + sg_df$value, base=log_base)
    }
    
    ##perform binning
    if (bin && log_base > 0) {
        sg_df$value <- factor(floor(sg_df$value))
    }
    sg_df$kind <- "e"

    ##include junctions if necessary
    if (j_incl && p_j > 0) {
        junc_x <- seq(min(min(ranges(gr_e))),
                      max(max(ranges(gr_e))),
                      length.out=p_j+2)
        junc_x <- junc_x[2:(p_j+1)]
        
        junc_y <- 2.25 * same_scale_n
        s_size <- .5
        junc_w <- width(range(gr_e)) / (2.5*(max(p_j, 5)+1))

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
            gg_j$value <- factor(floor(gg_j$value), levels=0:max(floor(gg_j$value)))
        }

        sg_df <- rbind(sg_df, gg_j)
    }

    if (y_flip) {
        temp <- sg_df$ymin
        sg_df$ymin <- -sg_df$ymax
        sg_df$ymax <- -temp
    }        

    return(sg_df)
}



#' construct base ggplot2 object for splicegrahm plot
#'
#' @param sg_df data.frame output from \code{sg_create}
#' @param gr_e \code{GenomicRanges} for exons
#' @param use_blk see \code{splicegrahm} documentation
#' @param j_incl see \code{splicegrahm} documentation
#' @param genomic see \code{splicegrahm} documentation
#' @param log_base see \code{splicegrahm} documentation
#' @param bin see \code{splicegrahm} documentation
#' @param highlight see \code{splicegrahm} documentation
#' @param n number of samples in \code{vals_e}, \code{vals_j}
#' @param p_j number of junctions in \code{gr_j}
#' @param iflip logical whether model will be on negative strand
#' @param mirror logical whether model should be flipped on
#'        vertical axis (defualt = FALSE)
#' @param same_scale_n number of samples that should be used to set vertical scaling of
#'        \code{splicegrahm2} plot (if \code{splicegrahm2} parameter
#'        \code{same_scale = FALSE}, then simply \code{n}) (default = n)
#' 
#' @keywords internal
#' @author Patrick Kimes
sg_drawbase <- function(sg_df, use_blk, j_incl, genomic, gr_e,
                        log_base, bin, n, highlight, p_j, iflip,
                        mirror = FALSE, same_scale_n = n) {

    pal <- plot_colors()
    hl_cols <- pal$col3
    
    ##add fake scale of ones for alpha plotting
    sg_df$ones <- factor(rep(1, length.out=nrow(sg_df)),
                         levels=seq(0, 1, .2))#rep(seq(0, 1, .2), length.out=nrow(sg_df))
    
    ##base of plot
    sg_obj <- ggplot(sg_df, aes(xmin=xmin, xmax=xmax,
                                ymin=ymin, ymax=ymax,
                                color=value, fill=value,
                                alpha=ones))

    ##add highlighting of cluster assignments if necessary
    if (!is.null(highlight)) {
        hl_tab <- c(1, table(highlight))
        hl_h <- cumsum(hl_tab)
        k <- length(hl_tab) - 1
        
        wrc <- width(range(gr_e))
        minrc <- min(c(start(gr_e), end(gr_e)))
        maxrc <- max(c(start(gr_e), end(gr_e)))
        
        for (i in 1:k) {
            sg_obj <- sg_obj +
                annotate("rect", size=.125, 
                         xmin=c(minrc-wrc*0.05, maxrc),
                         xmax=c(minrc, maxrc+wrc*0.05),
                         ymin=rep(hl_h[i], 2)+.25, ymax=rep(hl_h[i+1], 2)-.25,
                         fill=hl_cols[i], color=hl_cols[i], alpha=1)
        }
    }
    
    
    ##add basic plot structure
    sg_obj <- sg_obj + 
        geom_rect(size=.125 + 0.055) +
        scale_alpha_discrete("splicing", range=0:1,
                             labels=paste0(seq(0, 100, 20), "%"), drop=FALSE) + 
        xlab(ifelse(genomic, paste0("Genomic Coordinates, ", seqnames(gr_e[1])),
                    "non-genomic coordinates")) +
        theme_bw()

    
    ##adjust y-axis scaling
    sg_obj <- sg_obj +
        scale_y_continuous("", breaks=NULL,
                           limits=ifelse(c(mirror, mirror),
                               c(-(2.15+j_incl)*same_scale_n, 1),
                               c(-1, (2.15+j_incl)*same_scale_n)))
    
    
    ##frame exon heatmaps with box if not using black background
    if (use_blk) {
        sg_obj <- sg_obj +
            theme(panel.grid.minor.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.background = element_rect(fill="#3C3C3C"))
    } else {
        sg_obj <- sg_obj +
            annotate("rect", size = .125,
                     xmin = start(ranges(gr_e)) - .25,
                     xmax = end(ranges(gr_e)) + .25,
                     ymin = ifelse(mirror, -(n+1+.25), .75-.125),
                     ymax = ifelse(mirror, -(.75-.125), n+1+.25),
                     alpha = 1, color = "#3C3C3C", fill = NA)
    }

    ##add continuous or discrete color palette
    if (log_base > 0 && bin) {
        v_max <- length(levels(sg_df$value)) - 1
        sg_obj <- sg_obj + 
            scale_color_manual("expr", breaks=levels(sg_df$value), values=pal$col1(v_max+1), guide="none") +
            scale_fill_manual("expr", breaks=levels(sg_df$value), values=pal$col1(v_max+1),
                              labels=paste0("<", log_base^(1:(v_max+1))), drop=FALSE)
    } else {
        sg_obj <- sg_obj +
            scale_color_continuous("expr", low="#f7fbff", high="#08306b", guide="none") +
            scale_fill_continuous("expr", low="#f7fbff", high="#08306b")
    }


    ##adjust the size and shape of text labels
    sg_obj + guides(fill=guide_legend(
                        keywidth = 2,
                        keyheight = .5,
                        reverse=TRUE,
                        label.theme=element_text(size=rel(8), angle=0)),
                    alpha=guide_legend(
                        keywidth=2,
                        keyheight=.5,
                        reverse=TRUE,
                        override.aes=list(color="#08306b", fill=.rgb2hex(pal$col2(1))),
                        label.theme=element_text(size=rel(8), angle=0)))
    
}



#' add junction annotations to base ggplot2 object for splicegrahm plot
#'
#' @param sg_obj \code{ggplot} object created by \code{sg_drawbase}
#' @param sg_df data.frame output from \code{sg_create}
#' @param iflip logical values whether the plotting should be reversed
#' @param gr_e \code{GenomicRanges} for exons
#' @param gr_j \code{GenomicRanges} for junctions
#' @param vals_j matrix of junction coverages
#' @param n number of samples in \code{vals_e}, \code{vals_j}
#' @param p_j number of junctions in \code{gr_j}
#' @param j_incl see \code{splicegrahm} documentation
#' @param use_blk see \code{splicegrahm} documentation
#' @param highlight see \code{splicegrahm} documentation
#' @param mirror logical whether model should be flipped on
#'        vertical axis (defualt = FALSE)
#' @param same_scale_n number of samples that should be used to set vertical scaling of
#'        \code{splicegrahm2} plot (if \code{splicegrahm2} parameter
#'        \code{same_scale = FALSE}, then simply \code{n}) (default = n)
#' 
#' @keywords internal
#' @author Patrick Kimes
sg_drawjuncs <- function(sg_obj, sg_df, j_incl, use_blk, iflip,
                         gr_e, gr_j, vals_j, n, p_j, highlight,
                         mirror = FALSE, same_scale_n = n) {

    ##strand of junctions for arrow heads
    arrowhead <- ifelse(as.character(strand(gr_j)) == "-", "last", "first")
    
    ##color generators
    pal <- plot_colors()
    hl_cols <- pal$col3

    ##long list of letters
    ll_LETTERS <- long_letters()
    
    ##add splicing arrows to plot
    e_prop <- rowMeans(vals_j > 0) 
    w_prop <- width(ranges(gr_j)) / width(range(gr_e))

    ##whether to flip arrow plotting
    iud <- ifelse(mirror, -1, 1)
    
    for (j in 1:p_j) {
        circle1 <- .pseudoArc(xmin=start(ranges(gr_j))[j],
                              xmax=end(ranges(gr_j))[j],
                              ymin=n+1, height=n*sqrt(w_prop[j]))

        ##only include arrows if direction is known
        if (all(strand(gr_e) == "*")) {
            sg_obj <- sg_obj +
                annotate("path", size=.75,
                         x=circle1$x, y=iud*circle1$y,
                         color=.rgb2hex(pal$col2(e_prop[j])))
        } else {
            sg_obj <- sg_obj +
                annotate("path", size=.75,
                         x=circle1$x, y=iud*circle1$y,
                         color=.rgb2hex(pal$col2(1)), alpha = e_prop[j],
                         arrow=grid::arrow(length=grid::unit(.015, "npc"), ends=arrowhead[j]))
        }
    }


    ##add text labels for splicing arrows
    if (j_incl) {
        junc_x <- seq(min(min(ranges(gr_e))),
                      max(max(ranges(gr_e))),
                      length.out=p_j+2)
        junc_x <- junc_x[2:(p_j+1)]
        junc_y <- 2.25*same_scale_n
        s_size <- .5
        junc_w <- width(range(gr_e)) / (2.5*(max(p_j, 5)+1))
    
        abc <- ll_LETTERS[1:p_j]
        if (iflip) { abc <- rev(abc) }

        sg_obj <- sg_obj +
            annotate("text", size=3,
                     x=(start(ranges(gr_j)) + end(ranges(gr_j))) / 2,
                     y=iud*(n+3)*(1 + sqrt(w_prop)), vjust=0+mirror, 
                     label=abc,
                     color=ifelse(use_blk, "#F0F0F0", "#3C3C3C"))

        sg_obj <- sg_obj +
            annotate("text", size=3,
                     x=junc_x,
                     y=iud*rep(junc_y + (n+4)*(.025 + s_size), p_j),
                     vjust=0+mirror,
                     label=abc,
                     color=ifelse(use_blk, "#F0F0F0", "#3C3C3C"))


        ##add highlighting if necessary
        if (!is.null(highlight)) {
            hl_tab <- c(1, table(highlight))
            hl_h <- cumsum(hl_tab)
            k <- length(hl_tab) - 1
            
            for (i in 1:k) {
                sg_obj <- sg_obj +
                    annotate("rect", size = .125,
                             xmin=c(min(junc_x)-junc_w*2, max(junc_x)+junc_w+1),
                             xmax=c(min(junc_x)-junc_w-1, max(junc_x)+junc_w*2),
                             ymin=(rep(hl_h[i], 2)+.25)*s_size + junc_y,
                             ymax=(rep(hl_h[i+1], 2)-.25)*s_size + junc_y,
                             fill=hl_cols[i], color=hl_cols[i], alpha=1)
            }
        }

        ##add rectangles around junction heatmaps
        if (!use_blk) {
            sg_obj <- sg_obj +
                annotate("rect", size = .125,
                         xmin = junc_x - junc_w - 1, xmax = junc_x + junc_w + 1,
                         ymin=iud*((.75 - .125)*s_size + junc_y),
                         ymax=iud*((n+1+.25)*s_size + junc_y),
                         alpha=1, color=.rgb2hex(pal$col2(1)), fill=NA)
        }
    }

    sg_obj
}
