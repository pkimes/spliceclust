#' construct ggplot2 object for plotting
#'
#' @param gr_e \code{GenomicRanges} for exons
#' @param e_loads see \code{splicegralp} documentation
#'
#' @keywords internal
#' @author Patrick Kimes
sl_create <- function(gr_e, e_loads) {

    ##construct ggplot2 object
    sl_df <- data.frame(xmin=start(ranges(gr_e)),
                       xmax=end(ranges(gr_e)),
                       e_loads)
    sl_df <- reshape2::melt(sl_df, id.vars=c("xmin", "xmax"))
    sl_df$ymin <- -1
    sl_df$ymax <- 1

    sl_df
}
    


#' construct base ggplot2 object for splicegralp plot
#'
#' @param sl_df data.frame output from \code{sl_create}
#' @param gr_e \code{GenomicRanges} for exons
#' @param gr_j \code{GenomicRanges} for junctions
#' @param p_e number of exons in \code{gr_e}
#' @param p_e number of junctions in \code{gr_j}
#' @param genomic see \code{splicegralp} documentation
#' @param use_blk see \code{splicegralp} documentation
#' @param bin see \code{splicegralp} documentation
#' 
#' @keywords internal
#' @author Patrick Kimes
sl_drawbase <- function(sl_df, gr_e, gr_j, p_e, p_j, e_loads, j_loads, 
                        genomic, use_blk, load_lims, n_vec) {

    ##color generators
    pal <- plot_colors()
    crp <- pal$col4
    crp2 <- pal$col5
    
    ##strand of junctions for arrow heads
    arrowhead <- ifelse(as.character(strand(gr_j)) == "-", "last", "first")

    ##plot on genomic coordinates
    g_obj <- ggplot(sl_df, aes(xmin=xmin, xmax=xmax,
                               ymin=ymin, ymax=ymax,
                               color=value, fill=value))


    ## ##add horizontal line first
    ## g_obj <- g_obj + 
    ##     geom_hline(yintercept=0,
    ##                color=ifelse(use_blk, "#F0F0F0", "#3C3C3C"))

    
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
        scale_color_gradient2("Loading", limits=color_lim,
                              low="#053061", mid="#f7f7f7", high="#67001f",
                              guide="none") +
        scale_fill_gradient2("Loading", limits=color_lim,
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
            circle2 <- cbind(circle1, xmin=sl_df$xmin[1], xmax=sl_df$xmax[1],
                             ymin=sl_df$xmin[1], ymax=sl_df$ymax[1],
                             variable=paste0("Dir", ipc), value=0)

            ##only include arrows if direction is known
            if (all(strand(gr_e) == "*")) {
                g_obj <- g_obj +
                    geom_path(data=circle2, size=.75, aes(x=x, y=y),
                              color=.rgb2hex(crp2(j_col_scale[j, ipc])))
            } else {
                g_obj <- g_obj +
                    geom_path(data=circle2, size=.75, aes(x=x, y=y),
                              color=.rgb2hex(crp2(j_col_scale[j, ipc])),
                              arrow=grid::arrow(length=grid::unit(.015*n_vec, "npc"),
                                  ends=arrowhead[j]))
            }
        }
    }

    g_obj
}


