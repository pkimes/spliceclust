#' Splice Graphs Heatmaps
#'
#' Splicing graphs as disjoint heatmaps for the comparative visualization
#' of multiple splice graphs in a single figure
#'
#' @param cc a connected component GRanges object with exon and junction
#'        information
#' @param sort_sep a logical whether to sort each exon, junction separately
#'        (default = FALSE)
#' @param sort_idx an integer value specifying the order of the samples in
#'        each exon, see details for more information on all possible
#'        input (default = 1)
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
#' @export
#' @import ggplot2 RColorBrewer reshape2 grid
#' @author Patrick Kimes


library("ggbio")
library("GenomicRanges")
library("ggplot2")
library("RColorBrewer")
library("reshape2")
library("grid")


heatgraph <- function(cc, sort_sep = FALSE, sort_idx = 1,
                      log_base = 10, log_shift = 1, bin = TRUE,
                      genomic = TRUE, ex_use = 2/3, flip_neg = TRUE, 
                      j_incl = FALSE, highlight = NULL,
                      use_blk = FALSE, ...) {

    
    ##color generators
    crp <- colorRampPalette(c("#f7fbff", "#08306b"))
    crp2 <- colorRamp(c("#f7fbff", "#047760"))

    
    ##parse input connected component
    s_col <- grepl("^s", names(values(cc)))
    depth <- as.data.frame(values(cc)[s_col])


    ##split cc into exon and junction objects
    exon_row <- values(cc)$kind == "e"
    junc_row <- !exon_row

    
    ##dataset dimension
    n <- ncol(depth)
    p_e <- sum(exon_row)
    p_j <- sum(junc_row)

    
    ##change GRanges coordinates if non-genomic coordinates are desired
    if (!genomic) {
        dna_len <- width(range(cc))
        rna_len <- sum(width(ranges(cc))[exon_row])
        cc2 <- cc
        ranges(cc2) <- shift(ranges(cc2),
                             -start(ranges(cc2)[exon_row][1])+1)

        cc2_gaps <- distance(ranges(cc2)[exon_row][-p_e],
                             ranges(cc2)[exon_row][-1]) *
                                 (1 - (1-ex_use)/ex_use *
                                      rna_len/(dna_len-rna_len))

        eeb <- start(ranges(cc2)[exon_row])[-1]
        for (ii in (p_e-1):1) {
            start(ranges(cc2))[start(ranges(cc2)) >= eeb[ii]] <-
                start(ranges(cc2))[start(ranges(cc2)) >= eeb[ii]] - cc2_gaps[ii]
            end(ranges(cc2))[end(ranges(cc2)) >= eeb[ii]] <-
                end(ranges(cc2))[end(ranges(cc2)) >= eeb[ii]] - cc2_gaps[ii]
        }
        cc <- cc2
    }


    ##if negative strand ordering, flip order of cc
    if (flip_neg && all(as.vector(strand(cc))[exon_row] == '-')) {
        cc <- rev(cc)
        exon_row <- rev(exon_row)
        junc_row <- rev(junc_row)
        depth <- depth[(p_e+p_j):1, ]
    }

    
    ##determine order of samples
    if (sort_sep) {
        depth <- t(apply(depth, 1, sort))
    } else {
        if (sort_idx == 1) {
            idx <- order(depth[1, ])
        } else if (sort_idx == 2) {
            pca <- prcomp(t(depth[exon_row, ]))
            idx <- order(pca$x[, 2])
        } else {
            idx <- 1:n
        }
        depth <- depth[, idx]
    }

    
    ##strand of junctions for arrow heads
    arrowhead <- ifelse(as.character(strand(cc))[junc_row] == "+",
                        "first", "last")
    

    ##construct ggplot2 object
    gg_e <- data.frame(xmin=start(ranges(cc)[exon_row]),
                       xmax=end(ranges(cc)[exon_row]),
                       depth[exon_row, ])
    gg_e <- melt(gg_e, id.vars=c("xmin", "xmax"))
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
        
    
    ##plot on genomic coordinates
    g_obj <- ggplot(gg_e, aes(xmin=xmin, xmax=xmax,
                              ymin=ymin, ymax=ymax,
                              color=value, fill=value))


    ##add annotation colors if specified
    if (!is.null(highlight)) {
        hl_cols <- brewer.pal(8, "Set2")
        for (i in 1:n) {
            if (highlight[i] > 0)
                g_obj <- g_obj +
                    annotate("rect", xmin=-width(range(cc))*0.1, xmax=width(range(cc))*1.1,
                             ymin=i, ymax=i+1, fill=hl_cols[highlight[i]],
                             alpha=1/2)
        }
    }

    
    ##add basic plot structure
    g_obj <- g_obj + 
        geom_hline(yintercept=n/2, color=ifelse(use_blk, "#F0F0F0", "#3C3C3C")) +
        geom_rect() +
        scale_y_continuous(breaks=NULL, limits=c(0, (2.1+j_incl)*n)) +
        ylab("") +
        xlab(ifelse(genomic, paste0("Genomic Coordinates, ", seqnames(cc[1])),
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
                     xmin=start(ranges(cc))[exon_row]-2,
                     xmax=end(ranges(cc))[exon_row]+2,
                     ymin=0, ymax=n+1,
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
    w_t <- diff(range(as.vector(ranges(cc))))
    e_prop <- rowMeans(depth[junc_row, ] > 0) 
    for (j in 1:p_j) {
        w_prop <- width(ranges(cc)[junc_row][j]) / w_t
        circle1 <- .pseudoArc(xmin=start(ranges(cc))[junc_row][j],
                              xmax=end(ranges(cc))[junc_row][j],
                              ymin=n+2, height=1*n*sqrt(w_prop))
        g_obj <- g_obj +
            annotate("path", size=.75,
                     x=circle1$x, y=circle1$y,
                     color=.rgb2hex(crp2(e_prop[j])),
                     arrow=arrow(length=unit(.015, "npc"), ends=arrowhead[j]))
        if (j_incl) {
            g_obj <- g_obj +
                annotate("text", size=3,
                         x=(start(ranges(cc))[junc_row][j] + end(ranges(cc))[junc_row][j])/2,
                         y=n + 5 + 1*n*sqrt(w_prop), vjust=0, 
                         label=LETTERS[j], color=ifelse(use_blk, "#F0F0F0", "#3C3C3C"))
        }
    }


    ##include junctions if necessary
    if (j_incl) {

        junc_x <- seq(min(min(ranges(cc))),
                      max(max(ranges(cc))),
                      length.out=p_j+2)[2:(p_j+1)]

        ##flip order of plotting to match negative strand
        if (flip_neg && all(as.vector(strand(cc)[exon_row]) == "-")) {
            junc_x <- rev(junc_x)
        }
        
        junc_y <- 2.25*n
        s_size <- .5

        ##create data.frame
        gg_j <- data.frame(xmin=junc_x - 50, xmax=junc_x + 50, depth[junc_row, ])
        gg_j <- melt(gg_j, id.vars=c("xmin", "xmax"))
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
        gg_ej <- gg_e
        gg_ej$kind <- "e"
        gg_ej <- rbind(gg_ej, gg_j)
        
        g_obj <-
            g_obj %+% gg_ej +
            annotate("text", size=3,
                     x=junc_x, y=rep(junc_y + n*s_size + 5, p_j),
                     label=LETTERS[1:p_j], vjust=0,
                     color=ifelse(use_blk, "#F0F0F0", "#3C3C3C"))
        if (!use_blk) {
            g_obj <- g_obj +
                annotate("rect", size=.25,
                         xmin=junc_x - 51, xmax=junc_x + 51,
                         ymin=junc_y - 2, ymax=junc_y + n*s_size + 2,
                         alpha=1, color=.rgb2hex(crp2(1)), fill=NA)
        }
                             
    }


    ##plot with horizontal axis oriented on negative strand
    if (flip_neg && all(as.vector(strand(cc))[exon_row] == '-')) {
        g_obj <- g_obj + scale_x_reverse()
    }

    g_obj
}





## ################################################################################

##helper function for adding arcs to ggplot2 object via geom_path()
## modified from: http://stackoverflow.com/questions/6862742/
.pseudoArc <- function(xmin = 0, xmax = 10,
                       ymin = 0, height = 10, npoints = 100) {
    x_c  <- (xmax + xmin) / 2
    x_r <- (xmax - xmin) / 2
    tt <- seq(0, pi, length.out = npoints)
    
    xx <- x_c + x_r * cos(tt)
    yy <- ymin + height * sin(tt)
    return(data.frame(x = xx, y = yy))
}


.rgb2hex <- function(rgb) {
    if (is.matrix(rgb)) {
        apply(round(rgb), 1,
              function(x) { 
                  paste0("#", paste0(as.hexmode(x), collapse=""))
              })
    } else {
        paste0("#", paste0(as.hexmode(round(rgb)), collapse=""))
    }        
}
