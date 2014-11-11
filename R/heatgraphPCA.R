#' Splice Graphs PCA Plot
#'
#' PCA as colors on splice graph
#'
#' @param cc a connected component GRanges object with exon and junction
#'        information
#' @param npc a numeric value specifying number of PCs to plot (default = 3)
#' @param pc_sep a logical whether PCA should be performed on exon and junction
#'        coverage separately (default = TRUE)
#' @param ej_w a numeric vector of length two specifying the relative sum of squares
#'        for exon and junctions (default = c(1, 1))
#' @param log_base a numeric specifying the scale of the expression values at each exon,
#'        which 0 resulting in no log scaling being applied (default = 10)
#' @param log_shift a numeric specifying the shift to be used in the log transformation
#'        for handling 0 values (default = 1)
#' @param genomic a logical whether genomic coordinates should be used to
#'        plot the heatmap (default = TRUE)
#' @param ex_use a numeric specifying the proportion of the plot exons should occupy if
#'        non-genomic coordinate plotting is desired (default = 2/3)
#' @param flip_neg a logical whether to flip plotting of genes on negative strand
#'        to run left to right (default = TRUE)
#' @param use_blk a logical whether to use a black background (default = FALSE)
#' @param ... other parameters to be passed
#' 
#' @return
#' a ggplot2 plot showing the specified number of principal components
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


graphPCA <- function(cc, npc = 3, pc_sep = TRUE, ej_w = c(1, 1),
                     log_base = 10, log_shift = 1, genomic = TRUE, 
                     ex_use = 2/3, flip_neg = TRUE, use_blk = FALSE, ...) {

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

    
    ##color generators
    crp <- colorRampPalette(c("#053061", "#f7f7f7", "#67001f"))
    crp2 <- colorRamp(c("#053061", "#f7f7f7", "#67001f"))


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

    
    ## ##plot with horizontal axis oriented on negative strand
    ## if (flip_neg && all(as.vector(strand(cc))[exon_row] == '-')) {
    ##     t_range <- ranges(cc)
    ##     start(ranges(cc)) <- -end(t_range)
    ##     end(ranges(cc)) <- -start(t_range)
    ##     strand(cc) <- ifelse(as.vector(strand(cc)) == '-', '+', '-')
    ## }

    
    ##strand of junctions for arrow heads
    arrowhead <- ifelse(as.character(strand(cc))[junc_row] == "+",
                        "first", "last")


    ##transform if desired
    if (log_base > 0) {
        depth <- log(log_shift + depth, base=log_base)
    }


    ##compute PCA decomposition of vectors
    if (pc_sep) {
        pca_e <- prcomp(t(depth[exon_row, ]))
        pca_j <- prcomp(t(depth[junc_row, ]))

        basis_e <- pca_e$rotation[, 1:npc, drop=FALSE]
        basis_j <- pca_j$rotation[, 1:npc, drop=FALSE]
        pvar_e <- pca_e$sdev^2 / sum(pca_e$sdev^2)
        pvar_j <- pca_j$sdev^2 / sum(pca_j$sdev^2)
        pvar_ej <- paste(
            paste0("exon var expl = ", round(pvar_e, 3)),
            paste0("junc var expl = ", round(pvar_j, 3)), sep=", ")

    } else {
        ej_w <- ej_w / sum(ej_w)
        depth2 <- depth - matrix(rowMeans(depth), nrow=nrow(depth), ncol=ncol(depth))
        row_var <- apply(depth2, 1, sd)^2
        depth2[exon_row, ] <- depth2[exon_row, ] * sqrt(ej_w[1]/sum(row_var[exon_row]))
        depth2[junc_row, ] <- depth2[junc_row, ] * sqrt(ej_w[2]/sum(row_var[junc_row]))
        pca_ej <- prcomp(t(depth2))
        
        basis_e <- pca_ej$rotation[exon_row, 1:npc, drop=FALSE]
        basis_j <- pca_ej$rotation[junc_row, 1:npc, drop=FALSE]
        pvar_ej <- pca_ej$sdev^2 / sum(pca_ej$sdev^2)
        pvar_ej <- paste0("total var expl = ", round(pvar_ej, 3))
        
    }
    
    
    ##construct ggplot2 object
    gg_e <- data.frame(xmin=start(ranges(cc)[exon_row]),
                       xmax=end(ranges(cc)[exon_row]),
                       basis_e)
    gg_e <- melt(gg_e, id.vars=c("xmin", "xmax"))
    gg_e$ymin <- -1
    gg_e$ymax <- 1
    
    
    ##plot on genomic coordinates
    g_obj <- ggplot(gg_e, aes(xmin=xmin, xmax=xmax,
                              ymin=ymin, ymax=ymax,
                              color=value, fill=value))

    
    ##add basic plot structure
    g_obj <-
        g_obj + 
        geom_hline(yintercept=0, color=ifelse(use_blk, "#F0F0F0", "#3C3C3C")) +
        geom_rect() +
        facet_grid(variable~.) +
        scale_y_continuous(breaks=NULL, limits=c(-1, 3)) +
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
                     ymin=-1, ymax=1,
                     alpha=1, color="#3C3C3C", fill=NA)
    }

    
    ##add continuous color palette
    g_obj <- g_obj +
        scale_color_gradient2("expr", low="#053061", mid="#f7f7f7", high="#67001f", guide="none") +
        scale_fill_gradient2("expr", low="#053061", mid="#f7f7f7", high="#67001f")

    
    ##add splicing arrows to plot
    w_t <- diff(range(as.vector(ranges(cc))))
    for (j in 1:p_j) {
        for (ipc in 1:npc) {
            w_prop <- width(ranges(cc)[junc_row][j]) / w_t
            circle1 <- .pseudoArc(xmin=start(ranges(cc))[junc_row][j],
                                  xmax=end(ranges(cc))[junc_row][j],
                                  ymin=1, height=2*sqrt(w_prop))
            circle2 <- cbind(circle1, xmin=gg_e$xmin[1], xmax=gg_e$xmax[1],
                             ymin=gg_e$xmin[1], ymax=gg_e$ymax[1],
                             variable=paste0("PC", ipc), value=0)
            g_obj <- g_obj +
                geom_path(data=circle2, size=.75, aes(x=x, y=y),
                          color=.rgb2hex(crp2(basis_j[j, ipc]/2+.5)), 
                          arrow=arrow(length=unit(.015*npc, "npc"), ends=arrowhead[j]))
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
