#' SplicePCA: splicing graph PCA plot
#'
#' Splice graphs with PCA loadings shown along the exons and junctions to illustrate
#' expression patterns across a connected component. The function can be used to
#' illustrate, e.g. differences across groups/clusters or principal component
#' loadings. Note that PCA is performed on the log-transfored data by default.
#'
#' @param obj a \code{concomp} object with exon and junction information
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
#' @param txlist a GRangesList of transcripts or genes which should be queried and
#'        added to the plot if falling within the region of the connected component
#'        (default = NULL)
#' @param txdb a transcript database which can be used to query the transcript IDs
#'        identified from txlist (default = NULL)
#' @param orgdb a database that can be queried using keys obtained from \code{txdb}
#'        to determine corresponding gene symbols (default = NULL)
#' @param scores a logical whether to produce score plot rather than loadings plot
#'        (default = FALSE)
#' @param plot a logical whether to output a plot or the actually PCA analyses
#'        performed by \code{prcomp} (default = TRUE)
#' @param ... other parameters to be passed
#'
#' @details
#' \code{npc} cannot be larger than the rank of the centered data. When \code{pc_sep = TRUE},
#' this corresponds to \code{min(p_e, p_j, n-1)}, where \code{p_e}, \code{p_j}, \code{n}
#' correspond to the number of exons, junctions, and samples. When \code{pc_sep = TRUE},
#' this corresponds to \code{min(p_e+p_j, n-1)}.
#' 
#' @return
#' a ggplot2 plot showing the specified number of loadings
#'
#' @name splicepca
#' @export
#' @import ggplot2
#' @importFrom GGally ggpairs wrap
#' @author Patrick Kimes
NULL

.splicepca.concomp <- function(obj, npc = 3, pc_sep = TRUE, ej_w = c(1, 1),
                               log_base = 10, log_shift = 1,
                               genomic = TRUE, ex_use = 2/3,
                               flip_neg = TRUE, use_blk = FALSE,
                               txlist = NULL, txdb = NULL,
                               orgdb = NULL, scores = FALSE,
                               plot = TRUE, ...) {
    
    ##currently can't include gene models if not plotting on genomic scale
    if (!is.null(txlist) && !genomic)
        cat("ignoring txlist input since plotting on non-genomic scale. \n")

    
    ##unpack concomp
    vals_e <- eCoverage(obj)
    vals_j <- jCoverage(obj)
    
    ##dataset dimension
    p_e <- nrow(vals_e)
    p_j <- nrow(vals_j)
    n <- ncol(vals_e)

    ##check that npc is less than dim of data
    if (pc_sep == TRUE) {
        if (min(p_e, p_j, n-1) < 1) {
            warning("not enough dimension in data")
            return(NULL)
        } else if (npc > min(p_e, p_j, n-1)) {
            npc <- min(p_e, p_j, n-1)
            warning(paste("original npc larger than dim of data, using npc =", npc))
        }
    } else {
        if (min(p_e + p_j, n-1) < 1) {
            warning("not enough dimension in data")
            return(NULL)
        } else if (npc > min(p_e + p_j, n-1)) {
            npc <- min(p_e + p_j, n-1)
            warning(paste("original npc larger than dim of data, using npc =", npc))
        }
    }            

    
    ##transform if desired
    if (log_base > 0) {
        vals_e <- log(log_shift + vals_e, base=log_base)
        vals_j <- log(log_shift + vals_j, base=log_base)
    }


    ##compute PCA decomposition of vectors
    if (pc_sep) {
        pca_e <- prcomp(t(vals_e))
        pca_j <- prcomp(t(vals_j))

        e_loads <- pca_e$rotation[, 1:npc, drop=FALSE]
        j_loads <- pca_j$rotation[, 1:npc, drop=FALSE]

        e_scores <- pca_e$x[, 1:npc, drop=FALSE]
        j_scores <- pca_j$x[, 1:npc, drop=FALSE]
            
        ##flip junction scores/loadings if scores opposite of exon scores
        for (ipc in 1:npc) {
            if (cor(e_scores[, ipc], j_scores[, ipc]) < 0) {
                j_scores[, ipc] <- -j_scores[, ipc]
                j_loads[, ipc] <- -j_loads[, ipc]
            }
        }
        ej_scores <- data.frame(cbind(e_scores, j_scores))
        names(ej_scores) <- c(paste0("E_PC", 1:npc),
                              paste0("J_PC", 1:npc))
        
        ##generate text for % var expl
        pvar_e <- pca_e$sdev^2 / sum(pca_e$sdev^2)
        pvar_j <- pca_j$sdev^2 / sum(pca_j$sdev^2)
        pvar_e <- pvar_e[1:npc]
        pvar_j <- pvar_j[1:npc]
        if (scores) {
            pvar_ej <- c(
                paste0("EXON var expl = ", sprintf("%.2f", round(100*pvar_e, 2)), "%"),
                paste0("JUNC var expl = ", sprintf("%.2f", round(100*pvar_j, 2)), "%"))
        } else {
            pvar_ej <- paste(
                paste0("EXON var expl = ", sprintf("%.2f", round(100*pvar_e, 2)), "%"),
                paste0("JUNC var expl = ", sprintf("%.2f", round(100*pvar_j, 2)), "%"), sep="\n")
        }
            
        
    } else {
        ej_w <- ej_w / sum(ej_w)

        vals_e <- vals_e - matrix(rowMeans(vals_e), nrow=nrow(vals_e), ncol=ncol(vals_e))
        vals_j <- vals_j - matrix(rowMeans(vals_j), nrow=nrow(vals_j), ncol=ncol(vals_j))
        
        row_var_e <- apply(vals_e, 1, sd)^2
        vals_e <- vals_e * sqrt(ej_w[1]/sum(row_var_e))

        row_var_j <- apply(vals_j, 1, sd)^2
        vals_j <- vals_j * sqrt(ej_w[2]/sum(row_var_j))

        pca_ej <- prcomp(t(rbind(vals_e, vals_j)))
        
        e_loads <- pca_ej$rotation[1:p_e, 1:npc, drop=FALSE]
        j_loads <- pca_ej$rotation[-(1:p_e), 1:npc, drop=FALSE]

        ej_scores <- data.frame(pca_ej$x[, 1:npc, drop=FALSE])
        names(ej_scores) <- paste0("PC", 1:npc)
        
        pvar_ej <- pca_ej$sdev^2 / sum(pca_ej$sdev^2)
        pvar_ej <- paste0("TOTAL var expl = ", sprintf("%.2f", round(100*pvar_ej, 2)), "%")
        pvar_ej <- pvar_ej[1:npc]
    }
    
    
    ##return values if not interested in plots
    if (!plot) {
        if (pc_sep) {
            return(list(pca_e = pca_e, pca_j = pca_j))
        } else {
            return(pca_ej)
        }
    }
    
    
    ##create final scores or loadings plot
    if (scores) {
        ## create scores plot
        ## diag=list(continuous=wrap("densityDiag")) doesn't work accept params
        pl <- ggpairs(ej_scores, 
                      upper="blank", mapping=aes(color="type", fill="type"), 
                      lower=list(continuous=wrap("points", alpha=.1)),
                      diag=list(continuous=wrap("barDiag", alpha=.7,
                                    bins=30)))
        ## use bw theme
        pl <- pl + theme_bw()


        if (pc_sep) {
            ## change colors for exon PCA
            for (i in 1:npc) {
                for (j in 1:i) {
                    p <- pl[i, j]
                    p <- p + scale_color_manual(values="#08306b") +
                    scale_fill_manual(values="#08306b")
                    pl[i, j] <- p
                }
            }
            
            ## change colors for junction PCA
            for (i in (npc+1):(2*npc)) {
                for (j in (npc+1):i) {
                    p <- pl[i, j]
                    p <- p + scale_color_manual(values="#047760") +
                    scale_fill_manual(values="#047760")
                    pl[i, j] <- p
                }
            }

            ## change colors for cross PCA
            for (i in (npc+1):(2*npc)) {
                for (j in 1:npc) {
                    p <- pl[i, j]
                    p <- p + scale_color_manual(values="black") +
                        scale_fill_manual(values="black")
                    pl[i, j] <- p
                }
            }

        } else {
            ## change everything to black if plotting joint
            for (i in 1:ncol(ej_scores)) {
                for (j in 1:ncol(ej_scores)) {
                    p <- pl[i, j]
                    p <- p + scale_color_manual(values="black") +
                        scale_fill_manual(values="black")
                    pl[i, j] <- p
                }
            }
        }            

        ##return final plot
        pl
        
        
    } else {
        ## create splicing plot
        pl <- splicegralp(obj, e_loads, j_loads, load_lims = c(-1, 1), 
                          genomic = genomic, ex_use = ex_use,
                          flip_neg = flip_neg, use_blk = use_blk,
                          txlist = txlist, txdb = txdb, orgdb = orgdb)

        if (is(pl, "Tracks")) {
            
            ## add % var expl to each panel, determine where to add text
            xtxt <- 0
            if (flip_neg && all(strand(eRanges(obj)) == "-")) {
                xtxt <- max(pl@plot[[1]]$data$xmax)
            } else if (genomic) {
                xtxt <- min(pl@plot[[1]]$data$xmin)
            }

            ann_text <- data.frame(xmin = xtxt, xmax = xtxt, ymin = 0, ymax = 0,
                                   x = xtxt, y = 2.1, lab = pvar_ej, value = 0,
                                   variable = paste0("Dir", 1:length(pvar_ej)))

            new_geom <- geom_text(data = ann_text, 
                                  aes(x = x, y = y, label = lab),
                                  size = 3, family = "mono",
                                  color = ifelse(use_blk, "white", "black"),
                                  hjust = 0, vjust = 0)
            ## return final plot
            pl@plot[[1]] <- pl@plot[[1]] + new_geom
            pl@grobs[[1]] <- pl@grobs[[1]] + new_geom
            
        } else {
            ## add % var expl to each panel, determine where to add text
            xtxt <- 0
            if (flip_neg && all(strand(eRanges(obj)) == "-")) {
                xtxt <- max(pl$data$xmax)
            } else if (genomic) {
                xtxt <- min(pl$data$xmin)
            }

            ann_text <- data.frame(xmin = xtxt, xmax = xtxt, ymin = 0, ymax = 0,
                                   x = xtxt, y = 2.1, lab = pvar_ej, value = 0,
                                   variable = paste0("Dir", 1:length(pvar_ej)))

            ## return final plot
            pl <- pl + geom_text(data = ann_text, 
                                 aes(x = x, y = y, label = lab),
                                 size = 3, family = "mono",
                                 color = ifelse(use_blk, "white", "black"),
                                 hjust = 0, vjust = 0)
        }

        pl
    }
}

#' @keywords internal
#' @title splicepca method
#' @name splicepca
#' @aliases splicepca,concomp-method
setMethod("splicepca",
          signature(obj = "concomp"),
          function(obj, ...) .splicepca.concomp(obj, ...))

