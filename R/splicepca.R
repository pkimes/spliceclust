#' SplicePCA: splicing graph PCA plot
#'
#' Splice graphs with PCA loadings shown along the exons and junctions to illustrate
#' expression patterns across a connected component. The function can be used to
#' illustrate, e.g. differences across groups/clusters or principal component
#' loadings.
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
#' @name splicepca
#' @export
#' @import ggplot2 RColorBrewer reshape2 grid
#' @author Patrick Kimes
NULL

.splicepca.concomp <- function(obj, npc = 3, pc_sep = TRUE, ej_w = c(1, 1),
                               log_base = 10, log_shift = 1,
                               genomic = TRUE, ex_use = 2/3,
                               flip_neg = TRUE, use_blk = FALSE,
                               txlist = NULL, txdb = NULL,
                               orgdb = NULL, ...) {
    
    ##currently can't include gene models if not plotting on genomic scale
    if (!is.null(txlist) && !genomic)
        cat("ignoring txlist input since plotting on non-genomic scale. \n")

    
    ##unpack concomp
    vals_e <- exonValues(obj)
    vals_j <- juncValues(obj)
    
    ##dataset dimension
    p_e <- nrow(vals_e)
    p_j <- nrow(vals_j)


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

        ##flip junction loadings if scores opposite of exon scores
        for (ipc in 1:npc) {
            if (cor(pca_e$x[, ipc], pca_j$x[, ipc]) < 0)
                j_loads[, ipc] <- -j_loads[, ipc]
        }
        
        pvar_e <- pca_e$sdev^2 / sum(pca_e$sdev^2)
        pvar_j <- pca_j$sdev^2 / sum(pca_j$sdev^2)
        pvar_ej <- paste(
            paste0("exon var expl = ", round(pvar_e, 3)),
            paste0("junc var expl = ", round(pvar_j, 3)), sep=", ")

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

        pvar_ej <- pca_ej$sdev^2 / sum(pca_ej$sdev^2)
        pvar_ej <- paste0("total var expl = ", round(pvar_ej, 3))
        
    }
    
    splicegralp(obj, e_loads, j_loads, load_lims = c(-1, 1), 
                genomic = genomic, ex_use = ex_use,
                flip_neg = flip_neg, use_blk = use_blk,
                txlist = txlist, txdb = txdb, orgdb = orgdb)
}



#' @rdname splicepca
#' @aliases splicepca,concomp-method
setMethod("splicepca",
          signature(obj = "concomp"),
          function(obj, ...) .splicepca.concomp(obj, ...))
