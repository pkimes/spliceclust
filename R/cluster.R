#' Clustering method for concomp object
#'
#' Method returns cluster labels for \code{concomp} object 
#'
#' @param obj a \code{concomp} object
#' @param exon a logical whether to return labels based on exon expression
#'        (default = TRUE)
#' @param junc a logical whether to return labels based on junction expression
#'        (default = TRUE)
#' @param both a logical whether to return labels based on joint clustering
#'        (default = TRUE)
#' @param ej_w a numeric vector of length two specifying the relative sum of squares
#'        for exon and junctions to be used if both = TRUE (default = c(1, 1))
#' @param agg_f a function used to aggregate the expression for each sample
#'        to determine low vs. non-low expression (default = max)
#' @param min_cov a numeric value specifying the minimum threshold for agg_f of each
#'        sample to determine low vs. non-low expression (default = 1)
#' @param min_n an integer value specifying the minimum number of non-low samples
#'        needed for clustering (default = 10)
#' 
#' @return
#' cluster labels as vector of cluster indices with same length
#' as number of observations
#'
#' @name cluster
#' @export
#' @author Patrick Kimes
NULL

.cluster.concomp <- function(obj, exon = TRUE, junc = TRUE,
                             both = TRUE, ej_w = c(1, 1),
                             agg_f = max, min_cov = 1, min_n = 10) {

    if (!exon & !junc & !both)
        stop("at least one of exon, junc or both must be TRUE")

    
    ## unpack concomp
    vals_j <- juncValues(obj)
    vals_e <- exonValues(obj)


    if (ncol(vals_j) != ncol(vals_e))
        stop("number of samples fo juncValues and exonValues must match")

    
    ## initialize labels
    exon_labs <- rep(1, ncol(vals_e))
    junc_labs <- rep(1, ncol(vals_e))
    both_labs <- rep(1, ncol(vals_e))

    
    ## clustering with exons
    if (exon) {
        is_high <- (apply(vals_e, 2, agg_f) >= min_cov)
        if (sum(is_high) >= min_n) {
            modvals_e <- t(as.matrix(vals_e[, is_high, drop=FALSE]))
            modvals_e <- diag(mean(rowSums(modvals_e)) /
                                  rowSums(modvals_e)) %*% modvals_e
            modvals_e <- log10(modvals_e+1)
            modvals_e <- modvals_e - matrix(colMeans(modvals_e), byrow=TRUE,
                                            nrow=nrow(modvals_e),
                                            ncol=ncol(modvals_e))

            labs <- kmeans(modvals_e, 2, nstart=5)$cluster
            if (mean(labs == 1) < 0.5) { labs <- 3-labs }
            exon_labs[is_high] <- 1 + labs
        }
    }


    ## clustering with junctions
    if (junc) {
        is_high <- (apply(vals_j, 2, agg_f) >= min_cov)
        if (sum(is_high) >= min_n) {
            modvals_j <- t(as.matrix(vals_j[, is_high, drop=FALSE]))
            modvals_j <- diag(mean(rowSums(modvals_j)) /
                                  rowSums(modvals_j)) %*% modvals_j
            modvals_j  <- log10(modvals_j+1)
            modvals_j <- modvals_j - matrix(colMeans(modvals_j), byrow=TRUE,
                                            nrow=nrow(modvals_j),
                                            ncol=ncol(modvals_j))

            labs <- kmeans(modvals_j, 2, nstart=5)$cluster
            if (mean(labs == 1) < 0.5) { labs <- 3-labs }
            junc_labs[is_high] <- 1 + labs
        }
    }


    ## clustering with exon and junctions
    if (both) {
        ## ej_w <- ej_w / sum(ej_w)
        ## modvals_e <- vals_e - matrix(rowMeans(vals_e), nrow=nrow(vals_e), ncol=ncol(vals_e))
        ## modvals_j <- vals_j - matrix(rowMeans(vals_j), nrow=nrow(vals_j), ncol=ncol(vals_j))
        ## row_var_e <- apply(modvals_e, 1, sd)^2
        ## modvals_e <- modvals_e * sqrt(ej_w[1]/sum(row_var_e))
        ## row_var_j <- apply(modvals_j, 1, sd)^2
        ## modvals_j <- modvals_j * sqrt(ej_w[2]/sum(row_var_j))
        ## modvals_ej <- rbind(modvals_e, modvals_j)

        modvals_ej <- rbind(vals_e, vals_j)
        
        ## clustering with exon and junction
        is_high <- (apply(modvals_ej, 2, agg_f) >= min_cov)
        if (sum(is_high) >= min_n) {
            modvals_ej <- t(as.matrix(modvals_ej[, is_high, drop=FALSE]))
            modvals_ej <- diag(mean(rowSums(modvals_ej)) /
                                   rowSums(modvals_ej)) %*% modvals_ej
            modvals_ej <- log10(modvals_ej+1)
            modvals_ej <- modvals_ej - matrix(colMeans(modvals_ej), byrow=TRUE,
                                              nrow=nrow(modvals_ej),
                                              ncol=ncol(modvals_ej))

            labs <- kmeans(modvals_ej, 2, nstart=5)$cluster
            if (mean(labs == 1) < 0.5) { labs <- 3-labs }
            both_labs[is_high] <- 1 + labs
        }
        
    }
    
    ## report cluster labels
    list(exon_labs = exon_labs,
         junc_labs = junc_labs,
         both_labs = both_labs)
}



#' @rdname cluster
#' @aliases cluster,concomp-method
setMethod("cluster", signature(obj = "concomp"),
          function(obj, ...) .cluster.concomp(obj, ...))






.clustdiff <- function(obj, labs) {
    
    vals_e <- exonValues(obj)
    vals_j <- juncValues(obj)
    modvals_ej <- rbind(vals_e, vals_j)

    p_e <- nrow(vals_e)
    p_j <- nrow(vals_j)
    
    ## clustering with exon and junction
    modvals_ej <- t(as.matrix(modvals_ej[, labs > 1, drop=FALSE]))
    modvals_ej <- diag(mean(rowSums(modvals_ej)) /
                           rowSums(modvals_ej)) %*% modvals_ej
    modvals_ej <- log10(modvals_ej+1)
    modvals_ej <- modvals_ej - matrix(colMeans(modvals_ej), byrow=TRUE,
                                      nrow=nrow(modvals_ej),
                                      ncol=ncol(modvals_ej))

    sub_labs <- labs[labs > 1]
    vals_ej1 <- modvals_ej[sub_labs == 2, , drop=FALSE]
    vals_ej2 <- modvals_ej[sub_labs == 3, , drop=FALSE]
    
    load_ej <- colMeans(vals_ej1) - colMeans(vals_ej2)
    load_ej <- load_ej / sqrt(sum(load_ej^2))
    
    list(e_loads = matrix(load_ej[1:p_e], ncol=1),
         j_loads = matrix(load_ej[-(1:p_e)], ncol=1))
}
         
