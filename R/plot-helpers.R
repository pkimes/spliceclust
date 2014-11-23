

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


##helper for converting rgb triple to hex color
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


#' match concomp with gene or transcript name
#' 
#' Function matches concomp to nearest transcript and corresponding
#' gene names if both transcript list and organism db are provided. If
#' \code{orgdb} is not provided, the transcript IDs are returned rather than
#' the gene symbols
#'
#' @param obj a \code{concomp} object
#' @param txlist a list of transcripts, e.g. \code{exonsBy(txdb)}
#' @param txdb a transcript database, e.g. \code{TxDb.Hsapiens.UCSC.hg19.knownGene}
#' @param orgdb a organism database, e.g. \code{org.Hs.eg.db} (default = NULL)
#'
#' @return a vector of transcript or gene names
#'
#' @import GenomicRanges 
#' @export
#' @author Patrick Kimes
concomp2name <- function(obj, txlist, txdb, orgdb = NULL) {
    
    gr_e <- exons(obj)
    cand_idx <- unique(queryHits(findOverlaps(txlist, gr_e)))
    
    if (length(cand_idx) > 0 && !is.null(orgdb)) {
        gene_id <- select(txdb, keys=as.character(cand_idx),
                          columns="GENEID", keytype="TXID")$GENEID
        return(select(orgdb, keys=as.character(gene_id),
                      columns="SYMBOL", keytype="ENTREZID")$SYMBOL)
    } else {
        return(select(txdb, keys=as.character(cand_idx),
                      columns="TXNAME", keytype="TXID")$TXNAME)
    }
}

