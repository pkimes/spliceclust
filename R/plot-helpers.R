

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
#' Function matches \code{\link[=concomp-constructor]{concomp}} to nearest
#' transcript using \code{txdb} and \code{txlist}. This is further mapepd to
#' gene names if \code{orgdb} is specified. If \code{orgdb} is not provided,
#' transcript IDs are returned instead of gene symbols.
#'
#' @param obj a \code{\link[=concomp-constructor]{concomp}} object
#' @param txlist a list of transcripts, e.g. \code{exonsBy(txdb)}
#' @param txdb a transcript database, e.g. \code{TxDb.Hsapiens.UCSC.hg19.knownGene}
#'        (default = NULL)
#' @param orgdb a organism database, e.g. \code{org.Hs.eg.db} (default = NULL)
#'
#' @return a vector of transcript or gene names
#'
#' @import GenomicRanges 
#' @export
#' @author Patrick Kimes
concomp2name <- function(obj, txlist, txdb = NULL, orgdb = NULL) {
    
    gr_e <- exons(obj)
    cand_idx <- unique(queryHits(findOverlaps(txlist, gr_e)))
    
    if (is.null(txdb)) {
        return(cbind(cand_idx,
                     as.character(cand_idx)))
    } else {
        if (length(cand_idx) > 0 && !is.null(orgdb)) {
            gene_id <- select(txdb, keys=as.character(cand_idx),
                              columns="GENEID", keytype="TXID")$GENEID
            pair <- rep(NA, length(gene_id))
            pair[!is.na(gene_id)] <- select(orgdb, keys=as.character(gene_id[!is.na(gene_id)]),
                                            columns="SYMBOL", keytype="ENTREZID")$SYMBOL
            pair[is.na(pair)] <- "NA"
            return(cbind(cand_idx, pair))
    
        } else {
            pair <- select(txdb, keys=cand_idx,
                           columns="TXNAME", keytype="TXID")$TXNAME
            return(cbind(cand_idx, pair))
        }
    }
}


#' match matrix of diffsplice output with gene or transcript name
#' 
#' Function matches \code{data.frame} of diffsplice output to nearest
#' transcripts using \code{txdb} and \code{txlist}. This is further mapped to
#' gene names if \code{orgdb} is specified. If \code{orgdb} is not provided,
#' transcript IDs are returned instead of gene symbols.
#'
#' @param obj a \code{data.frame} loaded from \code{\link{readchr}}
#' @param ichr a numeric value specifying the chromosome number, e.g. 9
#' @param seqlen a numeric value specifying the sequence length corresponding to the
#'        chromosome specified in \code{ichr}, e.g. take from \code{seqlen(txdb)}
#' @param txlist a list of transcripts, e.g. \code{exonsBy(txdb)}
#' @param txdb a transcript database, e.g. \code{TxDb.Hsapiens.UCSC.hg19.knownGene}
#'        (default = NULL)
#' @param orgdb a organism database, e.g. \code{org.Hs.eg.db} (default = NULL)
#'
#' @return a \code{data.frame} with connected component IDs and overlapping
#'         transcript or gene names
#'
#' @import GenomicRanges
#' @export
#' @author Patrick Kimes
diffsplice2name <- function(obj, ichr, seqlen, txlist, txdb = NULL, orgdb = NULL) {
    
    chr <- obj[obj$kind == "e", ]
    chr$chr <- paste0("chr", ichr)
    
    chr_gl <- makeGRangesFromDataFrame(chr)
    seqlengths(chr_gl) <- seqlen
    
    chr_gl <- split(chr_gl, chr$gIdx)

    matches <- findOverlaps(txlist, chr_gl)
    cc_matches <- names(chr_gl)[subjectHits(matches)]
    tx_matches <- as.character(queryHits(matches))

    output <- data.frame("cc_id"=cc_matches,
                         "tx_id"=tx_matches)

    if (!is.null(txdb)) {
        tx_nm <- select(txdb, keys=tx_matches,
                        columns="TXNAME", keytype="TXID")$TXNAME
        output$tx_nm <- tx_nm

        if (length(tx_matches) > 0 && !is.null(orgdb)) {
            gene_id <- select(txdb, keys=tx_matches,
                              columns="GENEID", keytype="TXID")$GENEID
            pair <- rep(NA, length(gene_id))
            pair[!is.na(gene_id)] <- select(orgdb, keys=as.character(gene_id[!is.na(gene_id)]),
                                            columns="SYMBOL", keytype="ENTREZID")$SYMBOL
            output$gn_id <- gene_id
            output$gn_nm <- pair
        }
    }
    
    return(output)
}
