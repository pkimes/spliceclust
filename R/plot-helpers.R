#' helper function for adding arcs to ggplot2 object via geom_path()
#'
#' modified from: http://stackoverflow.com/questions/6862742/
#'
#' @param xmin numeric value of left end of arc
#' @param xmax numeric value of right end of arc
#' @param ymin numeric value of bottom of arc
#' @param height numeric value of arc height
#' @param npoints integer number of points to draw arc
#'
#' @name pseudoArc
#' @keywords internal
#' @author Patrick Kimes
.pseudoArc <- function(xmin = 0, xmax = 10,
                       ymin = 0, height = 10, npoints = 100) {
    x_c  <- (xmax + xmin) / 2
    x_r <- (xmax - xmin) / 2
    tt <- seq(0, pi, length.out = npoints)
    
    xx <- x_c + x_r * cos(tt)
    yy <- ymin + height * sin(tt)
    return(data.frame(x = xx, y = yy))
}



#' helper for converting rgb triple to hex color
#'
#' @param rgb matrix with three columns with RGB values in each row
#'
#' @return vector of hex colors with length same as number of rows in
#'         input matrix
#'
#' @name rgb2hex
#' @keywords internal
#' @author Patrick Kimes
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
    strand(gr_e) <- "*"
    cand_idx <- unique(queryHits(findOverlaps(txlist, gr_e)))
    
    if (is.null(txlist)) {
        cand_id <- as.character(cand_idx)
    } else {
        cand_id <- names(txlist)[cand_idx]
    }
    
    if (is.null(txdb)) {
        return(data.frame(cand_idx, cand_id, stringsAsFactors=FALSE))
    } else {
        if (length(cand_idx) > 0) {
            if (is.null(orgdb)) {
                suppressMessages(
                    pair <- select(txdb, keys=cand_id,
                                   columns="TXNAME", keytype="TXID")$TXNAME
                    )
                return(data.frame(cand_idx, pair,
                                  stringsAsFactors=FALSE))
            } else {
                suppressMessages(
                    gene_id <- select(txdb, keys=cand_id,
                                      columns="GENEID", keytype="TXID")$GENEID
                    )
                pair <- rep(NA, length(gene_id))
                ## ugly way to handle ENTREZID and orgdb mismatch
                output <- tryCatch({
                    suppressMessages(
                        pair[!is.na(gene_id)] <- select(orgdb, keys=as.character(gene_id[!is.na(gene_id)]),
                                                        columns="SYMBOL", keytype="ENTREZID")$SYMBOL
                        )
                }, error = function(err) {
                    print("there was an error")
                    return(suppressMessages(select(txdb, keys=cand_id,
                                                   columns="TXNAME", keytype="TXID")$TXNAME))
                })
                if (all(is.na(pair))) { pair <- output }
                pair[is.na(pair)] <- "NA"
                return(data.frame(cand_idx, pair, stringsAsFactors=FALSE))
            }
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

    output <- data.frame("cc_id"=cc_matches, "tx_id"=tx_matches,
                         stringsAsFactors=FALSE)

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



#' match gene name (symbol) to connected component index
#' 
#' Function matches a character string gene symbol to overlapping connected components
#' according to annotations in \code{txdb} and \code{txlist}. Note to self: R follows
#' copy-on-write semantics.
#'
#' @param name a \code{character} string specifying the gene of interest
#' @param obj a \code{GRangesList} of connected component ranges, each name corresponding
#'        to the 'gene index'
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
name2gidx <- function(name, obj, txlist, txdb, orgdb) {
    ##identify gene models
    eid <- select(orgdb, keys=name, columns="ENTREZID", keytype="SYMBOL")
    txid <- select(txdb, keys=eid$ENTREZID, columns=c("TXID", "TXNAME"), keytype="GENEID")
    tx <- txlist[as.character(txid$TXID)]
    uniontx <- reduce(unlist(tx))

    ##compare to GRangesList of connect components
    matches <- subjectHits(findOverlaps(uniontx, obj))
    names(obj)[unique(matches)]
}



#' function defining useful vector of letters
#'
#' @keywords internal
#' @author Patrick Kimes
long_letters <- function() {
    ll <- LETTERS
    for (k in 2:5)
        ll <- c(ll, do.call("paste0", rep(list(LETTERS), k)))
    return(ll)
}



#' function defining useful colors
#'
#' @keywords internal
#' @author Patrick Kimes
plot_colors <- function() {
    crp <- colorRampPalette(c("#f7fbff", "#08306b"))
    crp2 <- colorRamp(c("#f7fbff", "#047760"))
    hl_cols <- c("#646464", "#CA4942", "#5A3589")
    crp_3col <- colorRampPalette(c("#053061", "#f7f7f7", "#67001f"))
    crp_3col2 <- colorRamp(c("#053061", "#f7f7f7", "#67001f"))
    return(list(col1 = crp, col2 = crp2, col3 = hl_cols,
                col4 = crp_3col, col5 = crp_3col2))
}


#' function defining all of the sort orders
#'
#' @param sort_idx integer value, see documentation for \code{splicegrahm}
#' @param vals_e matrix of exon coverages
#' @param vals_j matrix of junction coverages
#' @param n number of samples
#'
#' @details
#' sort_idx can take values of either:
#' \itemize{
#' \item{\code{1}}: sort on first exon
#' \item{\code{2}}: sort on PC 2
#' \item{\code{3}}: sort on mean exon coverage
#' \item{\code{4}}: sort on mean exon log-coverage
#' \item{\code{5}}: rev sort on mean exon coverage
#' \item{\code{6}}: rev sorton mean exon log-coverage
#' }
#' 
#' @keywords internal
#' @author Patrick Kimes
sampl_sort <- function(sort_idx, vals_e, vals_j, n) {
    if (length(sort_idx) == 1) {
        if (sort_idx == 1) {
            idx <- order(vals_e[1, ])
        } else if (sort_idx == 2) {
            pca <- prcomp(t(vals_e))
            idx <- order(pca$x[, 2])
        } else if (sort_idx == 3) {
            idx <- order(colMeans(vals_e))
        } else if (sort_idx == 4) {
            idx <- order(colMeans(log10(vals_e+1)))
        } else if (sort_idx == 5) {
            idx <- order(colMeans(vals_e), decreasing=TRUE)
        } else if (sort_idx == 6) {
            idx <- order(colMeans(log10(vals_e+1)), decreasing=TRUE)
        } else {
            idx <- 1:n
        }
    } else if (length(sort_idx) == n) {
        idx <- sort_idx
    } else {
        idx <- 1:n
    }
    return(idx)
}



#' adjust ranges if necessary
#'
#' @param gr_e \code{GenomicRanges} for exons
#' @param gr_j \code{GenomicRanges} for junctions
#' @param tx_plot output form \code{find_annotations}
#' @param ex_use see \code{splicegrahm} documentation
#' @param gr_base \code{GenomicRanges} to use for adjusting gr_e, gr_j,
#'        if NULL, then just adjust based on \code{gr_e}
#'        (default = NULL)
#'
#' @keywords internal
#' @author Patrick Kimes
adj_ranges <- function(gr_e, gr_j, tx_plot, ex_use, gr_base = NULL) {
    annot_track <- NULL

    ##use gr_e if gr_base not provided
    if (is.null(gr_base)) { gr_base <- gr_e } 

    ##create 'union model' based on txlist and gr_base
    if (is.null(tx_plot)) {
        models <- gr_base
    } else {
        bb <- tx_plot
        mcols(bb) <- NULL
        models <- c(gr_base, bb)
    }
    strand(models) <- "*"
    models <- reduce(models)

    dna_len <- width(range(models))
    rna_len <- sum(width(ranges(models)))
    base_shift <- min(start(models))
    
    ##don't squish exons if not needed (should make shrink = 0)
    if (rna_len/dna_len >= ex_use) {
        if(!is.null(tx_plot)) {
            annot_track <- ggplot() +
                geom_alignment(tx_plot, gap.geom="arrow", aes(group=tx)) +
                    theme_bw()
        }
        return(list(gr_e = gr_e, gr_j = gr_j,
                    annot_track = annot_track,
                    genomic = TRUE))
    }
    
    ##shrinkage factor for introns
    shrink <- 1 - (1-ex_use)/ex_use * rna_len/(dna_len - rna_len)
    ir_i <- gaps(ranges(models))
    preintrons <- function(x) { sum(width(restrict(ir_i, end=x))) }
    
    if (!is.nan(shrink)) {
        ##shift exons and junctions
        e_shifts <- sapply(start(gr_e), preintrons) * shrink
        start(gr_e) <- start(gr_e) - e_shifts
        end(gr_e) <- end(gr_e) - e_shifts
        gr_e <- shift(gr_e, -base_shift+1)

        if (length(gr_j) > 0) {
            j_shifts1 <- sapply(start(gr_j), preintrons) * shrink
            j_shifts2 <- sapply(end(gr_j), preintrons) * shrink
            start(gr_j) <- start(gr_j) - j_shifts1
            end(gr_j) <- end(gr_j) - j_shifts2
            gr_j <- shift(gr_j, -base_shift+1)
        }
    }
    
    ##shift annotation models
    if (!is.null(tx_plot)) {
        if (!is.nan(shrink)) {
            m_shifts1 <- sapply(start(tx_plot), preintrons) * shrink
            m_shifts2 <- sapply(end(tx_plot), preintrons) * shrink
            start(tx_plot) <- start(tx_plot) - m_shifts1
            end(tx_plot) <- end(tx_plot) - m_shifts2
            tx_plot <- shift(tx_plot, -base_shift+1)
            
            annot_track <- ggplot() +
                geom_alignment(tx_plot, gap.geom="arrow", aes(group=tx)) +
                    theme_bw()
        }
    }
    
    ## return gr_e, gr_j, annot_track
    list(gr_e = gr_e, gr_j = gr_j, annot_track = annot_track, genomic = FALSE)
}



#' construct annotation track for concomp object
#'
#' @param obj \code{concomp} object
#' @param txlist see \code{splicegrahm} documentation
#' @param txdb see \code{splicegrahm} documentation
#' @param orgdb see \code{splicegrahm} documentation
#' @param eps see \code{splicegrahm} documentation
#'
#' @return
#' list containing a \code{GRangesList}, \code{tx_match}, and a
#' \code{ggplot} object, \code{annot_track}
#' 
#' @keywords internal
#' @author Patrick Kimes
find_annotations <- function(obj, txlist, txdb, orgdb, eps) {
    tx_plot <- NULL

    cands <- concomp2name(obj, txlist, txdb, orgdb)
    cand_idx <- cands[, 1]
    cand_names <- cands[, 2]

    if (length(cand_idx) > 0) {
        tx_match <- txlist[cand_idx]
        names(tx_match) <- make.unique(cand_names, sep = "_")
        
        ##take tx_match and compare if lies outside
        if (!is.null(eps)) {
            tx_ranges <- unlist(range(tx_match))
            obj_range <- range(exons(obj))
            tx_match <- tx_match[start(tx_ranges) > start(obj_range) - eps &
                                     end(tx_ranges) < end(obj_range) + eps]
        }
        
        ## unlist models to allow arrows in ggbio plot (bug?)
        if (length(tx_match) > 0) {
            tx_len <- sapply(tx_match, length)
            tx_plot <- unlist(tx_match)
            mcols(tx_plot)$tx <- rep(names(tx_match), times=tx_len)
        }
    }

    tx_plot
}





##
## alternatively, can shrink using e.g. maxgap = 100, with
## shrink.fun <- shrinkageFun(gaps(gr_e), max.gap = 100)
## new_gr_e <- shrink.fun(gr_e)
##  shrinkageFun uses findIntervals(x, y) to determine how many introns are to the left - clever
##

## alternative approach for find_annotations
##
## cands <- concomp2name(obj, txlist, txdb, orgdb)
## cand_idx <- cands[, 1]
## cand_names <- unique(cands[, 2])[1]
## eval(bquote(
##     annot_track <- ggplot(txdb) +
##         geom_alignment(which = genesymbol[.(cand_names)]) +
##             theme_bw()
##     ))
## list(tx_match = tx_match, annot_track = annot_track)



