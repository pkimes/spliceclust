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
    cand_idx <- unique(queryHits(findOverlaps(txlist, gr_e)))
    
    if (is.null(txdb)) {
        return(cbind(cand_idx,
                     as.character(cand_idx)))
    } else {
        if (length(cand_idx) > 0) {
            if (is.null(orgdb)) {
                pair <- select(txdb, keys=as.character(cand_idx),
                               columns="TXNAME", keytype="TXID")$TXNAME
                return(cbind(cand_idx, pair))
            } else {
                gene_id <- select(txdb, keys=as.character(cand_idx),
                                  columns="GENEID", keytype="TXID")$GENEID
                pair <- rep(NA, length(gene_id))
                pair[!is.na(gene_id)] <- select(orgdb, keys=as.character(gene_id[!is.na(gene_id)]),
                                                columns="SYMBOL", keytype="ENTREZID")$SYMBOL
                pair[is.na(pair)] <- "NA"
                return(cbind(cand_idx, pair))
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
    return(list(col1 = crp, col2 = crp2, col3 = hl_cols))
}



#' adjust ranges if necessary
#'
#' @param gr_e \code{GenomicRanges} for exons
#' @param gr_j \code{GenomicRanges} for junctions
#' @param dna_len total range of \code{gr_e}
#' @param rna_len total length of \code{gr_e}
#' @param ex_use see \code{splicegrahm} documentation
#' @param p_e number of exons in \code{gr_e}
#'
#' @keywords internal
#' @author Patrick Kimes
adj_ranges <- function(gr_e, gr_j, dna_len, rna_len, ex_use, p_e) {
    shift_1 <- -start(ranges(gr_e))[1]+1
    ranges(gr_e) <- shift(ranges(gr_e), shift_1)
    ranges(gr_j) <- shift(ranges(gr_j), shift_1)
    
    shrink_by <- 1 - (1 - ex_use)/ex_use * rna_len/(dna_len - rna_len)
    gps <- distance(ranges(gr_e)[-p_e], ranges(gr_e)[-1]) * shrink_by
    
    eeb <- start(ranges(gr_e))[-1]
    for (ii in (p_e-1):1) {
        i_adj <- start(ranges(gr_e)) >= eeb[ii]
        start(ranges(gr_e))[i_adj] <- start(ranges(gr_e))[i_adj] - gps[ii]
        
        i_adj <- end(ranges(gr_e)) >= eeb[ii]
        end(ranges(gr_e))[i_adj] <- end(ranges(gr_e))[i_adj] - gps[ii]
        
        i_adj <- start(ranges(gr_j)) >= eeb[ii]
        start(ranges(gr_j))[i_adj] <- start(ranges(gr_j))[i_adj] - gps[ii]
        
        i_adj <- end(ranges(gr_j)) >= eeb[ii]
        end(ranges(gr_j))[i_adj] <- end(ranges(gr_j))[i_adj] - gps[ii]
    }
    list(gr_e = gr_e, gr_j = gr_j)
}



#' function defining all of the sort orders
#'
#' @param sort_idx integer value, see documentation for \code{splicegrahm}
#' @param vals_e matrix of exon coverages
#' @param vals_j matrix of junction coverages
#' 
#' @keywords internal
#' @author Patrick Kimes
sampl_sort <- function(sort_idx, vals_e, vals_j) {
    if (length(sort_idx) == 1) {
        if (sort_idx == 1) {
            idx <- order(vals_e[1, ])
        } else if (sort_idx == 2) {
            pca <- prcomp(t(vals_e))
            idx <- order(pca$x[, 2])
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


#' construct annotation track for concomp object
#'
#' @param obj \code{concomp} object
#' @param txlist see \code{splicegrahm} documentation
#' @param txdb see \code{splicegrahm} documentation
#' @param orgdb see \code{splicegrahm} documentation
#'
#' @return
#' list containing a \code{GRangesList}, \code{tx_match}, and a
#' \code{ggplot} object, \code{annot_track}
#' 
#' @keywords internal
#' @author Patrick Kimes
find_annotations <- function(obj, txlist, txdb, orgdb) {
    tx_match <- NULL
    annot_track <- NULL

    cands <- concomp2name(obj, txlist, txdb, orgdb)
    cand_idx <- cands[, 1]
    cand_names <- cands[, 2]

    if (length(cand_idx) > 0) {
        tx_match <- txlist[cand_idx]
        names(tx_match) <- make.unique(cand_names)
        annot_track <- ggplot(tx_match) + geom_alignment() + theme_bw()            
    }
    list(tx_match = tx_match, annot_track = annot_track)
}


#' construct ggplot2 object for plotting
#'
#' @param gr_e \code{GenomicRanges} for exons
#' @param gr_j \code{GenomicRanges} for junctions
#' @param vals_e matrix of exon coverages
#' @param vals_j matrix of junction coverages
#' @param n number of samples in \code{vals_e}, \code{vals_j}
#' @param p_j number of junctions in \code{gr_j}
#' @param j_incl see \code{splicegrahm} documentation
#' @param log_base see \code{splicegrahm} documentation
#' @param log_shift see \code{splicegrahm} documentation
#' @param bin see \code{splicegrahm} documentation
#'
#' @keywords internal
#' @author Patrick Kimes
sg_create <- function(gr_e, gr_j, vals_e, vals_j, j_incl,
                      log_base, log_shift, bin, n, p_j) {

    ## if (!(loc %in% c(-1, 0, 1)))
    ##     stop("loc value is invalid, must be 1, 0, -1")
    
    gg_e <- data.frame(xmin=start(ranges(gr_e)),
                       xmax=end(ranges(gr_e)),
                       vals_e)
    gg_e <- reshape2::melt(gg_e, id.vars=c("xmin", "xmax"))
    gg_e$ymin <- as.numeric(gg_e$variable)
    gg_e$ymax <- gg_e$ymin + 1
    
    ##transform if desired
    if (log_base > 0) {
        gg_e$value <- log(log_shift + gg_e$value, base=log_base)
    }
    
    ##perform binning
    if (bin && log_base > 0) {
        gg_e$value <- factor(floor(gg_e$value))
    }
    
    ##include junctions if necessary
    if (j_incl) {
        junc_x <- seq(min(min(ranges(gr_e))),
                      max(max(ranges(gr_e))),
                      length.out=p_j+2)
        junc_x <- junc_x[2:(p_j+1)]
        
        junc_y <- 2.25*n
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
            gg_j$value <- factor(floor(gg_j$value))
        }

        gg_e$kind <- "e"
        gg_e <- rbind(gg_e, gg_j)
    }

    return(gg_e)
}



#' construct base ggplot2 object for splicegrahm plot
#'
#' @param sg_df data.frame output from \code{sg_create}
#' @param gr_e \code{GenomicRanges} for exons
#' @param n number of samples in \code{vals_e}, \code{vals_j}
#' @param use_blk see \code{splicegrahm} documentation
#' @param j_incl see \code{splicegrahm} documentation
#' @param genomic see \code{splicegrahm} documentation
#' @param log_base see \code{splicegrahm} documentation
#' @param bin see \code{splicegrahm} documentation
#' 
#' @keywords internal
#' @author Patrick Kimes
sg_drawbase <- function(sg_df, use_blk, j_incl, genomic, gr_e,
                        log_base, bin, n) {
    
    ##color generators
    pal <- plot_colors()

    sg_obj <- ggplot(sg_df, aes(xmin=xmin, xmax=xmax,
                              ymin=ymin, ymax=ymax,
                              color=value, fill=value))

    ##add horizontal line first
    sg_obj <- sg_obj + 
        geom_hline(yintercept=n/2,
                   color=ifelse(use_blk, "#F0F0F0", "#3C3C3C"))
    
    ##add basic plot structure
    sg_obj <- sg_obj + 
        geom_rect() +
        scale_y_continuous(breaks=NULL, limits=c(0, (2.15+j_incl)*n)) +
        ylab("") +
        xlab(ifelse(genomic, paste0("Genomic Coordinates, ", seqnames(gr_e[1])),
                    "non-genomic coordinates")) +
        theme_bw()

    ##frame exons if not using black background
    if (use_blk) {
        sg_obj <- sg_obj +
            theme(panel.grid.minor.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.background = element_rect(fill="#3C3C3C"))
    } else {
        sg_obj <- sg_obj +
            annotate("rect", size = .25,
                     xmin = start(ranges(gr_e)) - .5,
                     xmax = end(ranges(gr_e)) + .5,
                     ymin = 1, ymax = n+1,
                     alpha = 1, color = "#3C3C3C", fill = NA)
    }

    ##add continuous or discrete color palette
    if (log_base > 0 && bin) {
        v_max <- max(as.numeric(as.character(sg_df$value)))
        sg_obj <- sg_obj + 
            scale_color_manual("expr", breaks=0:v_max, values=pal$col1(v_max+1), guide="none") +
            scale_fill_manual("expr", breaks=0:v_max, values=pal$col1(v_max+1),
                              labels=paste0("<", log_base^(1:(v_max+1))))
    } else {
        sg_obj <- sg_obj +
            scale_color_continuous("expr", low="#f7fbff", high="#08306b", guide="none") +
            scale_fill_continuous("expr", low="#f7fbff", high="#08306b")
    }

    sg_obj
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
#' 
#' @keywords internal
#' @author Patrick Kimes
sg_drawjuncs <- function(sg_obj, sg_df, j_incl, use_blk, iflip,
                         gr_e, gr_j, vals_j, n, p_j) {

    ##strand of junctions for arrow heads
    arrowhead <- ifelse(as.character(strand(gr_j)) == "-", "last", "first")
    
    ##color generators
    pal <- plot_colors()

    ##long list of letters
    ll_LETTERS <- long_letters()
    
    ##add splicing arrows to plot
    e_prop <- rowMeans(vals_j > 0) 
    for (j in 1:p_j) {
        w_prop <- width(ranges(gr_j)[j]) / width(range(gr_e))
        circle1 <- .pseudoArc(xmin=start(ranges(gr_j))[j],
                              xmax=end(ranges(gr_j))[j],
                              ymin=n+1, height=1*n*sqrt(w_prop))

        ##only include arrows if direction is known
        if (all(strand(gr_e) == "*")) {
            sg_obj <- sg_obj +
                annotate("path", size=.75,
                         x=circle1$x, y=circle1$y,
                         color=.rgb2hex(pal$col2(e_prop[j])))
        } else {
            sg_obj <- sg_obj +
                annotate("path", size=.75,
                         x=circle1$x, y=circle1$y,
                         color=.rgb2hex(pal$col2(e_prop[j])),
                         arrow=grid::arrow(length=grid::unit(.015, "npc"), ends=arrowhead[j]))
        }
    }

    ##add text labels for splicing arrows
    if (j_incl) {
        junc_x <- seq(min(min(ranges(gr_e))),
                      max(max(ranges(gr_e))),
                      length.out=p_j+2)
        junc_x <- junc_x[2:(p_j+1)]
        junc_y <- 2.25*n
        s_size <- .5
        junc_w <- width(range(gr_e)) / (2.5*(max(p_j, 5)+1))
    
        abc <- ll_LETTERS[1:p_j]
        if (iflip) { abc <- rev(abc) }
        w_prop <- width(ranges(gr_j)) / width(range(gr_e))

        sg_obj <- sg_obj +
            annotate("text", size=3,
                     x=(start(ranges(gr_j)) + end(ranges(gr_j))) / 2,
                     y=(n+1)*(1 + sqrt(w_prop)), vjust=0, 
                     label=abc,
                     color=ifelse(use_blk, "#F0F0F0", "#3C3C3C"))

        sg_obj <- sg_obj +
            annotate("text", size=3,
                     x=junc_x, y=rep(junc_y + (n+1)*(.025 + s_size), p_j),
                     label=abc, vjust=0,
                     color=ifelse(use_blk, "#F0F0F0", "#3C3C3C"))

        ##add rectangles around 
        if (!use_blk) {
            sg_obj <- sg_obj +
                annotate("rect", size=.25,
                         xmin=junc_x - junc_w - 1, xmax=junc_x + junc_w + 1,
                         ymin=junc_y + s_size, ymax=junc_y + (n+1)*s_size,
                         alpha=1, color=.rgb2hex(pal$col2(1)), fill=NA)
        }   
    }

    sg_obj
}


## fix 1 to 2 groups in a single plot
## -- pre-parsing of data ...
## -- annotations: 
## -- sg_create: create sg_df separately and then combine
## -- sg_drawbase: 
## -- sg_drawjuncs: 
