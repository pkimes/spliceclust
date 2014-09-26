#' ---
#' title: "Chr9 Analysis"
#' output:
#'    html_document:
#'        toc: true
#'        fig_width: 2
#'        fig_height: 6
#'        fig_caption: true
#' ---

#+ echo=FALSE
knitr::opts_chunk$set(collapse=TRUE, comment="##")

    
#' ------------------------------------------------------------------------
#' ## Preliminary
#' 

#' load necessary packages
#+ hide=TRUE
library("ggplot2")
library("ggbio")
library("RColorBrewer")
library("GenomicRanges")
library("SplicingGraphs")
library("annotate")
library("BSgenome.Hsapiens.UCSC.hg19")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("org.Hs.eg.db")
library("PTAk")
library("igraph")
library("mvnmle")


#' set nicer colors
brew_set1 <- brewer.pal(9, "Set1")
brew_pastel2 <- brewer.pal(8, "Pastel2")
brew_accent <- brewer.pal(8, "Accent")
brew_greys <- brewer.pal(9, "Greys")


#' specify project directory and source functions
root <-"/Users/pkimes/Dropbox/Git/spliceclust/"
source(paste0(root, "R/readchr.R"))


#' load chromesome 9 data for 177 lusc samples
## complete datasest
chr9 <- readchr(paste0(root, "data/lusc/chr9_gene.txt"), 177)
## exon only dataset
chr9_e <- subset(chr9, kind == "e")



#' ------------------------------------------------------------------------
#' ## general chr9
#'

#' count up number of exons
exon_cnts <- table(chr9_e$gIdx)
summary(as.numeric(exon_cnts))


#' plot distribution of exon counts
#+ hide=TRUE, warning=FALSE, echo=FALSE
qplot(x=Freq, data=data.frame(exon_cnts), binwidth=I(1),
      color=I(brew_greys[5]), fill=I(brew_pastel2[1])) +
    theme_bw() +
    scale_y_continuous(trans=scales::log1p_trans(),
                       breaks=c(0, 10^(0:5))) +
    xlab("number of exons") +
    ylab("log(1 + #models) , real counts shown") +
    ggtitle("distribution of exon counts for connected components")



#' ------------------------------------------------------------------------
#' ## single exon genes
#'

#' subset on only single gene models
e1_gIdx <- names(which(exon_cnts == 1))
e1set <- subset(chr9_e, gIdx %in% e1_gIdx)
e1set$chr <- "chr9"


#' convert to `GRanges` object to use with `ggbio`
e1set_gr <- makeGRangesFromDataFrame(e1set,
                                     seqnames.field="chr",
                                     keep.extra.columns=TRUE)
seqlengths(e1set_gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)["chr9"]


#' plot location of single exon genes on chromosome
#+ hide=TRUE, warning=FALSE, echo=FALSE, fig.height=2, fig.width=10
## load hg19 cytoband information
data(ideoCyto, package = "biovizBase")
biovizBase::isIdeogram(ideoCyto$hg19)

##construct ideogram plot for comparison
gg_chr9 <- autoplot(subset(ideoCyto$hg19, seqnames(ideoCyto$hg19) == "chr9"),
                    layout="karyogram", cytoband=TRUE) +
                        guides(fill=FALSE)
fixed(gg_chr9) <- TRUE

##construct karyogram with chr19 with 
gg_e1 <- autoplot(e1set_gr, layout="karyogram", alpha=1/100)
fixed(gg_e1) <- TRUE

##combine plots
tracks("chr" = gg_chr9, "e1 genes" = gg_e1,
       heights=c(.5, .5), xlab="chr9 positions") +
    theme_tracks_sunset()


#' compute the average expression
#+ hide=TRUE, warning=FALSE, echo=FALSE
e1_exp <- as.matrix(mcols(e1set_gr)[paste0("s", 1:177)])
e1_expstat <- data.frame("med"=apply(e1_exp, 1, median),
                         "mean"=rowMeans(e1_exp),
                         "sd"=apply(e1_exp, 1, sd),
                         "mad"=apply(e1_exp, 1, mad))
qplot(x=med, y=mad,
      data=e1_expstat,
      color=I(brew_set1[2]), alpha=I(1/5)) +
    scale_y_continuous(trans=scales::log1p_trans(),
                       breaks=c(0, 10^(0:6))) +
    scale_x_continuous(trans=scales::log1p_trans(),
                       breaks=c(0, 10^(0:6))) +
    theme_bw() +
    xlab("median coverage across 177 samples") +
    ylab("MAD") +
    ggtitle("distribution of expression for components with single exon")


#' compute distribution of  exon length
#+ hide=TRUE, warning=FALSE, echo=FALSE
e1_widths <- data.frame("width"=width(ranges(e1set_gr)))
qplot(x=width, data=e1_widths,
      color=I(brew_greys[5]), fill=I(brew_pastel2[1])) +
    theme_bw() +
    scale_y_continuous(trans=scales::log1p_trans(),
                       breaks=c(0, 10^(0:5))) +
    xlab("length of exon (bp)") +
    ylab("log(1 + #models) , real counts shown") +
    ggtitle("distribution of exon lengths for components with single exon")


#' load transcript annotations and only look at chr9
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
isActiveSeq(txdb)[seqlevels(txdb)] <- FALSE
isActiveSeq(txdb)["chr9"] <- TRUE

##want to work at exon level aggregated by genes or transcripts
exbygene <- exonsBy(txdb, "gene")
exbytx <- exonsBy(txdb, "tx")


#' Calculate overlap with genes and transcripts
## overlap with genes
hits_g <- findOverlaps(exbygene, e1set_gr)
ucsc_ovlp_g <- select(org.Hs.eg.db, keys=names(exbygene[queryHits(hits_g)]),
                      columns="SYMBOL", keytype="ENTREZID")

## overlap with transcripts
hits_tx <- findOverlaps(exbytx, e1set_gr)
ucsc_ovlp_tx <- select(txdb, keys=names(exbytx[queryHits(hits_tx)]),
                       columns=c("GENEID", "TXNAME"), keytype="TXID")
ucsc_ovlp_tx$idx <- 1:nrow(ucsc_ovlp_tx)
temp <- select(org.Hs.eg.db, keys=unique(ucsc_ovlp_tx$GENEID),
               columns="SYMBOL", keytype="ENTREZID")
names(temp) <- c("GENEID", "SYMBOL")
ucsc_ovlp_tx <- merge(ucsc_ovlp_tx, temp, all.x=TRUE, sort=FALSE)
ucsc_ovlp_tx <- ucsc_ovlp_tx[order(ucsc_ovlp_tx$idx), ]

## 1. add distance to closest named genes
## 2. look at named genes for highest expression 1e genes


#' There were `r length(unique(queryHits(hits_g)))` unique UCSC genes that overlapped with
#' `r length(unique(subjectHits(hits_g)))` unique connected components. Total number of
#' overlaps was `r length(queryHits(hits_g))`.
#' 
#+ hide=TRUE, warning=FALSE, echo=FALSE, fig.width=10, fig.height=6
ovlp_table_g <- data.frame(table(ucsc_ovlp_g$SYMBOL))
ggplot(data=ovlp_table_g) +
    geom_text(aes(label=Var1, x=1:length(Var1), y=Freq)) +
    theme_bw() +
    xlab("") +
    ylab("# overlapping components") +
    ggtitle("UCSC KnownGenes overlapping single gene connected components")


#' Similar plot as above, but at the scale of transcripts.
#' 
#+ hide=TRUE, warning=FALSE, echo=FALSE, fig.width=10, fig.height=6
ovlp_table_tx <- data.frame(table(ucsc_ovlp_tx$SYMBOL))
ggplot(data=ovlp_table_tx) +
    geom_text(aes(label=Var1, x=1:length(Var1), y=Freq)) +
    theme_bw() +
    xlab("") +
    ylab("# overlapping components") +
    ggtitle("UCSC KnownGene transcripts overlapping single gene connected components")


#' We now add UCSC KnownGenes to above `tracks` plots. Note that these are not
#' real transcripts, just the 'union' transcripts for each gene constructed by
#' joining all reported isoforms. (note: need to specify `fixed() <- TRUE` for
#' ucsc track since plot includes ideogram.)
#+ hide=TRUE, warning=FALSE, echo=FALSE, fig.height=3, fig.width=10
gg_ucsc <- autoplot(exbygene)
fixed(gg_ucsc) <- TRUE

tracks("chr" = gg_chr9,
       "e1 genes" = gg_e1,
       "UCSC KG" = gg_ucsc,
       heights=c(1/3, 1/3, 1/3), xlab="chr9 positions") +
    theme_tracks_sunset()


#' Next, we focus on `r ovlp_table_g$Var1[which.max(ovlp_table_g$Freq)]`  which
#' overlapped with `r max(ovlp_table_g$Freq)` single exon connected components.
#+ fig.height=2, fig.width=10
top1 <- ovlp_table_g$Var1[which.max(ovlp_table_g$Freq)]
top1_tx <- subset(ucsc_ovlp_tx, ucsc_ovlp_tx$SYMBOL == top1)
table(top1_tx$TXNAME)

top1_txname <- names(which.max(table(top1_tx$TXNAME)))
top1_hits <- subset(hits_tx, ucsc_ovlp_tx$TXNAME == top1_txname)

gg_chr9zoom <- autoplot(exbytx[queryHits(top1_hits)])
fixed(gg_chr9zoom) <- FALSE

gg_e1zoom <- autoplot(e1set_gr[subjectHits(top1_hits)], alpha=1/5)
fixed(gg_e1zoom) <- FALSE

tracks("UCSC KG" = gg_chr9zoom,
       "e1 genes" = gg_e1zoom,
       heights=c(1/3, 1/3), xlab="chr9 positions") +
    theme_tracks_sunset()



#' ------------------------------------------------------------------------
#' ## two exon genes
#'

#' We can look at the number of exons in each connected component (same
#' as shown above)
head(table(exon_cnts))

#' We can also look at the number of exons and junctions in each connected
#' component. Note, ideally, no 2 e/j objects should exist.
ej_cnts <- table(chr9$gIdx)
head(table(ej_cnts))

#' looking at some genes which have 2 "exons" but no splicing
names(ej_cnts[ej_cnts == 2])[1:10]
chr9[chr9$gIdx == names(ej_cnts[ej_cnts == 2])[1], 1:6]



#' ------------------------------------------------------------------------
#' ## genes with greater than one exon
#'
        
#' also construct subset with more than a single exon
e2pset <- subset(chr9_e, !(gIdx %in% e1_gIdx))
e2pset$chr <- "chr9"

#' as with single exon genes, we construct a `GRangesList` object to use with `ggbio`
e2pset_gr <- makeGRangesFromDataFrame(e2pset,
                                      seqnames.field="chr",
                                      keep.extra.columns=TRUE)
seqlengths(e2pset_gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)["chr9"]


#' compute distribution of average expression over each exon and plot distributions
#' grouped by the number of exons each gene
#+ hide=TRUE, echo=FALSE, fig.height=5, fig.width=12
e_mean <- rowMeans(chr9_e[paste0("s", 1:177)])
e_mean <- data.frame(exp=e_mean,
                     n_exon=exon_cnts[chr9_e$gIdx]) 
e_median <- apply(chr9_e[paste0("s", 1:177)], 1, median)
e_median <- data.frame(exp=e_median,
                       n_exon=exon_cnts[chr9_e$gIdx]) 

ggplot(e_mean) + 
    geom_boxplot(aes(x=n_exon, y=exp, group=n_exon),
                 color=brew_set1[2], fill=brew_set1[3]) +
    theme_bw() +
    scale_y_continuous(trans=scales::log1p_trans(),
                       breaks=c(0, 10^(0:5))) +
    xlab("# exons") +
    ylab("log10 mean expr across 177 samples") +
    ggtitle("dist of exon level expr for chr9 components, grouped by number of exons")

ggplot(e_median) + 
    geom_boxplot(aes(x=n_exon, y=exp, group=n_exon),
                 color=brew_set1[2], fill=brew_set1[3]) +
    theme_bw() +
    scale_y_continuous(trans=scales::log1p_trans(),
                       breaks=c(0, 10^(0:5))) +
    xlab("# exons") +
    ylab("log10 median expr across 177 samples") +
    ggtitle("dist of expr for chr9 components, grouped by number of exons")



#' ------------------------------------------------------------------------
#' ## p16/p14 region
#'

#' ucsc browser for chr9:21,965,000-21,995,000
#' ![ucsc](images/ucsc_p16.png)  
#'

#' look in region around cdkn2a (p16/p14)  (on rev strand)
#'
p16_eid <- select(org.Hs.eg.db, keys="CDKN2A",
                  columns="ENTREZID", keytype="SYMBOL")
p16_txid <- select(txdb, keys=p16_eid$ENTREZID,
                   columns=c("TXID", "TXNAME"), keytype="GENEID")
p16_tx <- exbytx[as.character(p16_txid$TXID)]
p16_uniontx <- reduce(unlist(p16_tx))
bounds <- c(start(range(p16_uniontx)), end(range(p16_uniontx)))

cands1 <- (chr9_e$start > min(bounds)-10000 &
               chr9_e$gStart < max(bounds)) |
          (chr9_e$gStop < max(bounds)+10000 &
               chr9_e$gStop > min(bounds))
cand_g <- chr9_e$gIdx[cands1]

p16_set <- subset(chr9_e, gIdx %in% unique(cand_g))
p16_set$chr <- "chr9"

#' convert to `GRanges` object to use with `ggbio`
p16_gr <- makeGRangesFromDataFrame(p16_set,
                                   seqnames.field="chr",
                                   keep.extra.columns=TRUE)
seqlengths(p16_gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)["chr9"]

#' convert to `GRangesList` object to separate out groups of genes
p16_gl <- split(p16_gr, mcols(p16_gr)$gIdx)


#' plot
#+ fig.height=2, fig.width=10
gg_p16models <- autoplot(p16_tx)
fixed(gg_p16models) <- FALSE

gg_p16_cc <- autoplot(p16_gl)
fixed(gg_p16_cc) <- FALSE

tracks("UCSC" = gg_p16models,
       "concomp" = gg_p16_cc,
       heights=c(1/3, 1/3), xlab="chr9 positions") +
    theme_tracks_sunset()



#' ------------------------------------------------------------------------
#' ### Plotting expression and splicing at CDKN2A locus
#' 

#' From above, most likely connected component was `gene1791`
p16_ej <- subset(chr9, gIdx == "gene1791")
p16_ej$chr <- "chr9"

p16_gr1 <- makeGRangesFromDataFrame(p16_ej,
                                    seqnames.field="chr",
                                    keep.extra.columns=TRUE)
seqlengths(p16_gr1) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)["chr9"]
strand(p16_gr1) <- "-"

p16_gl1 <- split(p16_gr1, mcols(p16_gr1)$kind)


#' one approach to plotting expression with equal widths given to each
#' 'exon'.
#+ fig.height=6, fig.width=10

p16_gl1j <- p16_gl1$j
aa <- as.data.frame(findOverlaps(p16_gl1j, p16_gl1j, type="within"))

aa2 <- aggregate(aa$queryHits, list(aa$subjectHits), c)
aa2$len <- sapply(aa2$x, length)

aa2$h <- 0
for (i in sort(unique(aa2$len)))
    aa2$h[which(aa2$len == i)] <- sapply(aa2$x[which(aa2$len == i)],
              function(z) max(unlist(aa2$h[z]))+.3)

aa2$h2 <- aa2$h * (-1)^(1+aa2$h/.3)

mcols(p16_gl1j)$offset <- aa2$h
mcols(p16_gl1j)$offset2 <- 0.3*aa2$len
mcols(p16_gl1j)$offset3 <- aa2$h2

bb <-
    ggplot(p16_gl1$e) +
    geom_rect(fill="grey", color="grey30", size=.3) +
    geom_chevron(p16_gl1j, color="grey30", offset="offset3",
                 stat="identity", size=.5, alpha=1/2, aes(y=I(1.4))) +
    theme_alignment() +
    scale_x_continuous(breaks=1e4*
                           seq(floor(start(range(p16_gl1$e))/1e4),
                               ceiling(end(range(p16_gl1$e))/1e4)))

p16_te <- as.data.frame(mcols(p16_gl1$e)[paste0("s", 1:177)])
p16_te$eid <- 1:nrow(p16_te)
p16_te <- reshape2::melt(p16_te, id.var="eid")
p16_te$value <- log10(1+p16_te$value)

ggplot(p16_te, aes(xmin=eid-.4, xmax=eid+.4,
                   ymin=value, ymax=value+.05, group=variable)) +
    geom_rect(color="grey60", fill="grey30", alpha=1/5) + guides(color=FALSE) +
    geom_line(aes(x=eid, y=value), alpha=1/5) +
    theme_bw() +
    xlab("exon") + ylab("log10 expr")

#' adding annotations
#+ fig.height=6, fig.width=10
plotRangesLinkedToData(p16_gl1$e, stat.y=paste0("s", 1:177),
                       linetype=0, annotation=bb)


#' Alternative approach to plotting with no re-scaling.
#' Using genomic coordinates.
#+ fig.height=5, fig.width=10
p16_exp <- reshape2::melt(subset(p16_ej, kind=="e"),
                          id.vars=c("gIdx", "gStart", "gStop", "kind",
                              "start", "stop", "chr"))
p16_exp$log1p <- log10(p16_exp$value+1)

gg_logexp <-
    ggplot(p16_exp,
           aes(xmin=start, xmax=stop,
               ymin=log1p, ymax=log1p+.05)) +
    geom_rect(fill="white", color="grey30", alpha=1/5) +
    geom_line(aes(x=(start+stop)/2, y=log1p, group=variable), alpha=1/10) +
    theme_bw()

tracks(gg_logexp,
       autoplot(GRangesList(p16_gl1$e), gap.geom="arrow",
                fill="grey80", color="grey30") + theme_bw(),
       heights=c(3, 1))




#' looking at clustering along the exons at p16
#'
matplot(log10(1+p16_ej[p16_ej$kind == "e", -c(1:6, 184)]), type="l", lty=1)


#' PCA decomposition
pc <- prcomp(log10(1+t(p16_ej[p16_ej$kind == "e", -c(1:6, 184)])))

#' scores
plot(pc$x[, 1:2], pch=16, col=brew_set1[2])

#' loadings
#+ fig.height=3, fig.weight=10
qplot(data=reshape2::melt(pc$rotation[, 1:2]), x=rep(1:22, 2), y=value, color=Var2,
      geom="line") +
    theme_bw() +
    scale_x_continuous(breaks=1:22)

#' look at clusters
groups <- as.numeric(pc$x[, 1]>-1)*(1+as.numeric(pc$x[, 2]>0)) + 1

#' scores by clusters
plot(pc$x[, 1:2], pch=16, col=brew_set1[groups])


#' looking at clusters along genomic coordinates
#'
p16_exp$groups <- groups[as.numeric(substr(p16_exp$variable, 2, 5))]

ggplot(p16_exp,
       aes(xmin=start, xmax=stop,
           ymin=log1p, ymax=log1p+.05,
           fill=factor(groups))) +
    geom_rect(alpha=1/5) +
    geom_line(aes(x=(start+stop)/2, y=log1p,
                  group=variable, color=factor(groups)),
              alpha=1/5) +
    theme_bw()


#' looking at clusters with all exons having equal width,
#' onlying looking at first 11 exons
#+ fig.height=3, fig.width=10
p16_te$groups <- groups[as.numeric(substr(p16_te$variable, 2, 5))]

ggplot(subset(p16_te, eid <= 11),
       aes(xmin=eid-.4, xmax=eid+.4,
           ymin=value, ymax=value+.05,
           group=variable,
           fill=factor(groups),
           color=factor(groups))) +
    geom_rect(color="grey90", alpha=1/10) + guides(color=FALSE) +
    geom_line(aes(x=eid, y=value), alpha=1/2) +
    scale_x_continuous("exon", breaks=1:11) +
    theme_bw() +
    ylab("log10 expr")



#' ------------------------------------------------------------------------
#' ### Process p16 locus for graph-based analysis
#'

p16_v <- p16_ej[p16_ej$kind == "e", c("start", "stop")]
p16_e <- p16_ej[p16_ej$kind == "j", c("start", "stop")]

## determine explicit edgeset
p16_e$a <- sapply(p16_e$start, function(x) which(p16_v$stop == x))
p16_e$b <- sapply(p16_e$stop, function(x) which(p16_v$start == x))
p16_el1 <- as.matrix(p16_e[c("a", "b")])

## determine 'edges' that exist between consecutive exons
conseq <- which((p16_v$start[-1] - p16_v$stop[-22]) == 1)
p16_el2 <- cbind("a"=conseq, "b"=conseq+1)



#' ------------------------------------------------------------------------
#' ### Plot splicing graph for Sample 2
#'

p16_vw2 <- p16_ej[p16_ej$kind == "e", "s2"]
p16_ew2 <- p16_ej[p16_ej$kind == "j", "s2"]

p16_adj2 <- matrix(0, nrow(p16_v), nrow(p16_v))
diag(p16_adj2) <- p16_vw2
p16_adj2[p16_el1] <- p16_ew2
p16_adj2[p16_el2] <- apply(matrix(diag(p16_adj2)[p16_el2], ncol=2), 1, min)

g2 <- graph.adjacency(p16_adj2, weighted=TRUE, diag=FALSE)

edge_h <- 1/apply(get.edgelist(g2), 1, diff)
edge_h[edge_h == 1] <- 1e-10


#' example of an adjacency matrix
get.adjacency(graph.adjacency(p16_adj2))


#' nicer first draft of a splicing graph
#+ fig.width=12, fig.height=8
plot(g2, layout=cbind(1:22, 0),
     edge.curved=7*edge_h,
     edge.arrow.size=.3,
     edge.width=log2(E(g2)$weight+1),
     vertex.shape="rectangle", vertex.label=NA,
     vertex.size=4.5, vertex.size2=2)



#' ------------------------------------------------------------------------
#' ### Analyze 3-dimensional tensor of adjacency matrices
#'

#' #### construct adjacency tensor
adj_array <- array(0, dim=c(177, nrow(p16_v), nrow(p16_v)))
for (i in 1:177) {
    diag(adj_array[i, , ]) <- p16_ej[p16_ej$kind == "e", paste0("s", i)]
    adj_array[cbind(i, p16_el1)] <- p16_ej[p16_ej$kind == "j", paste0("s", i)]
    adj_array[cbind(i, p16_el2)] <-
        apply(matrix(diag(adj_array[i, , ])[p16_el2], ncol=2), 1, min)
}


#' Load SigFuge labels from previous data as a sanity check
p16_labs_sf <- read.table("p16_labs_sf.txt", header=TRUE)$clusters


#' #### Tucker PCA decomposition
tucker_pca <- PCAn(log10(1+adj_array), dim=c(3, 3, 3))
tucker_scores <- as.data.frame(t(tucker_pca[[1]]$v[1:2, ]))
tucker_scores$labs <- 0 + (tucker_scores$V1 > 0.05)*(1+(tucker_scores$V2 > 0))

#' Tucker decomposition (scores)
#+ fig.height=8, fig.weight=8
qplot(data=tucker_scores,
      x=V1, y=V2, color=factor(labs)) +
    theme_bw() +
        guides(color=FALSE)

##compare with SigFuge labels
table(tucker_scores$labs, p16_labs_sf)


#' #### PARAFAC decomposition
parafac_pca <- CANDPARA(log10(1+adj_array), dim=3)
parafac_scores <- as.data.frame(t(parafac_pca[[1]]$v[1:2, ] * parafac_pca[[3]]$d[1:2]))
parafac_scores$labs <- 0 + (parafac_scores$V1 > 3)*(1+(parafac_scores$V2 > 0))

#' PARAFAC decomposition (scores)
#+ fig.height=8, fig.weight=8
qplot(data=parafac_scores,
      x=V1, y=V2, color=factor(labs)) +
    theme_bw() +
    guides(color=FALSE)

##compare with SigFuge labels
table(parafac_scores$labs, p16_labs_sf)


#' #### PARAFAC decomposition without splices
adj_array_diag <- array(0, dim=c(177, nrow(p16_v), nrow(p16_v)))
for (i in 1:177) {
    diag(adj_array_diag[i, , ]) <- diag(adj_array[i, , ])
}
parafac_pca_diag <- CANDPARA(log10(1+adj_array_diag), dim=3)
parafac_scores_diag <- as.data.frame(t(parafac_pca_diag[[1]]$v[1:2, ] *
                                           parafac_pca_diag[[3]]$d[1:2]))
parafac_scores_diag$labs <- 0 + (parafac_scores_diag$V1 > 1.5)*(1+(parafac_scores_diag$V2 > 0))

#' PARAFAC decomposition with only exons (scores)
#+ fig.height=8, fig.weight=8
qplot(data=parafac_scores_diag,
      x=V1, y=V2, color=factor(labs)) +
    theme_bw() +
    guides(color=FALSE)

##compare with SigFuge labels
table(parafac_scores_diag$labs, p16_labs_sf)


#' #### PARAFAC decomposition without exon expr
adj_array_offdiag <- adj_array
for (i in 1:177) {
    diag(adj_array_offdiag[i, , ]) <- 0
}
parafac_pca_offdiag <- CANDPARA(log10(1+adj_array_offdiag), dim=3)
parafac_scores_offdiag <- as.data.frame(t(parafac_pca_offdiag[[1]]$v[1:2, ] *
                                           parafac_pca_offdiag[[3]]$d[1:2]))
parafac_scores_offdiag$labs <- 0 + (parafac_scores_offdiag$V1 > 1.5)*(1+(parafac_scores_offdiag$V2 > 0))

#' PARAFAC decomposition with only splices (scores)
#+ fig.height=8, fig.weight=8
qplot(data=parafac_scores_offdiag,
      x=V1, y=V2, color=factor(labs)) +
    theme_bw() +
    guides(color=FALSE)

##compare with SigFuge labels
table(parafac_scores_offdiag$labs, p16_labs_sf)






#' ------------------------------------------------------------------------
#' ### Last Updated
#'

Sys.time()



#' ------------------------------------------------------------------------
#' ### Session Information
#'

sessionInfo()
