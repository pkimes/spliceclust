#' ---
#' title: "LUSC chr9 (CDKN2A) Analysis"
#' output:
#'   html_document:
#'     toc: true
#'     fig_width: 12
#'     fig_height: 6
#'     fig_caption: false
#'   md_document:
#'     toc: true
#'     fig_width: 12
#'     fig_height: 6
#' ---

#+ echo=FALSE, message=FALSE, warning=FALSE, hide=TRUE
knitr::opts_chunk$set(collapse=TRUE, comment="##",
                      echo=FALSE, warning=FALSE,
                      error=FALSE, message=FALSE,
                      hide=TRUE)
start_time <- proc.time()


## load necessary packages
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


##define nicer colors
brew_set1 <- brewer.pal(9, "Set1")
brew_pastel2 <- brewer.pal(8, "Pastel2")
brew_accent <- brewer.pal(8, "Accent")
brew_greys <- brewer.pal(9, "Greys")


##specify project directory and source functions
root <-"/Users/pkimes/Dropbox/Git/spliceclust/"
source(paste0(root, "R/readchr.R"))


## load chromesome 9 data for 177 lusc samples
## complete datasest
chr9 <- readchr(paste0(root, "data/lusc/chr9_gene.txt"), 177)
## exon only dataset
chr9_e <- subset(chr9, kind == "e")





#'
#' <br> <br> <br> <hr> <br>
#' 
#' ## Chr9 Summary
#'

#'
#' <br> <br>
#' **# of exons in each CC**
#+ hide=FALSE
exon_cnts <- table(chr9_e$gIdx)
summary(as.numeric(exon_cnts))


#'
#' <br> <br>
#' **# of CC total: `r length(exon_cnts)`**
#' <br>
#' **# of CC w/ one exon: `r sum(exon_cnts==1)`**
#'


#'
#' <br> <br>
#' **# of CCs by # of exons**
#+ hide=FALSE
qplot(x=Freq, data=data.frame(exon_cnts), binwidth=I(1),
      color=I(brew_greys[5]), fill=I(brew_pastel2[1])) +
    theme_bw() +
    scale_y_continuous(trans=scales::log1p_trans(),
                       breaks=c(0, 10^(0:5))) +
    xlab("number of exons") +
    ylab("log(1 + #models) , real counts shown") +
    ggtitle("distribution of exon counts for connected components")


#'
#' <br> <br>
#' **mean exon expr. by # of exons**
#+ hide=FALSE, fig.height=5, fig.width=12
e_mean <- rowMeans(chr9_e[paste0("s", 1:177)])
e_mean <- data.frame(exp=e_mean,
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


#'
#' <br> <br>
#' **median exon expr. by # of exons**
#+ hide=FALSE, fig.height=5, fig.width=12
e_median <- apply(chr9_e[paste0("s", 1:177)], 1, median)
e_median <- data.frame(exp=e_median,
                       n_exon=exon_cnts[chr9_e$gIdx]) 
ggplot(e_median) + 
    geom_boxplot(aes(x=n_exon, y=exp, group=n_exon),
                 color=brew_set1[2], fill=brew_set1[3]) +
    theme_bw() +
    scale_y_continuous(trans=scales::log1p_trans(),
                       breaks=c(0, 10^(0:5))) +
    xlab("# exons") +
    ylab("log10 median expr across 177 samples") +
    ggtitle("dist of expr for chr9 components, grouped by number of exons")





#'
#' <br> <br> <br> <hr> <br>
#' 
#' ## Single exon CCs
#'

##subset on only single exon CCs
e1_gIdx <- names(which(exon_cnts == 1))
e1set <- subset(chr9_e, gIdx %in% e1_gIdx)
e1set$chr <- "chr9"

##convert to `GRanges` object to use with `ggbio`
e1set_gr <- makeGRangesFromDataFrame(e1set,
                                     seqnames.field="chr",
                                     keep.extra.columns=TRUE)
seqlengths(e1set_gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)["chr9"]

##load hg19 cytoband information
data(ideoCyto, package = "biovizBase")
biovizBase::isIdeogram(ideoCyto$hg19)

##construct ideogram plot for comparison
gg_chr9 <- autoplot(subset(ideoCyto$hg19, seqnames(ideoCyto$hg19) == "chr9"),
                    layout="karyogram", cytoband=TRUE) +
                        guides(fill=FALSE)
fixed(gg_chr9) <- TRUE

##construct karyogram with chr9
gg_e1 <- autoplot(e1set_gr, layout="karyogram", alpha=1/100)
fixed(gg_e1) <- TRUE


#'
#' <br> <br>
#' **location on chr9**
#+ hide=FALSE, fig.height=2, fig.width=10
tracks("chr" = gg_chr9, "e1 genes" = gg_e1,
       heights=c(.5, .5), xlab="chr9 positions") +
    theme_tracks_sunset()


##compute statistics
e1_exp <- as.matrix(mcols(e1set_gr)[paste0("s", 1:177)])
e1_expstat <- data.frame("med"=apply(e1_exp, 1, median),
                         "mean"=rowMeans(e1_exp),
                         "sd"=apply(e1_exp, 1, sd),
                         "mad"=apply(e1_exp, 1, mad))


#'
#' <br> <br>
#' **median vs. MAD expr. of single exon CCs**
#+ hide=FALSE, fig.height=8, fig.width=8
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


##compute distribution of exon length
e1_widths <- data.frame("width"=width(ranges(e1set_gr)))
qplot(x=width, data=e1_widths,
      color=I(brew_greys[5]), fill=I(brew_pastel2[1])) +
    theme_bw() +
    scale_y_continuous(trans=scales::log1p_trans(),
                       breaks=c(0, 10^(0:5))) +
    xlab("length of exon (bp)") +
    ylab("log(1 + #models) , real counts shown") +
    ggtitle("distribution of exon lengths for components with single exon")


##load transcript annotations and only look at chr9
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
isActiveSeq(txdb)[seqlevels(txdb)] <- FALSE
isActiveSeq(txdb)["chr9"] <- TRUE

##want to work at exon level aggregated by genes or transcripts
exbygene <- exonsBy(txdb, "gene")
exbytx <- exonsBy(txdb, "tx")


##calculate overlap with genes and transcripts
##overlap with genes
hits_g <- findOverlaps(exbygene, e1set_gr)
ucsc_ovlp_g <- select(org.Hs.eg.db, keys=names(exbygene[queryHits(hits_g)]),
                      columns="SYMBOL", keytype="ENTREZID")

##overlap with transcripts
hits_tx <- findOverlaps(exbytx, e1set_gr)
ucsc_ovlp_tx <- select(txdb, keys=names(exbytx[queryHits(hits_tx)]),
                       columns=c("GENEID", "TXNAME"), keytype="TXID")
ucsc_ovlp_tx$idx <- 1:nrow(ucsc_ovlp_tx)
temp <- select(org.Hs.eg.db, keys=unique(ucsc_ovlp_tx$GENEID),
               columns="SYMBOL", keytype="ENTREZID")
names(temp) <- c("GENEID", "SYMBOL")
ucsc_ovlp_tx <- merge(ucsc_ovlp_tx, temp, all.x=TRUE, sort=FALSE)
ucsc_ovlp_tx <- ucsc_ovlp_tx[order(ucsc_ovlp_tx$idx), ]


#'
#' <br> <br>
#' There were `r length(unique(queryHits(hits_g)))` unique UCSC genes that overlapped with
#' `r length(unique(subjectHits(hits_g)))` unique CCs. Total number of
#' overlaps was `r length(queryHits(hits_g))`.
#' 
#+ hide=FALSE, fig.width=10, fig.height=6
ovlp_table_g <- data.frame(table(ucsc_ovlp_g$SYMBOL))
ggplot(data=ovlp_table_g) +
    geom_text(aes(label=Var1, x=1:length(Var1), y=Freq)) +
    theme_bw() +
    xlab("") +
    ylab("# overlapping components") +
    ggtitle("UCSC KnownGenes overlapping single gene connected components")


##similar plot as above, but at the scale of transcripts.
ovlp_table_tx <- data.frame(table(ucsc_ovlp_tx$SYMBOL))
ggplot(data=ovlp_table_tx) +
    geom_text(aes(label=Var1, x=1:length(Var1), y=Freq)) +
    theme_bw() +
    xlab("") +
    ylab("# overlapping components") +
    ggtitle("UCSC KnownGene transcripts overlapping single gene connected components")


#'
#' <br> <br>
#' We now add UCSC KnownGenes to above `tracks` plots. Note that these are not
#' real transcripts, just the 'union' transcripts for each gene constructed by
#' joining all reported isoforms. (note: need to specify `fixed() <- TRUE` for
#' ucsc track since plot includes ideogram.)
#+ hide=FALSE, fig.height=3, fig.width=10
gg_ucsc <- autoplot(exbygene)
fixed(gg_ucsc) <- TRUE

tracks("chr" = gg_chr9,
       "e1 genes" = gg_e1,
       "UCSC KG" = gg_ucsc,
       heights=c(1/3, 1/3, 1/3), xlab="chr9 positions") +
    theme_tracks_sunset()


#'
#' <br> <br>
#' Next, we focus on `r ovlp_table_g$Var1[which.max(ovlp_table_g$Freq)]`  which
#' overlapped with `r max(ovlp_table_g$Freq)` single exon CCs.
#+ hide=FALSE, fig.height=2, fig.width=10
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





#'
#' <br> <br> <br> <hr> <br> 
#'
#' ## Two exon CCs
#'

#'
#' <br> <br>
#' **# of exons in each CC**
#+ hide=FALSE
head(table(exon_cnts))


#'
#' <br> <br>
#' **# of exons/junctions in each CC**
#+ hide=FALSE
ej_cnts <- table(chr9$gIdx)
head(table(ej_cnts))


#'
#' <br> <br>
#' **some CCs with 2 exons but no junctions**
#+ hide=FALSE
names(ej_cnts[ej_cnts == 2])[1:10]
chr9[chr9$gIdx == names(ej_cnts[ej_cnts == 2])[1], 1:6]





#'
#' <br> <br> <br> <hr> <br> 
#'
#' ## p16/p14 region
#'

#'
#' <br> <br>
#' **ucsc browser for chr9:21,960,000-22,090,000**
#' ![ucsc](images/ucsc_p16.png)  
#'


#'
#' <br> <br> <br> <hr> <br>
#' 
#' ### CC gene models
#'

#look in region around cdkn2a (p16/p14)  (on rev strand)
p16_eid <- select(org.Hs.eg.db, keys="CDKN2A",
                  columns="ENTREZID", keytype="SYMBOL")
p16_txid <- select(txdb, keys=p16_eid$ENTREZID,
                   columns=c("TXID", "TXNAME"), keytype="GENEID")
p16_tx <- exbytx[as.character(p16_txid$TXID)]
p16_uniontx <- reduce(unlist(p16_tx))
bounds <- c(start(range(p16_uniontx)), end(range(p16_uniontx)))

##take all CCs in region 10k up/down of bounds of gene model
cands1 <- (chr9_e$start > min(bounds)-10000) &
          (chr9_e$stop < max(bounds)+10000)
cand_g <- chr9_e$gIdx[cands1]

##find all CC gene indices in region
p16_set <- subset(chr9_e, gIdx %in% unique(cand_g))
p16_set$chr <- "chr9"

##convert to `GRanges` object to use with `ggbio`
p16_gr <- makeGRangesFromDataFrame(p16_set,
                                   seqnames.field="chr",
                                   keep.extra.columns=TRUE)
seqlengths(p16_gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)["chr9"]

##convert to `GRangesList` object to separate out groups of genes
p16_gl <- split(p16_gr, mcols(p16_gr)$gIdx)

##plot all CDKN2A UCSC gene models
gg_p16models <- autoplot(p16_tx)
fixed(gg_p16models) <- FALSE

##plot all candidate CCs
gg_p16_cc <- autoplot(p16_gl)
fixed(gg_p16_cc) <- FALSE

#'
#' <br> <br>
#' **compare UCSC CDKN2A models vs. CCs in region**
#+ hide=FALSE, fig.height=2, fig.width=10
tracks("UCSC" = gg_p16models,
       "concomp" = gg_p16_cc,
       heights=c(1/3, 1/3), xlab="chr9 positions") +
    theme_tracks_sunset()


#'
#' <br> <br>
#' **include nearby genes**
#+ hide=FALSE, fig.height=2, fig.width=10
nearby_eid <- select(org.Hs.eg.db,
                      keys=c("CDKN2A", "CDKN2B", "CDKN2B-AS1"),
                  columns="ENTREZID", keytype="SYMBOL")
nearby_txid <- select(txdb, keys=nearby_eid$ENTREZID,
                   columns=c("TXID", "TXNAME"), keytype="GENEID")
nearby_tx <- exbytx[as.character(nearby_txid$TXID)]

gg_nearby_models <- autoplot(nearby_tx)
fixed(gg_nearby_models) <- FALSE
tracks("UCSC CDKN2A" = gg_p16models,
       "UCSC nearby" = gg_nearby_models,
       "concomp" = gg_p16_cc,
       heights=c(1/4, 2/4, 1/4), xlab="chr9 positions") +
    theme_tracks_sunset()


#'
#' <br> <br>
#' **most likely connected component was `gene1791`**
p16_ej <- subset(chr9, gIdx == "gene1791")
p16_ej$chr <- "chr9"

p16_gr1 <- makeGRangesFromDataFrame(p16_ej,
                                    seqnames.field="chr",
                                    keep.extra.columns=TRUE)
seqlengths(p16_gr1) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)["chr9"]
strand(p16_gr1) <- "-"

p16_gl1 <- split(p16_gr1, mcols(p16_gr1)$kind)





#'
#' <br> <br> <br> <hr> <br>
#' 
#' ### Plotting exon expression (one CC)
#' 

#'
#' <br> <br>
#'  **equal width exons**
#+ hide=FALSE, fig.height=6, fig.width=10
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



#'
#' <br> <br>
#' **adding annotations**
#+ hide=FALSE, fig.height=6, fig.width=10
plotRangesLinkedToData(p16_gl1$e, stat.y=paste0("s", 1:177),
                       linetype=0, annotation=bb)





#'
#' <br> <br> <br> <hr> <br>
#' 
#' ### PCA exon expression (one CC)
#' 

#'
#' <br> <br>
#' looking at clustering along the exons at p16
#+ hide=FALSE
matplot(log10(1+p16_ej[p16_ej$kind == "e", -c(1:6, 184)]), type="l", lty=1)


##PCA decomposition (exons separately)
pc <- prcomp(log10(1+t(p16_ej[p16_ej$kind == "e", -c(1:6, 184)])))

#'
#' <br> <br>
#' exon PCA decomposition (scores)
#+ hide=FALSE, fig.width=8, fig.height=8
plot(pc$x[, 1:2], pch=16, col=brew_set1[2])


#'
#' <br> <br>
#' exon PCA decomposition (loadings)
#+ hide=FALSE, fig.height=3, fig.weight=10
qplot(data=reshape2::melt(pc$rotation[, 1:2]), x=rep(1:22, 2), y=value, color=Var2,
      geom="line") +
    theme_bw() +
    scale_x_continuous(breaks=1:22)


#'
#' <br> <br>
#' look at clusters
#+ hide=FALSE
groups <- as.numeric(pc$x[, 1]>-1)*(1+as.numeric(pc$x[, 2]>0)) + 1


#'
#' <br> <br>
#' scores by clusters
#+ hide=FALSE, fig.width=8, fig.height=8
plot(pc$x[, 1:2], pch=16, col=brew_set1[groups])



#'
#' <br> <br>
#' expression by clusters
#+ hide=FALSE, fig.height=3, fig.width=10
p16_te$cluster <- rep(groups, each=max(p16_te$eid))
ggplot(p16_te, aes(xmin=eid-.4, xmax=eid+.4,
                   ymin=value, ymax=value+.05, group=variable,
                   color=as.factor(cluster), fill=as.factor(cluster))) +
    geom_rect(color="grey90", alpha=1/10) + guides(color=FALSE) +
    geom_line(aes(x=eid, y=value), alpha=1/3) +
    guides(color=FALSE, fill=FALSE) +
    theme_bw() +
    xlab("exon") + ylab("log10 expr")



#'
#' <br> <br>
#' looking at clusters with all exons having equal width,
#' onlying looking at first 11 exons
#+ hide=FALSE, fig.height=3, fig.width=10
ggplot(subset(p16_te, eid <= 11),
       aes(xmin=eid-.4, xmax=eid+.4,
           ymin=value, ymax=value+.05,
           group=variable,
           fill=as.factor(cluster),
           color=as.factor(cluster))) +
    geom_rect(color="grey90", alpha=1/10) + guides(color=FALSE) +
    geom_line(aes(x=eid, y=value), alpha=1/3) +
    scale_x_continuous("exon", breaks=1:11) +
    theme_bw() +
    ylab("log10 expr")




#'
#' <br> <br> <br> <hr> <br>
#' 
#' ### Plotting junction coverage (one CC)
#' 

##PCA decomposition (junctions separately)
pc_j <- prcomp(log10(1+t(p16_ej[p16_ej$kind == "j", -c(1:6, 184)])))


#'
#' <br> <br>
#' look at clusters
#+ hide=FALSE
groups <- as.numeric(pc_j$x[, 1]>-1)*(1+as.numeric(pc_j$x[, 2]>0)) + 1


#'
#' <br> <br>
#' scores by clusters
#+ hide=FALSE, fig.width=8, fig.height=8
plot(pc_j$x[, 1:2], pch=16, col=brew_set1[groups])



#'
#' <br> <br>
#' expression by clusters
#+ hide=FALSE, fig.height=3, fig.width=10
p16_tj <- as.data.frame(mcols(p16_gl1$j)[paste0("s", 1:177)])
p16_tj$eid <- 1:nrow(p16_tj)
p16_tj <- reshape2::melt(p16_tj, id.var="eid")
p16_tj$value <- log10(1+p16_tj$value)
p16_tj$cluster <- rep(groups, each=max(p16_tj$eid))

ggplot(p16_tj, aes(xmin=eid-.4, xmax=eid+.4,
                   ymin=value, ymax=value+.05, group=variable,
                   color=as.factor(cluster), fill=as.factor(cluster))) +
    geom_rect(color=NA, alpha=1/5) + guides(color=FALSE) +
    guides(color=FALSE, fill=FALSE) +
    theme_bw() +
    xlab("exon") + ylab("log10 expr")

ggplot(p16_tj, aes(xmin=eid-.4, xmax=eid+.4,
                   ymin=value, ymax=value+.05, group=variable,
                   color=as.factor(cluster), fill=as.factor(cluster))) +
    geom_rect(color="grey90", alpha=1/10) + guides(color=FALSE) +
    geom_line(aes(x=eid, y=value), alpha=1/3) +
    guides(color=FALSE, fill=FALSE) +
    theme_bw() +
    xlab("exon") + ylab("log10 expr")







#'
#' <br> <br> <br> <hr> <br>
#' 
#' ### Graph-based view (one CC)
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




#'
#' <br> <br> <br> <hr> <br> 
#'
#' ## R Environment
#'

#' 
#' ### Running Time
#+ hide=FALSE
proc.time() - start_time


#'
#' ### Last Updated
#+ hide=FALSE
Sys.time()


#'
#' ### Session Information
#+ hide=FALSE
sessionInfo()

