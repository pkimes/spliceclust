#' ---
#' title: "LUSC chr19 (KLK12) Analysis"
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


## load chromesome 19 data for 177 lusc samples
## complete datasest
chr19 <- readchr(paste0(root, "data/lusc/chr19_gene.txt"), 177)
## exon only dataset
chr19_e <- subset(chr19, kind == "e")





#'
#' <br> <br> <br> <hr> <br>
#' 
#' ## Chr19 Summary
#'

#'
#' <br> <br>
#' **# of exons in each CC**
#+ hide=FALSE
exon_cnts <- table(chr19_e$gIdx)
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
e_mean <- rowMeans(chr19_e[paste0("s", 1:177)])
e_mean <- data.frame(exp=e_mean,
                     n_exon=exon_cnts[chr19_e$gIdx]) 

ggplot(e_mean) + 
    geom_boxplot(aes(x=n_exon, y=exp, group=n_exon),
                 color=brew_set1[2], fill=brew_set1[3]) +
    theme_bw() +
    scale_y_continuous(trans=scales::log1p_trans(),
                       breaks=c(0, 10^(0:5))) +
    xlab("# exons") +
    ylab("log10 mean expr across 177 samples") +
    ggtitle("dist of exon level expr for chr19 components, grouped by number of exons")


#'
#' <br> <br>
#' **median exon expr. by # of exons**
#+ hide=FALSE, fig.height=5, fig.width=12
e_median <- apply(chr19_e[paste0("s", 1:177)], 1, median)
e_median <- data.frame(exp=e_median,
                       n_exon=exon_cnts[chr19_e$gIdx]) 
ggplot(e_median) + 
    geom_boxplot(aes(x=n_exon, y=exp, group=n_exon),
                 color=brew_set1[2], fill=brew_set1[3]) +
    theme_bw() +
    scale_y_continuous(trans=scales::log1p_trans(),
                       breaks=c(0, 10^(0:5))) +
    xlab("# exons") +
    ylab("log10 median expr across 177 samples") +
    ggtitle("dist of expr for chr19 components, grouped by number of exons")





#'
#' <br> <br> <br> <hr> <br>
#' 
#' ## Single exon CCs
#'

##subset on only single exon CCs
e1_gIdx <- names(which(exon_cnts == 1))
e1set <- subset(chr19_e, gIdx %in% e1_gIdx)
e1set$chr <- "chr19"

##convert to `GRanges` object to use with `ggbio`
e1set_gr <- makeGRangesFromDataFrame(e1set,
                                     seqnames.field="chr",
                                     keep.extra.columns=TRUE)
seqlengths(e1set_gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)["chr19"]

##load hg19 cytoband information
data(ideoCyto, package = "biovizBase")
biovizBase::isIdeogram(ideoCyto$hg19)

##construct ideogram plot for comparison
gg_chr19 <- autoplot(subset(ideoCyto$hg19, seqnames(ideoCyto$hg19) == "chr19"),
                    layout="karyogram", cytoband=TRUE) +
                        guides(fill=FALSE)
fixed(gg_chr19) <- TRUE

##construct karyogram with chr19
gg_e1 <- autoplot(e1set_gr, layout="karyogram", alpha=1/100)
fixed(gg_e1) <- TRUE


#'
#' <br> <br>
#' **location on chr19**
#+ hide=FALSE, fig.height=2, fig.width=10
tracks("chr" = gg_chr19, "e1 genes" = gg_e1,
       heights=c(.5, .5), xlab="chr19 positions") +
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


##load transcript annotations and only look at chr19
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
isActiveSeq(txdb)[seqlevels(txdb)] <- FALSE
isActiveSeq(txdb)["chr19"] <- TRUE

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

tracks("chr" = gg_chr19,
       "e1 genes" = gg_e1,
       "UCSC KG" = gg_ucsc,
       heights=c(1/3, 1/3, 1/3), xlab="chr19 positions") +
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

gg_chr19zoom <- autoplot(exbytx[queryHits(top1_hits)])
fixed(gg_chr19zoom) <- FALSE

gg_e1zoom <- autoplot(e1set_gr[subjectHits(top1_hits)], alpha=1/5)
fixed(gg_e1zoom) <- FALSE

tracks("UCSC KG" = gg_chr19zoom,
       "e1 genes" = gg_e1zoom,
       heights=c(1/3, 1/3), xlab="chr19 positions") +
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
ej_cnts <- table(chr19$gIdx)
head(table(ej_cnts))


#'
#' <br> <br>
#' **some CCs with 2 exons but no junctions**
#+ hide=FALSE
names(ej_cnts[ej_cnts == 2])[1:10]
chr19[chr19$gIdx == names(ej_cnts[ej_cnts == 2])[1], 1:6]





#'
#' <br> <br> <br> <hr> <br> 
#'
#' ## KLK12 region
#'


#'
#' <br> <br> <br> <hr> <br>
#' 
#' ### CC gene models
#'

#look in region around KLK12  (on rev strand)
klk12_eid <- select(org.Hs.eg.db, keys="KLK12",
                  columns="ENTREZID", keytype="SYMBOL")
klk12_txid <- select(txdb, keys=klk12_eid$ENTREZID,
                   columns=c("TXID", "TXNAME"), keytype="GENEID")
klk12_tx <- exbytx[as.character(klk12_txid$TXID)]
klk12_uniontx <- reduce(unlist(klk12_tx))
bounds <- c(start(range(klk12_uniontx)), end(range(klk12_uniontx)))

##take all CCs in region 10k up/down of bounds of gene model
cands1 <- (chr19_e$start > min(bounds)-10000) &
          (chr19_e$stop < max(bounds)+10000)
cand_g <- chr19_e$gIdx[cands1]

##find all CC gene indices in region
klk12_set <- subset(chr19_e, gIdx %in% unique(cand_g))
klk12_set$chr <- "chr19"

##convert to `GRanges` object to use with `ggbio`
klk12_gr <- makeGRangesFromDataFrame(klk12_set,
                                   seqnames.field="chr",
                                   keep.extra.columns=TRUE)
seqlengths(klk12_gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)["chr19"]

##convert to `GRangesList` object to separate out groups of genes
klk12_gl <- split(klk12_gr, mcols(klk12_gr)$gIdx)

##plot all KLK12 UCSC gene models
gg_klk12models <- autoplot(klk12_tx)
fixed(gg_klk12models) <- FALSE

##plot all candidate CCs
gg_klk12_cc <- autoplot(klk12_gl)
fixed(gg_klk12_cc) <- FALSE

#'
#' <br> <br>
#' **compare UCSC KLK12 models vs. CCs in region**
#+ hide=FALSE, fig.height=2, fig.width=10
tracks("UCSC" = gg_klk12models,
       "concomp" = gg_klk12_cc,
       heights=c(1/3, 1/3), xlab="chr19 positions") +
    theme_tracks_sunset()


#'
#' <br> <br>
#' **include nearby genes (KLK9, KLK10, KLK11, KLK12, KLK13)**
#+ hide=FALSE, fig.height=2, fig.width=10
nearby_eid <- select(org.Hs.eg.db,
                      keys=paste0("KLK", 9:13),
                  columns="ENTREZID", keytype="SYMBOL")
nearby_txid <- select(txdb, keys=nearby_eid$ENTREZID,
                   columns=c("TXID", "TXNAME"), keytype="GENEID")
nearby_tx <- exbytx[as.character(nearby_txid$TXID)]

gg_nearby_models <- autoplot(nearby_tx)
fixed(gg_nearby_models) <- FALSE
tracks("UCSC KLK12" = gg_klk12models,
       "UCSC nearby" = gg_nearby_models,
       "concomp" = gg_klk12_cc,
       heights=c(1/4, 2/4, 1/4), xlab="chr19 positions") +
    theme_tracks_sunset()


#'
#' <br> <br>
#' **most likely connected component was `gene9317`**
klk12_ej <- subset(chr19, gIdx == "gene9317")
klk12_ej$chr <- "chr19"

klk12_gr1 <- makeGRangesFromDataFrame(klk12_ej,
                                    seqnames.field="chr",
                                    keep.extra.columns=TRUE)
seqlengths(klk12_gr1) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)["chr19"]
strand(klk12_gr1) <- "-"

klk12_gl1 <- split(klk12_gr1, mcols(klk12_gr1)$kind)





#'
#' <br> <br> <br> <hr> <br>
#' 
#' ### Plotting exon expression (one CC)
#' 

#'
#' <br> <br>
#'  **equal width exons**
#+ hide=FALSE, fig.height=6, fig.width=10
klk12_gl1j <- klk12_gl1$j
aa <- as.data.frame(findOverlaps(klk12_gl1j, klk12_gl1j, type="within"))

aa2 <- aggregate(aa$queryHits, list(aa$subjectHits), c)
aa2$len <- sapply(aa2$x, length)

aa2$h <- 0
for (i in sort(unique(aa2$len)))
    aa2$h[which(aa2$len == i)] <- sapply(aa2$x[which(aa2$len == i)],
              function(z) max(unlist(aa2$h[z]))+.3)

aa2$h2 <- aa2$h * (-1)^(1+aa2$h/.3)

mcols(klk12_gl1j)$offset <- aa2$h
mcols(klk12_gl1j)$offset2 <- 0.3*aa2$len
mcols(klk12_gl1j)$offset3 <- aa2$h2

bb <-
    ggplot(klk12_gl1$e) +
    geom_rect(fill="grey", color="grey30", size=.3) +
    geom_chevron(klk12_gl1j, color="grey30", offset="offset3",
                 stat="identity", size=.5, alpha=1/2, aes(y=I(1.4))) +
    theme_alignment() +
    scale_x_continuous(breaks=1e4*
                           seq(floor(start(range(klk12_gl1$e))/1e4),
                               ceiling(end(range(klk12_gl1$e))/1e4)))

klk12_te <- as.data.frame(mcols(klk12_gl1$e)[paste0("s", 1:177)])
klk12_te$eid <- 1:nrow(klk12_te)
klk12_te <- reshape2::melt(klk12_te, id.var="eid")
klk12_te$value <- log10(1+klk12_te$value)

ggplot(klk12_te, aes(xmin=eid-.4, xmax=eid+.4,
                   ymin=value, ymax=value+.05, group=variable)) +
    geom_rect(color="grey60", fill="grey30", alpha=1/5) + guides(color=FALSE) +
    geom_line(aes(x=eid, y=value), alpha=1/5) +
    theme_bw() +
    xlab("exon") + ylab("log10 expr")



#'
#' <br> <br>
#' **adding annotations**
#+ hide=FALSE, fig.height=6, fig.width=10
a <- plotRangesLinkedToData(klk12_gl1$e,
                            stat.y=paste0("s", 1:177),
                            linetype=0, annotation=bb)
b <- autoplot(klk12_gl1$e, color="blue", fill="red", alpha=1/5) 

##clean up top plot
a@grobs[[1]]@ggplot <-
    a@grobs[[1]]@ggplot +
        geom_line(aes(x=x.new, y=value), alpha=1/5) + 
            aes(group=.ggbio.group, alpha=1/10) +
                scale_colour_manual(values=rep("black", 177)) +
                    scale_y_log10() +
                        theme(legend.position="none")
##replace annotation plot
a@grobs[[3]] <- b

##print to screen
a




#'
#' <br> <br> <br> <hr> <br>
#' 
#' ### PCA exon expression (one CC)
#' 

#'
#' <br> <br>
#' looking at clustering along the exons at KLK12
#+ hide=FALSE
matplot(log10(1+klk12_ej[klk12_ej$kind == "e", -c(1:6, 184)]), type="l", lty=1)


##PCA decomposition (exons separately)
pc <- prcomp(log10(1+t(klk12_ej[klk12_ej$kind == "e", -c(1:6, 184)])))

#'
#' <br> <br>
#' exon PCA decomposition (scores)
#+ hide=FALSE, fig.width=8, fig.height=8
plot(pc$x[, 1:2], pch=16, col=brew_set1[2])


#'
#' <br> <br>
#' exon PCA decomposition (loadings)
#+ hide=FALSE, fig.height=3, fig.weight=10
qplot(data=reshape2::melt(pc$rotation[, 1:2]), x=rep(1:ncol(pc$rotation), 2), y=value, color=Var2,
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
klk12_te$cluster <- rep(groups, each=max(klk12_te$eid))
ggplot(klk12_te, aes(xmin=eid-.4, xmax=eid+.4,
                   ymin=value, ymax=value+.05, group=variable,
                   color=as.factor(cluster), fill=as.factor(cluster))) +
    geom_rect(color="grey90", alpha=1/10) + guides(color=FALSE) +
    geom_line(aes(x=eid, y=value), alpha=1/3) +
    guides(color=FALSE, fill=FALSE) +
    theme_bw() +
    xlab("exon") + ylab("log10 expr")



## #'
## #' <br> <br>
## #' **adding annotations**
## #+ hide=FALSE, fig.height=6, fig.width=10
## a <- plotRangesLinkedToData(klk12_gl1$e,
##                             stat.y=paste0("s", 1:177),
##                             linetype=0, annotation=bb)
## b <- autoplot(klk12_gl1$e, color="blue", fill="red", alpha=1/5) 

## ##clean up top plot
## a@grobs[[1]]@ggplot <-
##     a@grobs[[1]]@ggplot +
##         geom_line(aes(x=x.new, y=value), alpha=1/5) + 
##             aes(group=.ggbio.group, alpha=1/10) +
##                 scale_colour_manual(values=rep("black", 177)) +
##                     scale_y_log10() +
##                         theme(legend.position="none")
## ##replace annotation plot
## a@grobs[[3]] <- b

## ##print to screen
## a



## #'
## #' <br> <br>
## #' looking at clusters with all exons having equal width,
## #' onlying looking at first 11 exons
## #+ hide=FALSE, fig.height=3, fig.width=10
## ggplot(subset(klk12_te, eid <= 11),
##        aes(xmin=eid-.4, xmax=eid+.4,
##            ymin=value, ymax=value+.05,
##            group=variable,
##            fill=as.factor(cluster),
##            color=as.factor(cluster))) +
##     geom_rect(color="grey90", alpha=1/10) + guides(color=FALSE) +
##     geom_line(aes(x=eid, y=value), alpha=1/3) +
##     scale_x_continuous("exon", breaks=1:11) +
##     theme_bw() +
##     ylab("log10 expr")




#'
#' <br> <br> <br> <hr> <br>
#' 
#' ### Plotting junction coverage (one CC)
#' 

##PCA decomposition (junctions separately)
pc_j <- prcomp(log10(1+t(klk12_ej[klk12_ej$kind == "j", -c(1:6, 184)])))


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
klk12_tj <- as.data.frame(mcols(klk12_gl1$j)[paste0("s", 1:177)])
klk12_tj$eid <- 1:nrow(klk12_tj)
klk12_tj <- reshape2::melt(klk12_tj, id.var="eid")
klk12_tj$value <- log10(1+klk12_tj$value)
klk12_tj$cluster <- rep(groups, each=max(klk12_tj$eid))

ggplot(klk12_tj, aes(xmin=eid-.4, xmax=eid+.4,
                   ymin=value, ymax=value+.05, group=variable,
                   color=as.factor(cluster), fill=as.factor(cluster))) +
    geom_rect(color=NA, alpha=1/5) + guides(color=FALSE) +
    guides(color=FALSE, fill=FALSE) +
    theme_bw() +
    xlab("exon") + ylab("log10 expr")

ggplot(klk12_tj, aes(xmin=eid-.4, xmax=eid+.4,
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

klk12_v <- klk12_ej[klk12_ej$kind == "e", c("start", "stop")]
klk12_e <- klk12_ej[klk12_ej$kind == "j", c("start", "stop")]

## determine explicit edgeset
klk12_e$a <- sapply(klk12_e$start, function(x) which(klk12_v$stop == x))
klk12_e$b <- sapply(klk12_e$stop, function(x) which(klk12_v$start == x))
klk12_el1 <- as.matrix(klk12_e[c("a", "b")])

## determine 'edges' that exist between consecutive exons
conseq <- which((klk12_v$start[-1] - klk12_v$stop[-nrow(klk12_v)]) == 1)
klk12_el2 <- cbind("a"=conseq, "b"=conseq+1)


#'
#' <br> <br>
#' **Plot splicing graph for Sample 7**
#+ hide=FALSE, fig.width=8, fig.height=14

klk12_vw2 <- klk12_ej[klk12_ej$kind == "e", "s7"]
klk12_ew2 <- klk12_ej[klk12_ej$kind == "j", "s7"]

klk12_adj2 <- matrix(0, nrow(klk12_v), nrow(klk12_v))
diag(klk12_adj2) <- klk12_vw2
klk12_adj2[klk12_el1] <- klk12_ew2
klk12_adj2[klk12_el2] <- apply(matrix(diag(klk12_adj2)[klk12_el2], ncol=2), 1, min)

g2 <- graph.adjacency(klk12_adj2, weighted=TRUE, diag=FALSE)

edge_h <- 1/apply(get.edgelist(g2), 1, diff)
edge_h[edge_h == 1] <- 1e-10

plot(g2, layout=cbind(1:13, 0),
     edge.curved=7*edge_h,
     edge.arrow.size=.3,
     edge.width=log2(E(g2)$weight+1),
     vertex.shape="rectangle", vertex.label=NA,
     vertex.size=5.5, vertex.size2=log2(V(g2)$w+1))



#'
#' <br> <br>
#' **example of adjacency matrix**
#+ hide=FALSE, fig.width=8, fig.height=14
get.adjacency(graph.adjacency(klk12_adj2))






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

