#' ---
#' title: "Research Update"
#' author: "Patrick Kimes"
#' output:
#'   ioslides_presentation:
#'     fig_width: 4
#'     fig_height: 4
#'     smaller: true
#'     widescreen: true
#' ---

#+ hide=TRUE, echo=FALSE
knitr::opts_chunk$set(collapse=TRUE, comment="##",
                      echo=FALSE, warning=FALSE, error=FALSE,
                      hide=TRUE, message=FALSE)
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

## set nicer colors
brew_set1 <- brewer.pal(9, "Set1")
brew_pastel2 <- brewer.pal(8, "Pastel2")
brew_accent <- brewer.pal(8, "Accent")
brew_greys <- brewer.pal(9, "Greys")

## specify project directory and source functions
root <-"/Users/pkimes/Dropbox/Git/spliceclust/"
source(paste0(root, "R/readchr.R"))

## load chromesome 6 data for 177 lusc samples
## complete datasest
chr9 <- readchr(paste0(root, "data/lusc/chr9_gene.txt"), 177)
## exon only dataset
chr9_e <- subset(chr9, kind == "e")


#'
#' # Overview
#'


#' 
#' ## Context {.build}
#'
#' <div class="columns-2">
#' <img src="images/examples_CDKN2A.png" width="350px" />
#'  
#' - Building on `SigFuge`
#'     - per-gene analysis
#'     - cluster by per-base expression
#' -
#' -
#' - extend per-gene clustering
#'     - better annotations
#'     - incorporate splicing
#' </div>
#'


#'
#' ## Connected Components {.flexbox .vcenter}
#'
#' <img src="images/concomp1.png" width="600px" />
#' 
#' - joint work with Yin Hu, Jan Prins
#' -
#' - data driven annotations
#' - graphs: nodes (exons), edges (splice junctions)
#' -
#' - requires specifying parameters
#'     - `thresh_junctionfilter_num_present_samples` (8)
#'     - `thresh_junctionfilter_present_support` (5)
#'     - `coverageThreshold_exon` (.5)
#'     - `coverageThreshold_intron` (10)
#'


#'
#' ## Connected Components {.flexbox .vcenter}
#'
#+ hide=FALSE
chr9[107:116, 1:9]



#'
#' #LUSC chr9 Analysis
#'


#'
#' ## Lengths of connected components {.flexbox .vcenter}
#' 
#+ hide=FALSE, fig.height=5, fig.width=7
exon_cnts <- table(chr9_e$gIdx)
qplot(x=Freq, data=data.frame(exon_cnts), binwidth=I(1),
      color=I(brew_greys[5]), fill=I(brew_pastel2[1])) +
    theme_bw() +
    scale_y_continuous(trans=scales::log1p_trans(),
                       breaks=c(0, 10^(0:5))) +
    xlab("number of exons") +
    ylab("log(1 + #models) , real counts shown") +
    ggtitle("distribution of exon counts for connected components")



#'
#' ## Distribution of median exon expression
#' 
#' compute distribution of average expression over each exon and plot distributions
#' grouped by the number of exons each gene
#+ hide=FALSE, echo=FALSE, fig.height=5, fig.width=10
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
#' ## Location of single exon components {.flexbox .vcenter}
#' 
## subset on only single gene models
e1_gIdx <- names(which(exon_cnts == 1))
e1set <- subset(chr9_e, gIdx %in% e1_gIdx)
e1set$chr <- "chr9"
e1set_gr <- makeGRangesFromDataFrame(e1set,
                                     seqnames.field="chr",
                                     keep.extra.columns=TRUE)
seqlengths(e1set_gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)["chr9"]
## plot location of single exon genes on chromosome
## load hg19 cytoband information
data(ideoCyto, package = "biovizBase")
biovizBase::isIdeogram(ideoCyto$hg19)
## construct ideogram plot for comparison
gg_chr9 <- autoplot(subset(ideoCyto$hg19, seqnames(ideoCyto$hg19) == "chr9"),
                    layout="karyogram", cytoband=TRUE) +
                        guides(fill=FALSE)
fixed(gg_chr9) <- TRUE
## construct karyogram with chr19 with 
gg_e1 <- autoplot(e1set_gr, layout="karyogram", alpha=1/100)
fixed(gg_e1) <- TRUE
#'
#+ hide=FALSE, fig.height=2, fig.width=10
##combine plots
tracks("chr" = gg_chr9, "e1 genes" = gg_e1,
       heights=c(.5, .5), xlab="chr9 positions") +
    theme_tracks_sunset()


#'
#' ## ...contrast with UCSC KnownGenes {.flexbox .vcenter}
#' 
## load transcript annotations and only look at chr9
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
isActiveSeq(txdb)[seqlevels(txdb)] <- FALSE
isActiveSeq(txdb)["chr9"] <- TRUE
## want to work at exon level aggregated by genes or transcripts
exbygene <- exonsBy(txdb, "gene")
exbytx <- exonsBy(txdb, "tx")
##contrast with UCSC
gg_ucsc <- autoplot(exbygene)
fixed(gg_ucsc) <- TRUE
#+ hide=FALSE, fig.height=3, fig.width=10
tracks("chr" = gg_chr9,
       "e1 genes" = gg_e1,
       "UCSC KG" = gg_ucsc,
       heights=c(1/3, 1/3, 1/3), xlab="chr9 positions") +
    theme_tracks_sunset()


#'
#' ## ...contrast with UCSC KnownGenes
#' 
## Calculate overlap with genes and transcripts
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
#' - `r length(queryHits(hits_g))` overlaps
#' - overlap consisting of:
#'     - `r length(unique(queryHits(hits_g)))` UCSC genes
#'     - `r length(unique(subjectHits(hits_g)))` connected components
#' 
#' <div class="centered">
#+ hide=FALSE, fig.width=8, fig.height=4
ovlp_table_g <- data.frame(table(ucsc_ovlp_g$SYMBOL))
ggplot(data=ovlp_table_g) +
    geom_text(aes(label=Var1, x=1:length(Var1), y=Freq)) +
    theme_bw() +
    xlab("") +
    ylab("# overlapping components") +
    ggtitle("UCSC KnownGenes overlapping single gene connected components")
#' </div>


#' 
#' ## Example of overlap
#' 
#' <div class="centered">
#+ hide=FALSE, fig.height=2, fig.width=10
top1 <- ovlp_table_g$Var1[which.max(ovlp_table_g$Freq)]
top1_tx <- subset(ucsc_ovlp_tx, ucsc_ovlp_tx$SYMBOL == top1)

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
#' </div>
#' 
#' - Gene: `r ovlp_table_g$Var1[which.max(ovlp_table_g$Freq)]`
#' -  overlapped with `r max(ovlp_table_g$Freq)` single exon connected components.


#'
#' # CDKN2A (p16/p14) locus
#'


#'
#' ## UCSC Genome Browser view {.flexbox .vcenter}
#'
#' <img src="images/ucsc_p16.png" width="600px" />  
#'   
#' - ucsc browser for chr9:21,965,000-21,995,000
#'


#'
#' ## Connected components in region {.flexbox .vcenter}
#' 
## look in region around cdkn2a (p16/p14)  (on rev strand)
p16_eid <- select(org.Hs.eg.db, keys="CDKN2A",
                  columns="ENTREZID", keytype="SYMBOL")
p16_txid <- select(txdb, keys=p16_eid$ENTREZID,
                   columns=c("TXID", "TXNAME"), keytype="GENEID")
p16_tx <- exbytx[as.character(p16_txid$TXID)]
p16_uniontx <- reduce(unlist(p16_tx))
bounds <- c(start(range(p16_uniontx)), end(range(p16_uniontx)))

cands1 <- (chr9_e$start > min(bounds) - 10000) &
          (chr9_e$stop < max(bounds) + 10000)
cand_g <- chr9_e$gIdx[cands1]

p16_set <- subset(chr9_e, gIdx %in% unique(cand_g))
p16_set$chr <- "chr9"

## convert to `GRanges` object to use with `ggbio`
p16_gr <- makeGRangesFromDataFrame(p16_set,
                                   seqnames.field="chr",
                                   keep.extra.columns=TRUE)
seqlengths(p16_gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)["chr9"]

## convert to `GRangesList` object to separate out groups of genes
p16_gl <- split(p16_gr, mcols(p16_gr)$gIdx)

gg_p16models <- autoplot(p16_tx)
fixed(gg_p16models) <- FALSE

gg_p16_cc <- autoplot(p16_gl)
fixed(gg_p16_cc) <- FALSE

#+ hide=FALSE, fig.height=2, fig.width=10
tracks("UCSC" = gg_p16models,
       "concomp" = gg_p16_cc,
       heights=c(1/3, 1/3), xlab="chr9 positions") +
    theme_tracks_sunset()



#'
#' ## Exon level expression for gene1791 {.flexbox .vcenter}
#' 
## From above, most likely connected component was `gene1791`
p16_ej <- subset(chr9, gIdx == "gene1791")
p16_ej$chr <- "chr9"
p16_gr1 <- makeGRangesFromDataFrame(p16_ej,
                                    seqnames.field="chr",
                                    keep.extra.columns=TRUE)
seqlengths(p16_gr1) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)["chr9"]
strand(p16_gr1) <- "-"
p16_gl1 <- split(p16_gr1, mcols(p16_gr1)$kind)

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

#+ hide=FALSE, fig.height=5, fig.width=10
ggplot(p16_te, aes(xmin=eid-.4, xmax=eid+.4,
                   ymin=value, ymax=value+.05, group=variable)) +
    geom_rect(color="grey60", fill="grey30", alpha=1/5) + guides(color=FALSE) +
    geom_line(aes(x=eid, y=value), alpha=1/5) +
    theme_bw() +
    xlab("exon") + ylab("log10 expr")



#'
#' ## ...adding annotations {.flexbox .vcenter}
#'
#+ hide=FALSE, fig.height=6, fig.width=10
plotRangesLinkedToData(p16_gl1$e, stat.y=paste0("s", 1:177),
                       linetype=0, annotation=bb)


#'
#' ## Genomic coordinates {.flexbox .vcenter}
#' 
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

#+ hide=FALSE, fig.height=5, fig.width=10
tracks(gg_logexp,
       autoplot(GRangesList(p16_gl1$e), gap.geom="arrow",
                fill="grey80", color="grey30") + theme_bw(),
       "UCSC" = gg_p16models,
       "concomp" = gg_p16_cc,
       heights=c(3, 1, 1, 1))



#'
#' # Data Object 
#' 


#'
#' ## Data Object {.build}
#'
#' > - have only really looked at exon information
#' > - also have coverage of junctions
#' > -
#' > - how to study?
#' >     - separately and aggregate?
#' >     - graph structure?
#' > - 
#' > - example using `gene1791` component
#' 


## PCA decomposition (exon and splicing separately)
pc <- prcomp(log10(1+t(p16_ej[p16_ej$kind == "e", -c(1:6, 184)])))
pc_j <- prcomp(log10(1+t(p16_ej[p16_ej$kind == "j", -c(1:6, 184)])))


#'
#' ## Standard PCA, exon expression {.flexbox .vcenter}
#'
#+ hide=FALSE, fig.height=5.5, fig.width=5.5
plot(pc$x[, 1:2], pch=16, col=brew_set1[2])


#'
#' ## Standard PCA, exon expression {.flexbox .vcenter}
#'
#+ hide=FALSE, fig.height=3, fig.width=10
qplot(data=reshape2::melt(pc$rotation[, 1:2]), x=rep(1:22, 2), y=value, color=Var2,
      geom="line") +
    theme_bw() +
    scale_x_continuous(breaks=1:22)


#' look at clusters
groups <- as.numeric(pc$x[, 1]>-1)*(1+as.numeric(pc$x[, 2]>0)) + 1


#'
#' ## Standard PCA, exon expression {.flexbox .vcenter}
#'
#+ hide=FALSE, fig.height=5.5, fig.width=5.5
plot(pc$x[, 1:2], pch=16, col=brew_set1[groups])


#'
#' ## Standard PCA, junction coverage {.flexbox .vcenter}
#'
#+ hide=FALSE, fig.height=5.5, fig.width=5.5
plot(pc_j$x[, 1:2], pch=16, col=brew_set1[2],
     main="PCA scores for splicing junctions")


#'
#' ## Graph based approaches
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

p16_vw2 <- p16_ej[p16_ej$kind == "e", "s2"]
p16_ew2 <- p16_ej[p16_ej$kind == "j", "s2"]

p16_adj2 <- matrix(0, nrow(p16_v), nrow(p16_v))
diag(p16_adj2) <- p16_vw2
p16_adj2[p16_el1] <- p16_ew2
p16_adj2[p16_el2] <- apply(matrix(diag(p16_adj2)[p16_el2], ncol=2), 1, min)

g2 <- graph.adjacency(p16_adj2, weighted=TRUE, diag=FALSE)

edge_h <- 1/apply(get.edgelist(g2), 1, diff)
edge_h[edge_h == 1] <- 1e-10

#+ hide=FALSE, fig.width=12, fig.height=8
plot(g2, layout=cbind(1:22, 0),
     edge.curved=7*edge_h,
     edge.arrow.size=.3,
     edge.width=log2(E(g2)$weight+1),
     vertex.shape="rectangle", vertex.label=NA,
     vertex.size=4.5, vertex.size2=2)


#'
#' ## Graph based approaches
#'
#+ hide=FALSE
get.adjacency(graph.adjacency(p16_adj2))[1:15, 1:15]
#'
#' - example of an adjacency matrix for ONE sample
#' - 177 samples leads to a 3-dimensional tensor
#'     - tensor factorizations
#'     - spectral methods


#'
#' ## Ongoing work
#'
#' - understanding connected components better
#'     - hard rules
#'     - looking at mismatch density, multimapping, ...
#' -
#' - developing suitable analysis methods
#'     - clustering graph objects
#'     - explore alternative formulations
#' 












#' 
#' # `R` Environment
#'


#' 
#' ## `R` Environment
#'
#' Running time:
#+ echo=TRUE, hide=FALSE
proc.time() - start_time
#'
#' Compiled on:
#+ echo=TRUE, hide=FALSE
Sys.time()


#' 
#' ## `R` Environment
#'
#' Session Information:
#+ echo=TRUE, hide=FALSE
sessionInfo()

