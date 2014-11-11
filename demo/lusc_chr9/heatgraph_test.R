#' ---
#' title: "Splice Graph Heatmaps (CDKN2A)"
#' output:
#'   html_document:
#'     toc: true
#'     fig_width: 12
#'     fig_height: 5
#'     fig_caption: false
#'   md_document:
#'     toc: true
#'     fig_width: 12
#'     fig_height: 5
#' ---

#+ echo=FALSE, message=FALSE, warning=FALSE, hide=TRUE
knitr::opts_chunk$set(collapse=TRUE, comment="##",
                      echo=TRUE, warning=FALSE,
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


##specify project directory and source functions
root <-"/Users/pkimes/Dropbox/Git/spliceclust/"
source(paste0(root, "R/readchr.R"))
source(paste0(root, "R/heatgraph.R"))
source(paste0(root, "R/heatgraphPCA.R"))


## load CDKN2A component
load(paste0(root, "demo/lusc_chr9/p16_gr1.rdata"))


#'
#' basic plot
#' 
heatgraph(p16_gr1)


#'
#' use non-genomic spacing
#' 
heatgraph(p16_gr1, genomic=FALSE)


#'
#' use non-genomic spacing
#' 
heatgraph(p16_gr1, genomic=FALSE, ex_use=9/10)


#'
#' use black background
#' 
heatgraph(p16_gr1, genomic=FALSE, use_blk=TRUE)


#'
#' include splicing information
#' 
heatgraph(p16_gr1, genomic=FALSE, j_incl=TRUE)


#'
#' include splicing information (black background)
#' 
heatgraph(p16_gr1, genomic=FALSE, use_blk=TRUE, j_incl=TRUE)


#'
#' plot on continuous scale
#' 
heatgraph(p16_gr1, bin=FALSE, genomic=FALSE, j_incl=TRUE)


#'
#' bin using different log base
#' 
heatgraph(p16_gr1, log_base=2, genomic=FALSE, j_incl=TRUE)


#'
#' sort each box separately
#' 
heatgraph(p16_gr1, sort_sep=TRUE, genomic=FALSE, j_incl=TRUE)


#'
#' add annotations
#'
##lbls <- c(rep(1, 50), rep(2, 30), rep(3, 40),
##          rep(1, 20), rep(3, 10), rep(2, 27))
##heatgraph(p16_gr1, genomic=FALSE, use_blk = FALSE,
##          highligh=lbls)



#'
#' plot PCA of connected component with exon and junctions considered separately
#'
graphPCA(p16_gr1)


#'
#' plot PCA not using genomic coordinates
#'
graphPCA(p16_gr1, genomic=FALSE)


#'
#' compute PCA using exon and junction values jointly
#'
graphPCA(p16_gr1, pc_sep=FALSE, genomic=FALSE)



