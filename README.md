

spliceclust [![Build Status](https://travis-ci.org/pkimes/spliceclust.svg?branch=master)](https://travis-ci.org/pkimes/spliceclust)
=======================

## Contents
1. [Introduction](#intro)
2. [SpliceGraHM Examples](#splicegrahm)
3. [SplicePCA Examples](#splicepca)
4. [SplicePCP Examples](#splicepcp)
4. [SpliceGraHM2 Examples](#splicegrahm2)


## <a name="intro"></a> Introduction
___Note that this document is still under construction.___  

This package may be used to plot exon and splice junction coverage across a large cohort
of RNA-seq samples. The plotting approach is based on the idea of transposing
expression heatmaps on splicing diagrams. The plots are generated using the `ggplot2` and
`ggbio` packages.



## <a name="splicegrahm"></a> SpliceGraHM Examples

First we must construct a `concomp` object from a `GRangesList` object with exon and junction
boundaries and corresponding coverage values. For illustration purposes, we use a previously
constructed `GRanges` object, `klk12_gr1`, containing annotation information for the KLK12 gene
locus along with coverage for 177 lung squamous cell carcinoma samples.  

The exon and junction boundaries, and the corresponding kind label are given by:  


```r
ranges(klk12_gr1)
```

```
## IRanges of length 25
##         start      end width names
## [1]  51532348 51532468   121 33127
## [2]  51532468 51532598   131 33128
## [3]  51532469 51532597   129 33129
## [4]  51532598 51532713   116 33130
## [5]  51532713 51534044  1332 33131
## ...       ...      ...   ...   ...
## [21] 51537896 51538051   156 33147
## [22] 51537896 51538062   167 33148
## [23] 51537897 51538040   144 33149
## [24] 51538051 51538061    11 33150
## [25] 51538062 51538261   200 33151
```

```r
mcols(klk12_gr1)$kind
```

```
##  [1] "e" "j" "e" "e" "j" "e" "j" "j" "j" "e" "e" "j" "e" "j" "j" "j" "j"
## [18] "e" "e" "e" "j" "j" "e" "e" "e"
```

In addition to `kind`, the `klk12_gr1` metadata columns also include the gene name (`gIdx`),
gene boundaries (`gStart`, `gStop`), and coverage values for the 177 samples.  


```r
mcols(klk12_gr1)[1:5, 1:6]
```

```
## DataFrame with 5 rows and 6 columns
##          gIdx    gStart     gStop        kind        s1        s2
##   <character> <numeric> <numeric> <character> <numeric> <numeric>
## 1    gene9317  51532348  51538261           e   2.60331   5.32231
## 2    gene9317  51532348  51538261           j   0.00000   0.00000
## 3    gene9317  51532348  51538261           e   2.59690  10.37210
## 4    gene9317  51532348  51538261           e   2.15517   7.28448
## 5    gene9317  51532348  51538261           j   0.00000   7.00000
```

To construct the `concomp` (connected component) object for analysis, we first convert the
`GRanges` object into a `GRangesList` object of length 2, corresponding to exon and junction
information.  


```r
klk12_gl <- split(klk12_gr1, mcols(klk12_gr1)$kind)
klk12_cc <- concomp(klk12_gl)
```


We first demonstrate the default and __basic SpliceGraHM__ (Splice Graph Heat Map) plotting procedure.


```r
splicegrahm(klk12_cc)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

In the plot above, each box arranged horizontally corresponds to a contiguous exonic region
along the genome. Each box is colored by 177 horizontal lines, showing the expression level
for the 177 samples being analyzed. Note that the exons are plotted along genomic coordinates, and
log expression is shown in color.  

In addition to the boxes, the plot contains arrows which correspond to a splicing events with
sufficient support in the data (e.g. at least 8 samples, each with at least 5 reads spanning the splice
junction). The arrows are colored such that darker arrows were present in a higher proportion of the
177 samples (scale shown on right).  

The SpliceGraHM name is in reference to the fact that after removing the spacing between each exon
the plot simply reduces to a __standard heatmap__ of expression along a single gene. The corresponding
figure after removing the intronic gaps is shown below. The rectangle of color is a heatmap with rows
and columns corresponding to samples and exons, respectively.  


```r
splicegrahm(klk12_cc, genomic = FALSE, ex_use = 1, log_base = 2)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

It is possible to show the coverage of each splice junction using a similar convention with rows
corresponding to samples, and with color being used for coverage.


```r
splicegrahm(klk12_cc, genomic = FALSE, j_incl = TRUE)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

If annotation information is available, e.g. from UCSC KnownGenes, these can be passed to the
function to add an additional track to the plot. The appropriate `GRangesList` object to be
passed is illustrated in the following example.


```r
suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene"))
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
isActiveSeq(txdb)[seqlevels(txdb)] <- FALSE
isActiveSeq(txdb)[paste0("chr", 1:22)] <- TRUE
exbytx <- exonsBy(txdb, "tx")

splicegrahm(klk12_cc, genomic = TRUE, j_incl = TRUE, txlist = exbytx)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

If gene names are desired, the following can be used to match the transcript ID
in `txdb` against gene symbols (e.g. in `org.Hs.eg.db`).


```r
suppressPackageStartupMessages(library("org.Hs.eg.db"))

splicegrahm(klk12_cc, genomic = TRUE, j_incl = TRUE, txlist = exbytx,
            txdb = txdb, orgdb = org.Hs.eg.db)
```

```
## 'select()' returned 1:1 mapping between keys and columns
```

```
## 'select()' returned many:1 mapping between keys and columns
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)

The `splicegrahm` function can now also plot gene models in non-genomic space
with an additional parameter `eps`. The `eps` parameter determines how far up/down
from the connected component to lookfor overlapping gene models. If `eps = NULL`,
all overlapping gene models are included. If `eps = 1000`, only overlapping gene
models which are fully contained within 1000bp of the connected component
are included.


```r
splicegrahm(klk12_cc, genomic = FALSE, j_incl = TRUE, txlist = exbytx,
            txdb = txdb, orgdb = org.Hs.eg.db, eps = NULL)
```

```
## 'select()' returned 1:1 mapping between keys and columns
```

```
## 'select()' returned many:1 mapping between keys and columns
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)


```r
splicegrahm(klk12_cc, genomic = FALSE, j_incl = TRUE, txlist = exbytx,
            txdb = txdb, orgdb = org.Hs.eg.db, eps = 0)
```

```
## 'select()' returned 1:1 mapping between keys and columns
```

```
## 'select()' returned many:1 mapping between keys and columns
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)


```r
splicegrahm(klk12_cc, genomic = FALSE, j_incl = TRUE, txlist = exbytx,
            txdb = txdb, orgdb = org.Hs.eg.db, eps = 1000)
```

```
## 'select()' returned 1:1 mapping between keys and columns
```

```
## 'select()' returned many:1 mapping between keys and columns
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)


## <a name="splicepca"></a> SplicePCA Examples

Principal Component Analysis (PCA) is a popular exploratory analysis tool for studying low rank
structure in high-dimensional datasets by projecting the data along the directions of greatest
variation. These directions of greatest variation are often referred to as the PC "loading"
directions. The `splicepca` function is written to visualize the loading vectors of splicing
data.  

As an example, the following code can be used to visualize the first 3 PC loadings in the
KLK12 dataset from above.


```r
splicepca(klk12_cc, npc = 3)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png)

In the above PC loadings plot, red and blue are used to denote the magnitude of the PC loadings
for each exon and splice junction. From the plot above, we see that the greatest source of variation
in expression at the KLK12 gene locus is the overall expression of the gene. However, the second
PC loading shows interesting behavior where the primary source of variation is attributable to the
expression of the central exon. Scrolling up a little to the UCSC KnownGenes shown included above,
we see that the presence or absence of the central exon actually corresponds to differing isoforms
annotated to this gene. Note that above the PCs are computed separately for exons and junctions
 information. As such, we obtain separate PC "scores" for the exon and junction PC loadings.
 

```r
## splicepca(klk12_cc, npc = 3, scores = TRUE)
```

It is also possible to perform the PCA analysis using the concatenated exon and junction information by
setting the `pc_sep` parameter to `FALSE`, and specifying the relative "weight" of each with `ej_w`. The
exon and junction data are rescaled to each have sum-of-squares equal to the values specified by `ej_w`.
In the following example, we use equal weights for the two data sources.


```r
splicepca(klk12_cc, npc = 3, pc_sep = FALSE, ej_w = c(1, 1))
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15-1.png)


## <a name="splicepcp"></a> SplicePCP Examples

All above plots have used color to represent expression. However, low frequency or outlier events may
become lost in this particular view. To handle this problem, we also provide the option to plot the
splicing objects with vertical height (the y-axis) corresponding to expression. The name of the function is
a reference to [_parallel coordinates plots_](pcp). The same _KLK12_ data is shown below using the
`splicepcp` function.


```r
splicepcp(klk12_cc, genomic = TRUE, txlist = exbytx, txdb = txdb, orgdb = org.Hs.eg.db)
```

```
## 'select()' returned 1:1 mapping between keys and columns
```

```
## 'select()' returned many:1 mapping between keys and columns
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)

The plot includes 3 tracks:
  1. log-expression values for each exon or region of an exon
  2. the complete gene model and splice junctions
  3. annotated transcripts from the UCSC KnownGene database.

Currently, the function is being rewritten to also include junction coverage using a
parallel coordinates plot in a separate track.


## <a name="splicegrahm2"></a> SpliceGraHM2 Examples

To make direct comparison of two populations possible, we have created `splicegrahm2` which
draws two `splicegrahm` plots in a single figure. Consider, for example, the task of comparing
the behavior of splicing between the LUSC samples and head and neck squamous cell carcinoma (HNSC)
samples. We have loaded a connected component for 280 HNSC samples as `klk12_cc_hnsc`.


```r
splicegrahm2(klk12_cc, klk12_cc_hnsc, genomic=FALSE)
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17-1.png)

When gene models are also desired, they are placed in a central track for easy comparison to
the two `splicegrahm`s.


```r
splicegrahm2(klk12_cc, klk12_cc_hnsc, genomic=FALSE,
             txlist = exbytx, txdb = txdb, orgdb = org.Hs.eg.db)
```

```
## 'select()' returned 1:1 mapping between keys and columns
```

```
## 'select()' returned many:1 mapping between keys and columns
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18-1.png)

Alternatively, the bottom plot can be drawn with the same orientation as the top plot by
setting `mirror=FALSE`.


```r
splicegrahm2(klk12_cc, klk12_cc_hnsc, genomic=FALSE, mirror=FALSE)
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-1.png)

This plotting scheme may be useful for comparing two subpopulations from a single dataset.
Suppose we perform clustering on the LUSC dataset using the `cluster` function implemented
in the `spliceclust` package. Note that the two clusters differ in their usage of the
central exon, with a clear difference in expression and splicing at this exon.


```r
lbls <- cluster(klk12_cc)

##create connected component with just cluster 2 samples
klk12_cc2 <- klk12_cc
exonValues(klk12_cc2) <- exonValues(klk12_cc)[, lbls$junc_labs == 2]
juncValues(klk12_cc2) <- juncValues(klk12_cc)[, lbls$junc_labs == 2]

##create connected component with just cluster 3 samples
klk12_cc3 <- klk12_cc
exonValues(klk12_cc3) <- exonValues(klk12_cc)[, lbls$junc_labs == 3]
juncValues(klk12_cc3) <- juncValues(klk12_cc)[, lbls$junc_labs == 3]

splicegrahm2(klk12_cc2, klk12_cc3, genomic=FALSE, sort_idx1=6, sort_idx2=6)
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20-1.png)

While not modified above, the option ``same_scale`` can be used to specify whether
the two plots should be drawn using the same or separate scale along the y-axis.
This can be helpful when the two populations are of substantially different sizes.



## <a name="sessioninfo"></a> Session Information


```r
sessionInfo()
```

```
## R version 3.2.2 (2015-08-14)
## Platform: x86_64-apple-darwin15.0.0 (64-bit)
## Running under: OS X 10.11.3 (El Capitan)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] org.Hs.eg.db_3.2.3                     
##  [2] RSQLite_1.0.0                          
##  [3] DBI_0.3.1                              
##  [4] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
##  [5] GenomicFeatures_1.22.13                
##  [6] spliceclust_0.1.2                      
##  [7] BSgenome.Hsapiens.UCSC.hg19_1.4.0      
##  [8] BSgenome_1.38.0                        
##  [9] rtracklayer_1.30.2                     
## [10] Biostrings_2.38.4                      
## [11] XVector_0.10.0                         
## [12] GGally_1.0.1                           
## [13] annotate_1.48.0                        
## [14] XML_3.98-1.3                           
## [15] AnnotationDbi_1.32.3                   
## [16] Biobase_2.30.0                         
## [17] GenomicRanges_1.22.4                   
## [18] GenomeInfoDb_1.6.3                     
## [19] IRanges_2.4.7                          
## [20] S4Vectors_0.8.11                       
## [21] RColorBrewer_1.1-2                     
## [22] ggbio_1.18.5                           
## [23] BiocGenerics_0.16.1                    
## [24] ggplot2_2.0.0                          
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.3                biovizBase_1.18.0         
##  [3] lattice_0.20-33            Rsamtools_1.22.0          
##  [5] digest_0.6.9               plyr_1.8.3                
##  [7] futile.options_1.0.0       acepack_1.3-3.3           
##  [9] evaluate_0.8               BiocInstaller_1.20.1      
## [11] zlibbioc_1.16.0            rpart_4.1-10              
## [13] labeling_0.3               devtools_1.10.0           
## [15] splines_3.2.2              BiocParallel_1.4.3        
## [17] stringr_1.0.0              foreign_0.8-66            
## [19] RCurl_1.95-4.7             biomaRt_2.26.1            
## [21] munsell_0.4.2              nnet_7.3-12               
## [23] SummarizedExperiment_1.0.2 gridExtra_2.0.0           
## [25] roxygen2_5.0.1             Hmisc_3.17-1              
## [27] reshape_0.8.5              withr_1.0.1               
## [29] GenomicAlignments_1.6.3    bitops_1.0-6              
## [31] grid_3.2.2                 RBGL_1.46.0               
## [33] xtable_1.8-2               gtable_0.1.2              
## [35] magrittr_1.5               formatR_1.2.1             
## [37] scales_0.3.0               graph_1.48.0              
## [39] stringi_1.0-1              reshape2_1.4.1            
## [41] latticeExtra_0.6-28        futile.logger_1.4.1       
## [43] Formula_1.2-1              lambda.r_1.1.7            
## [45] tools_3.2.2                dichromat_2.0-0           
## [47] OrganismDbi_1.12.1         survival_2.38-3           
## [49] colorspace_1.2-6           cluster_2.0.3             
## [51] memoise_1.0.0              knitr_1.12.3              
## [53] VariantAnnotation_1.16.4
```


[pcp]: http://en.wikipedia.org/wiki/Parallel_coordinates
