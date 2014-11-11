---
title: "Splice Graph Heatmaps (KLK12)"
output:
  html_document:
    toc: true
    fig_width: 12
    fig_height: 5
    fig_caption: false
  md_document:
    toc: true
    fig_width: 12
    fig_height: 5
---




basic plot



```r
heatgraph(klk12_gr1)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png) 


basic plot, don't flip on strand



```r
heatgraph(klk12_gr1, flip_neg=FALSE)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) 


use non-genomic spacing



```r
heatgraph(klk12_gr1, genomic=FALSE)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 


use black background



```r
heatgraph(klk12_gr1, genomic=FALSE, use_blk=TRUE)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png) 


include splicing information



```r
heatgraph(klk12_gr1, genomic=FALSE, j_incl=TRUE)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png) 


include splicing information (black background)



```r
heatgraph(klk12_gr1, genomic=FALSE, use_blk=TRUE, j_incl=TRUE)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png) 


include splicing information (black background) without flipping



```r
heatgraph(klk12_gr1, genomic=FALSE, use_blk=TRUE, j_incl=TRUE, flip_neg=FALSE)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png) 


plot on continuous scale



```r
heatgraph(klk12_gr1, bin=FALSE, genomic=FALSE, j_incl=TRUE)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png) 


bin using different log base



```r
heatgraph(klk12_gr1, log_base=2, genomic=FALSE, j_incl=TRUE)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png) 


sort each box separately



```r
heatgraph(klk12_gr1, sort_sep=TRUE, genomic=FALSE, j_incl=TRUE)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png) 


add annotations



```r
##lbls <- c(rep(1, 50), rep(2, 30), rep(3, 40),
##          rep(1, 20), rep(3, 10), rep(2, 27))
##heatgraph(klk12_gr1, genomic=FALSE, use_blk = FALSE,
##          highligh=lbls)
```


plot PCA of connected component with exon and junctions considered separately



```r
graphPCA(klk12_gr1)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png) 


plot PCA not using genomic coordinates



```r
graphPCA(klk12_gr1, genomic=FALSE)
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png) 


compute PCA using exon and junction values jointly



```r
graphPCA(klk12_gr1, pc_sep=FALSE, genomic=FALSE)
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15-1.png) 

