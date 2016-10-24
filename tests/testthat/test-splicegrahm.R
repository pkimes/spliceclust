context("splicegrahm plotting method")
library("GenomicRanges")
library("ggplot2")

## ##############################################################################
## construct example dataset for test cases

## simple gene model overlapping KLK12 on chr19
ichr <- "chr19"
iseqlengths <- 59128983
igStart <- 51029000 + c(0, 2000, 5000, 1000, 1000, 3000)
igStop <- 51029000 + c(1000, 3000, 6000, 2000, 5000, 5000)
ikind <- c(rep("e", 3), rep("j", 3))
istr <- "-"
    
## samples w/ exon/junc coverage rank reversed
samp_cov <- sapply(c(s=0:19),
                   function(x) c(rep(2^((x+1)/2), 3), rep(2^((20-x)/2), 3)),
                   simplify=FALSE)

## construct concomp
simple_df <- data.frame(chr=ichr, seqlengths=iseqlengths,
                        gStart=igStart, gStop=igStop,
                        kind=ikind,
                        strand=istr,
                        samp_cov)
simple_cc <- concomp(simple_df)


## ##############################################################################
## actual test cases

test_that("splicegrahm default call works", {
    ## default plot with simple dataset
    plt <- splicegrahm(simple_cc)
    bld <- ggplot_build(plt)

    ## check that ggplot is produced
    expect_is(plt, "ggplot")

    ## check ranges of plot 
    expect_equal(bld$panel$ranges[[1]]$x.range[1], 51028700)
    expect_equal(bld$panel$ranges[[1]]$x.range[2], 51035300)

    ## tabulate the layer geoms
    layers_base <- sapply(plt$layer, function(x) class(x$geom)[1])

    ## check one GeomPaths per junction
    expect_equal(sum(layers_base == "GeomPath"), length(juncs(simple_cc)))
    
    ## check GeomRects for exon heatmap, exon outlines
    expect_equal(sum(layers_base == "GeomRect"), 2)
})


test_that("splicegrahm accepts j_incl to include junction coverage plots", {
    ## default plot with simple dataset
    plt <- splicegrahm(simple_cc, j_incl = TRUE)
    bld <- ggplot_build(plt)

    ## check that ggplot is produced
    expect_is(plt, "ggplot")

    ## tabulate the layer geoms
    layers_base <- sapply(plt$layer, function(x) class(x$geom)[1])

    ## check one GeomPaths per junction
    expect_equal(sum(layers_base == "GeomPath"), length(juncs(simple_cc)))
    
    ## check GeomRects for exon/junc heatmap, exon outlines, junc outlines
    expect_equal(sum(layers_base == "GeomRect"), 3)

    ## check two GeomText for junction labels
    expect_equal(sum(layers_base == "GeomText"), 2)

    ## check ranges of splice junctions (flip for negative strand)
    gp_idx <- which(layers_base == "GeomPath")
    gp_starts <- (-1) * sapply(gp_idx, function(idx) max(bld$data[[idx]]$x))
    gp_ends <- (-1) * sapply(gp_idx, function(idx) min(bld$data[[idx]]$x))
    expect_equal(start(juncs(simple_cc)), gp_starts)
    expect_equal(end(juncs(simple_cc)), gp_ends)

    ## check number of rects matchs #exons, #junctions, #heatmap boxes
    gr_idx <- which(layers_base == "GeomRect")
    gr_rows <- sort(sapply(gr_idx, function(idx) nrow(bld$data[[idx]])))
    expect_equal(gr_rows, sort(c(length(exons(simple_cc)), length(juncs(simple_cc)),
                                 prod(dim(exonValues(simple_cc))) +
                                     prod(dim(juncValues(simple_cc))))))

    ## check heatmap boxes placed in correct positions
    exon_frame_color <- "#3C3C3C"
    gr_dat <- do.call(rbind, bld$data[gr_idx])
    gr_dat <- gr_dat[gr_dat$colour == exon_frame_color, ]
    gr_starts <- (-1) * gr_dat$xmin
    gr_ends <- (-1) * gr_dat$xmax
    expect_equal(start(exons(simple_cc)), gr_starts)
    expect_equal(end(exons(simple_cc)), gr_ends)

    ## add tests for junction boxes?
    ## add tests for heatmap placement?
})


test_that("splicegrahm accepts sort_sep to sort exons, juncs separately", {
    ## default plot with simple dataset
    plt <- splicegrahm(simple_cc)
    bld <- ggplot_build(plt)
    
})


test_that("splicegrahm accepts sort_idx to maually specify sort order", {
    ## default plot with simple dataset
    plt <- splicegrahm(simple_cc)
    bld <- ggplot_build(plt)
})

#' #' @param obj a \code{concomp} object
#' #' @param j_incl a logical whether to include heatmaps for junctions
#' #'        (default = FALSE)
#' #' @param sort_sep a logical whether to sort each exon, junction separately
#' #'        (default = FALSE)
#' #' @param sort_idx an integer value specifying the order of the samples in
#' #'        each exon, see details for more information on all possible
#' #'        input, if length is n, then this ordering is used (default = 1)
#' 

test_that("splicegrahm accepts log_base, log_shift to scale heatmap colorscale", {
    ## default plot with simple dataset
    plt <- splicegrahm(simple_cc)
    bld <- ggplot_build(plt)
})


test_that("splicegrahm accepts bin FALSE to plot on continuous colorscale", {
    ## default plot with simple dataset
    plt <- splicegrahm(simple_cc)
    bld <- ggplot_build(plt)
})

#' #' @param log_base a numeric specifying the scale of the binning for the
#' #'        plotting of expression values at each exon, which 0 resulting in no long
#' #'        scaling being applied (default = 10)
#' #' @param log_shift a numeric specifying the shift to be used in the log transformation
#' #'        for handling 0 values (default = 1)
#' #' @param bin a logical whether to bin the values for easier viewing (default = TRUE)

test_that("splicegrahm accepts genomic, ex_use input to adjust x-axis scaling", {
    ## default plot with simple dataset
    plt <- splicegrahm(simple_cc)
    bld <- ggplot_build(plt)
})

#' #' @param genomic a logical whether genomic coordinates should be used to
#' #'        plot the heatmap (default = TRUE)
#' #' @param ex_use a numeric specifying the proportion of the plot exons should occupy if
#' #'        non-genomic coordinate plotting is desired (default = 2/3)

test_that("splicegrahm accepts flip_neg to adjust whether stranded-ness matters", {
    ## default plot with simple dataset
    plt <- splicegrahm(simple_cc)
    bld <- ggplot_build(plt)
})

#' #' @param flip_neg a logical whether to flip plotting of genes on negative strand
#' #'        to run left to right (default = TRUE)

test_that("splicegrahm accepts highlight input to add group annotations", {
    ## default plot with simple dataset
    plt <- splicegrahm(simple_cc)
    bld <- ggplot_build(plt)
})

#' #' @param highlight a vector of labels to highlight samples in groups or clusters
#' #'        (default = NULL)

test_that("splicegrahm accepts use_blk input to invert background color", {
    ## default plot with simple dataset
    plt <- splicegrahm(simple_cc)
    bld <- ggplot_build(plt)
})

#' #' @param use_blk a logical whether to use a black background (default = FALSE)

test_that("splicegrahm accepts eps, txlist, txdb, orgdb input to add gene annotations", {
    ## default plot with simple dataset
    plt <- splicegrahm(simple_cc)
    bld <- ggplot_build(plt)
})

#' #' @param eps a numeric value specifying the number of base pairs around \code{obj} to look
#' #'        for overlapping gene models, if eps = NULL, then all overlapping gene models are
#' #'        included (default = 1e4)
#' #' @param txlist a GRangesList of transcripts or genes which should be queried and
#' #'        added to the plot if falling within the region of the connected component
#' #'        (default = NULL)
#' #' @param txdb a transcript database which can be used to query the transcript IDs
#' #'        identified from txlist (default = NULL)
#' #' @param orgdb a database that can be queried using keys obtained from \code{txdb}
#' #'        to determine corresponding gene symbols (default = NULL)

test_that("splicegrahm accepts title parameter to add title to plot", {
    ## default plot with simple dataset
    plt <- splicegrahm(simple_cc)
    bld <- ggplot_build(plt)
})

#' #' @param title a character string title printed at the top of plot (default = "")
