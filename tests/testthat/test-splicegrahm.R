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
    expect_silent(plt <- splicegrahm(simple_cc))
    
    ## check that ggplot is produced
    expect_is(plt, "ggplot")
    bld <- ggplot_build(plt)

    ## check ranges of plot 
    expect_equal((-1) * bld$panel$ranges[[1]]$x.range[2], 51028700)
    expect_equal((-1) * bld$panel$ranges[[1]]$x.range[1], 51035300)

    ## determine layer geoms
    layers_base <- sapply(plt$layer, function(x) class(x$geom)[1])

    ## check one GeomPaths per junction
    expect_equal(sum(layers_base == "GeomPath"), length(juncs(simple_cc)))
    
    ## check GeomRects for exon heatmap, exon outlines
    expect_equal(sum(layers_base == "GeomRect"), 2)
})


test_that("splicegrahm accepts j_incl to include junction coverage plots", {
    ## plot with simple dataset
    expect_silent(plt <- splicegrahm(simple_cc, j_incl = TRUE))

    ## check that ggplot is produced
    expect_is(plt, "ggplot")
    bld <- ggplot_build(plt)

    ## determine layer geoms
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
    ## plots with simple dataset
    expect_silent(plt_sep <- splicegrahm(simple_cc, j_incl = TRUE, sort_sep = TRUE))
    expect_silent(plt_same <- splicegrahm(simple_cc, j_incl = TRUE, sort_sep = FALSE))

    ## check that ggplot is produced
    expect_is(plt_sep, "ggplot")
    expect_is(plt_same, "ggplot")
    bld_sep <- ggplot_build(plt_same)
    bld_same <- ggplot_build(plt_same)

    ## determine layer geoms
    layers_sep <- sapply(plt_sep$layer, function(x) class(x$geom)[1])
    layers_same <- sapply(plt_same$layer, function(x) class(x$geom)[1])

    ## check that layers are still the same
    expect_equal(layers_sep, layers_same)
    
    ## check that order of heatmaps is reversed for junctions
    dat_e_same <- plt_same$data[plt_same$data$kind == "e", ]
    dat_e_sep <- plt_sep$data[plt_sep$data$kind == "e", ]
    dat_j_same <- plt_same$data[plt_same$data$kind == "j", ]
    dat_j_sep <- plt_sep$data[plt_sep$data$kind == "j", ]
    expect_equivalent(dat_e_same, dat_e_sep)
    expect_equivalent(dat_j_same[, c("ymin", "ymax")],
                      dat_j_sep[, c("ymin", "ymax")])
    expect_equivalent(dat_j_same[, c("value", "variable")],
                      dat_j_sep[nrow(dat_j_sep):1, c("value", "variable")])
})


test_that("splicegrahm accepts sort_idx to maually specify sort order", {
    ## plot with simple dataset
    sidx <- c(11:20, 1:10)
    expect_silent(plt <- splicegrahm(simple_cc, j_incl = TRUE, sort_idx = sidx))

    ## check that ggplot is produced
    expect_is(plt, "ggplot")
    bld <- ggplot_build(plt)

    ## determine layer geoms
    layers_base <- sapply(plt$layer, function(x) class(x$geom)[1])

    ## check that specified n-vector of order is used
    dat_e <- unique(plt$data[plt$data$kind == "e", c("variable", "ymin")])
    dat_j <- unique(plt$data[plt$data$kind == "j", c("variable", "ymin")])
    dat_e <- as.character(dat_e$variable[order(dat_e$ymin)])
    dat_j <- as.character(dat_j$variable[order(dat_j$ymin)])
    expect_equal(dat_e, paste0("s", sidx))
    expect_equal(dat_j, paste0("s", sidx))
})


test_that("splicegrahm accepts log_base, log_shift to scale heatmap colorscale", {
    ## plot with simple dataset
    expect_silent(plt_b2s1 <- splicegrahm(simple_cc, j_incl=TRUE, log_base=2, log_shift=1))
    expect_silent(plt_b2s0 <- splicegrahm(simple_cc, j_incl=TRUE, log_base=2, log_shift=0))
    expect_silent(plt_b10s1 <- splicegrahm(simple_cc, j_incl=TRUE, log_base=10, log_shift=1))
    expect_silent(plt_b10s0 <- splicegrahm(simple_cc, j_incl=TRUE, log_base=10, log_shift=0))
    ## plot without log scaling (log_base=0)
    expect_silent(plt_b0s1 <- splicegrahm(simple_cc, j_incl=TRUE, log_base=0, log_shift=1))
    expect_silent(plt_b0s0 <- splicegrahm(simple_cc, j_incl=TRUE, log_base=0, log_shift=0))

    ## check that ggplot is produced
    expect_is(plt_b2s1, "ggplot")
    expect_is(plt_b2s0, "ggplot")
    expect_is(plt_b10s1, "ggplot")
    expect_is(plt_b10s0, "ggplot")
    expect_is(plt_b0s1, "ggplot")
    expect_is(plt_b0s0, "ggplot")

    ## check that all scales are still discrete
    expect_is(plt_b2s0$scales$scales[[4]], "ScaleDiscrete")
    expect_is(plt_b2s1$scales$scales[[4]], "ScaleDiscrete")
    expect_is(plt_b10s0$scales$scales[[4]], "ScaleDiscrete")
    expect_is(plt_b10s1$scales$scales[[4]], "ScaleDiscrete")
    ## check that scales are changed to continuous if non-log plotting
    expect_is(plt_b0s0$scales$scales[[4]], "ScaleContinuous")
    expect_is(plt_b0s1$scales$scales[[4]], "ScaleContinuous")

    ## check that non-log scaling ignores shift param
    expect_equal(plt_b0s0$data, plt_b0s1$data)
    
    ## check that scales are changed
    expect_equal(plt_b2s1$scales$scales[[4]]$labels,
                 paste0("<", 2^(1:11))) 
    expect_equal(plt_b2s0$scales$scales[[4]]$labels,
                 paste0("<", 2^(1:11))) 
    expect_equal(plt_b10s1$scales$scales[[4]]$labels,
                 paste0("<", 10^(1:4))) 
    expect_equal(plt_b10s0$scales$scales[[4]]$labels,
                 paste0("<", 10^(1:4))) 
})


test_that("splicegrahm accepts bin FALSE to plot on continuous colorscale", {
    ## default plot with simple dataset
    expect_silent(plt_nobin <- splicegrahm(simple_cc, j_incl= TRUE, bin = FALSE, log_base=2))
    expect_silent(plt_wbin <- splicegrahm(simple_cc, j_incl= TRUE, bin = TRUE, log_base=2))

    ## check that ggplots are produced
    expect_is(plt_nobin, "ggplot")
    expect_is(plt_wbin, "ggplot")
    bld_nobin <- ggplot_build(plt_nobin)
    bld_wbin <- ggplot_build(plt_wbin)

    ## check that color scales is discrete only if binning
    expect_is(plt_nobin$scales$scales[[4]], "ScaleContinuous")
    expect_is(plt_wbin$scales$scales[[4]], "ScaleDiscrete")
    
    ## check that value for color is factor only if binning
    expect_false(is.factor(plt_nobin$data$value))
    expect_true(is.factor(plt_wbin$data$value))
})


test_that("splicegrahm accepts genomic, ex_use input to adjust x-axis scaling", {
    ## plot with simple dataset
    expect_silent(plt_gen <- splicegrahm(simple_cc, genomic = TRUE))
    expect_silent(plt_non <- splicegrahm(simple_cc, genomic = FALSE))
    expect_silent(plt_non_10 <- splicegrahm(simple_cc, genomic = FALSE, ex_use = 1.0))

    ## plot with ex_use less than exon proportion in genomic space
    expect_warning(plt_non_01 <- splicegrahm(simple_cc, genomic = FALSE, ex_use = 0.1),
                   paste0("Using 'genomic = TRUE' since exons occupy more than specified 'ex_use' ",
                          "proportion of plot in genomic coordinates. No need to squish."))

    ## check that ggplot is produced
    expect_is(plt_gen, "ggplot")
    expect_is(plt_non, "ggplot")
    expect_is(plt_non_10, "ggplot")
    expect_is(plt_non_01, "ggplot")
    bld_gen <- ggplot_build(plt_gen)
    bld_non <- ggplot_build(plt_non)
    bld_non_10 <- ggplot_build(plt_non_10)
    bld_non_01 <- ggplot_build(plt_non_01)

    ## check that coordinate range shrinks in expected proportions
    width_gen <- diff(bld_gen$panel$ranges[[1]]$x.range)
    width_non <- diff(bld_non$panel$ranges[[1]]$x.range)
    width_non_10 <- diff(bld_non_10$panel$ranges[[1]]$x.range)
    expect_lt(width_non, width_gen)
    expect_lt(width_non_10, width_gen)
    expect_equal(width_non_10 / width_non, 2/3, tolerance=0.001)
    
    ## check that bad (small) ex_use creates same plot as genomic = TRUE
    expect_equal(bld_gen$data, bld_non_01$data)
    expect_equal(bld_gen$panel, bld_non_01$panel)
    expect_equal(plt_gen[names(plt_gen) != "plot_env"],
                 plt_non_01[names(plt_non_01) != "plot_env"])

    ## add tests for when txlist specified
})


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
