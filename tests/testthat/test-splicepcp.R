context("splicepcp plotting method")
library("GenomicRanges")
library("ggplot2")

## load annotations packages
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("BSgenome.Hsapiens.UCSC.hg19")
library("org.Hs.eg.db")


## ##############################################################################
## prepare reference hg19/hg37 genome

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
isActiveSeq(txdb)[seqlevels(txdb)] <- FALSE
isActiveSeq(txdb)["chr19"] <- TRUE
exbytx <- exonsBy(txdb, "tx")


## ##############################################################################
## construct example dataset for test cases

## simple gene model overlapping KLK12 on chr19
ichr <- "chr19"
iseqlengths <- 59128983
igStart <- 51532000 + c(0, 2000, 5000, 1000, 1000, 3000)
igStop <- 51532000 + c(1000, 3000, 6000, 2000, 5000, 5000)
ikind <- c(rep("e", 3), rep("j", 3))
istr <- "-"

## samples w/ exon/junc coverage random negative binomial 
samp_cov <- sapply(c(s=0:19),
                   function(x) rnbinom(6, size=1.5, pro=.005) + 1,
                   simplify=FALSE)

## construct concomp
simple_df <- data.frame(chr=ichr, seqlengths=iseqlengths,
                        gStart=igStart, gStop=igStop,
                        kind=ikind,
                        strand=istr,
                        samp_cov)
simple_cc <- concomp(simple_df)

## fix corresponding reference genomes
seqinfo(eRanges(simple_cc)) <- seqinfo(exbytx)[seqlevels(eRanges(simple_cc))]
seqinfo(jRanges(simple_cc)) <- seqinfo(exbytx)[seqlevels(jRanges(simple_cc))]


## ##############################################################################
## actual test cases

test_that("splicepcp default call works", {
    ## default plot with simple dataset
    expect_silent(plt <- splicepcp(simple_cc))

    ## check that returned plot includes 2 ggplots
    expect_is(plt, "Tracks")
    expect_length(plt@plot, 2)
    expect_is(plt@plot[[1]], "ggplot")
    expect_is(plt@plot[[2]], "ggplot")
    
    bld_p1 <- ggplot_build(plt@plot[[1]])
    bld_p2 <- ggplot_build(plt@plot[[2]])

    ## determine layer geoms
    layers_p1 <- sapply(bld_p1$plot$layers, function(x) class(x$geom)[1])
    layers_p2 <- sapply(bld_p2$plot$layers, function(x) class(x$geom)[1])

    ## check that parallel coord plot includes segments for coverage
    expect_equal(sum(layers_p1 == "GeomSegment"), 1)

    ## check that gene model includes one GeomPaths is drawn per junction
    expect_equal(sum(layers_p2 == "GeomPath"), length(jRanges(simple_cc)))

    ## check that gene model includes two GeomRects for exons and outlines
    expect_equal(sum(layers_p2 == "GeomRect"), 2)
})


test_that("splicepcp accepts log_base, log_shift to scale y-axis plotting", {
    ## default plot with simple dataset, default log_base=10, log_shift=1
    plt_base <- splicepcp(simple_cc)

    ## plot with log scaling and shifting
    expect_silent(plt_b2s1 <- splicepcp(simple_cc, log_base=2, log_shift=1))
    expect_silent(plt_b2s0 <- splicepcp(simple_cc, log_base=2, log_shift=0))
    ## plot without log scaling (log_base=0)
    expect_silent(plt_b0s1 <- splicepcp(simple_cc, log_base=0, log_shift=1))
    expect_silent(plt_b0s0 <- splicepcp(simple_cc, log_base=0, log_shift=0))

    ## check that Tracks plot is produced
    expect_is(plt_b2s1, "Tracks")
    expect_is(plt_b2s0, "Tracks")
    expect_is(plt_b0s1, "Tracks")
    expect_is(plt_b0s0, "Tracks")

    ## check that all gene models (track 2) are the equivalent
    expect_equivalent(plt_base@grobs[[2]], plt_b2s1@grobs[[2]])
    expect_equivalent(plt_base@grobs[[2]], plt_b2s0@grobs[[2]])
    expect_equivalent(plt_base@grobs[[2]], plt_b0s1@grobs[[2]])
    expect_equivalent(plt_base@grobs[[2]], plt_b0s0@grobs[[2]])

    ## check that y-axis labels changed in pcp (track 1)
    expect_equal(plt_base@plot[[1]]$labels$y, "log10 expression")
    expect_equal(plt_b2s1@plot[[1]]$labels$y, "log2 expression")
    expect_equal(plt_b2s0@plot[[1]]$labels$y, "log2 expression")
    expect_equal(plt_b0s1@plot[[1]]$labels$y, "expression")
    expect_equal(plt_b0s0@plot[[1]]$labels$y, "expression")

    ## check that plotted data is actually changed in pcp (track 1)
    bld_base_p1 <- ggplot_build(plt_base@plot[[1]])
    bld_b2s1_p1 <- ggplot_build(plt_b2s1@plot[[1]])
    bld_b2s0_p1 <- ggplot_build(plt_b2s0@plot[[1]])
    bld_b0s1_p1 <- ggplot_build(plt_b0s1@plot[[1]])
    bld_b0s0_p1 <- ggplot_build(plt_b0s0@plot[[1]])

    ## check that log_shift not used when log_base = 0 
    expect_equal(bld_b0s1_p1$data[[2]]$y, bld_b0s0_p1$data[[2]]$y)

    ## check that log_base, log_shift transform log_base = 0 (raw) values
    expect_equal(bld_base_p1$data[[2]]$y, log10(1 + bld_b0s0_p1$data[[2]]$y))
    expect_equal(bld_b2s1_p1$data[[2]]$y, log2(1 + bld_b0s0_p1$data[[2]]$y))
    expect_equal(bld_b2s0_p1$data[[2]]$y, log2(0 + bld_b0s0_p1$data[[2]]$y))
})


test_that("splicepcp accepts genomic, ex_use input to adjust x-axis scaling", {
    ## plot with simple dataset
    expect_silent(plt_gen <- splicepcp(simple_cc, genomic = TRUE))
    expect_silent(plt_non <- splicepcp(simple_cc, genomic = FALSE))
    expect_silent(plt_non_10 <- splicepcp(simple_cc, genomic = FALSE, ex_use = 1.0))

    ## plot with ex_use less than exon proportion in genomic space
    expect_warning(plt_non_01 <- splicepcp(simple_cc, genomic = FALSE, ex_use = 0.1),
                   paste0("Using 'genomic = TRUE' since exons occupy more than specified 'ex_use' ",
                          "proportion of plot in genomic coordinates. No need to squish."))

    ## check that Tracks plot is produced
    expect_is(plt_gen, "Tracks")
    expect_is(plt_non, "Tracks")
    expect_is(plt_non_10, "Tracks")
    expect_is(plt_non_01, "Tracks")
    bld_gen_p1 <- ggplot_build(plt_gen@plot[[1]])
    bld_gen_p2 <- ggplot_build(plt_gen@plot[[2]])
    bld_non_p1 <- ggplot_build(plt_non@plot[[1]])
    bld_non_p2 <- ggplot_build(plt_non@plot[[2]])
    bld_non_10_p1 <- ggplot_build(plt_non_10@plot[[1]])
    bld_non_10_p2 <- ggplot_build(plt_non_10@plot[[2]])
    bld_non_01_p1 <- ggplot_build(plt_non_01@plot[[1]])
    bld_non_01_p2 <- ggplot_build(plt_non_01@plot[[2]])

    ## check that coordinate range shrinks in expected proportions
    width_gen_p1 <- diff(bld_gen_p1$layout$panel_scales_x[[1]]$range$range)
    width_non_p1 <- diff(bld_non_p1$layout$panel_scales_x[[1]]$range$range)
    width_non_10_p1 <- diff(bld_non_10_p1$layout$panel_scales_x[[1]]$range$range)
    expect_lt(width_non_p1, width_gen_p1)
    expect_lt(width_non_10_p1, width_gen_p1)
    expect_equal(width_non_10_p1 / width_non_p1, 2/3, tolerance=0.001)
    
    ## check that bad (small) ex_use creates same plot as genomic = TRUE
    expect_equal(bld_gen_p1$data, bld_non_01_p1$data)
    expect_equal(bld_gen_p1$panel, bld_non_01_p1$panel)
    expect_equal(bld_gen_p1$plot[names(bld_gen_p1$plot) != "plot_env"],
                 bld_non_01_p1$plot[names(bld_non_01_p1$plot) != "plot_env"])

    ## check that coordinates are same for both tracks
    width_gen_p2 <- diff(bld_gen_p2$layout$panel_scales_x[[1]]$range$range)
    width_non_p2 <- diff(bld_non_p2$layout$panel_scales_x[[1]]$range$range)
    expect_equal(width_gen_p1, width_gen_p2)
    expect_equal(width_non_p1, width_non_p2)
    
    ## add tests for when txlist specified
})


test_that("splicepcp accepts flip_neg to adjust whether stranded-ness matters", {
    ## only including basic tests - flip_neg should be refactored out
    expect_silent(plt_neg_noflip <- splicepcp(simple_cc, flip_neg=FALSE))
    expect_silent(plt_neg_flip <- splicepcp(simple_cc, flip_neg=TRUE))
    expect_silent(plt_neg_annot <- splicepcp(simple_cc, flip_neg=FALSE, txlist=exbytx))

    expect_is(plt_neg_noflip, "Tracks")
    expect_is(plt_neg_flip, "Tracks")
    expect_is(plt_neg_annot, "Tracks")
    bld_neg_noflip <- ggplot_build(plt_neg_noflip@plot[[1]])
    bld_neg_flip <- ggplot_build(plt_neg_flip@plot[[1]])
    bld_neg_annot <- ggplot_build(plt_neg_annot@plot[[1]])

    ## check that coord flipped w/ neg strand
    expect_equal(bld_neg_noflip$layout$panel_scales_x[[1]]$range$range,
                 (-1) * rev(bld_neg_flip$layout$panel_scales_x[[1]]$range$range))
    
    ## check that coord not flipped w/ neg strand and neg annot if FALSE
    expect_equal(bld_neg_noflip$layout$panel_scales_x[[1]]$range$range,
                 bld_neg_annot$layout$panel_scales_x[[1]]$range$range)
})


test_that("splicepcp accepts highlight input to color plot by groups", {
    ## check that string labels
    hl_s1 <- rep("A", each=20)
    hl_s2 <- rep(c("A", "B"), each=10)
    hl_s10 <- rep(letters[1:10], 2)
    hl_f2 <- factor(hl_s2)
    hl_i2 <- rep(1:2, each=10)
    hl_short <- rep("A", 10)
    hl_long <- rep("A", 35)

    ## simple plot w/ no highlighting
    expect_silent(plt_base <- splicepcp(simple_cc))

    ## check that any number of highlight labels are accepted
    expect_silent(plt_s1 <- splicepcp(simple_cc, highlight=hl_s1))
    expect_silent(plt_s2 <- splicepcp(simple_cc, highlight=hl_s2))

    ## check that too many levels returns an error
    expect_error(plt_s10 <- splicepcp(simple_cc, highlight=hl_s10),
                 "highlight currently only support up to 9 unique groups.")

    ## check that factors can be specified
    expect_silent(plt_f2 <- splicepcp(simple_cc, highlight=hl_f2))

    ## check that integers can be specified
    expect_silent(plt_i2 <- splicepcp(simple_cc, highlight=hl_i2))
    
    ## chech that incorrect length vector returns error
    expect_error(plt_short <- splicepcp(simple_cc, highlight=hl_short),
                 "highlight must be a vector of length n.")
    expect_error(plt_long <- splicepcp(simple_cc, highlight=hl_long),
                 "highlight must be a vector of length n.")

    ## check that specified labels are used
    bld_s2_p1 <- ggplot_build(plt_s2@plot[[1]])
    expect_equivalent(unique(bld_s2_p1$plot$data[, c("variable", "hl")]),
                      data.frame(variable=paste0("s", 1:20), hl=hl_s2))
    expect_equal(length(unique(bld_s2_p1$data[[2]]$colour)),
                 length(unique(hl_s2)))
})


test_that("splicepcp only includes concomp gene model when imodel=TRUE", {
    ## default plot with simple dataset
    plt_base <- splicepcp(simple_cc)
    bld_base_p1 <- ggplot_build(plt_base@plot[[1]])
    bld_base_p2 <- ggplot_build(plt_base@plot[[2]])

    ## check that plot can be created w/ imodel = FALSE
    expect_silent(plt_nomodel <- splicepcp(simple_cc, imodel=FALSE))

    ## check that plot is a single ggplot
    expect_is(plt_nomodel, "ggplot")
    bld_nomodel <- ggplot_build(plt_nomodel)
    
    ## check that plot w/ imodel=FALSE matches pcp plot (panel 1) in default plot
    ## - for some reason expect_equal on full bld_* obj fails to run
    expect_equal(bld_nomodel$data, bld_base_p1$data)
    expect_equal(bld_nomodel$panel, bld_base_p1$panel)
    expect_equal(bld_nomodel$plot[names(bld_nomodel$plot) != "plot_env"],
                 bld_base_p1$plot[names(bld_base_p1$plot) != "plot_env"])

    ## check that Tracks plot still created if txlist specified
    expect_silent(plt_nom_annot <- splicepcp(simple_cc, imodel=FALSE, txlist=exbytx))
    expect_is(plt_nom_annot, "Tracks")
})


test_that("splicepcp accepts txlist, txdb, orgdb input to add gene annotations", {
    ## only including basic tests - parsing tested more w/ splicegrahm cases
    plt_base <- splicepcp(simple_cc)
    bld_base_p1 <- ggplot_build(plt_base@plot[[1]])
    bld_base_p2 <- ggplot_build(plt_base@plot[[2]])

    ## create a plot with a transcript annotations track
    expect_silent(plt_tx <- splicepcp(simple_cc, txlist=exbytx,
                                      txdb=txdb, orgdb=org.Hs.eg.db))
    expect_is(plt_tx, "Tracks")
    expect_true(length(plt_tx@plot) == 3)
    expect_is(plt_tx@plot[[1]], "ggplot")
    expect_is(plt_tx@plot[[2]], "ggplot")
    expect_is(plt_tx@plot[[3]], "ggplot")
    
    bld_tx_p1 <- ggplot_build(plt_tx@plot[[1]])
    bld_tx_p2 <- ggplot_build(plt_tx@plot[[2]])
    bld_tx_p3 <- ggplot_build(plt_tx@plot[[3]])
    
    ## check that top 2 tracks are same as basic plot (at least data and layers)
    expect_equivalent(bld_tx_p1$data, bld_base_p1$data)
    expect_equivalent(bld_tx_p1$plot$layers, bld_base_p1$plot$layers)
    expect_equivalent(bld_tx_p2$data, bld_base_p2$data)
    expect_equivalent(bld_tx_p2$plot$layers, bld_base_p2$plot$layers)

    ## create a plot with only transcripts fully contained in concomp range (eps=0)
    expect_silent(plt_tx_e0 <- splicepcp(simple_cc, txlist=exbytx, eps=0))
    expect_is(plt_tx_e0, "Tracks")
    bld_tx_e0_p3 <- ggplot_build(plt_tx_e0@plot[[3]]) 

    ## check that eps=0 only keeps fully contained transcript models (just one here)
    expect_equal(bld_tx_e0_p3$layout$panel_scales_y[[1]]$labels, "70043")
})
