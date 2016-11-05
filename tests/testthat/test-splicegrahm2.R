context("splicegrahm2 plotting method")
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
    
## samples w/ exon/junc coverage rank reversed
samp_cov1 <- sapply(c(s=0:19),
                    function(x) c(rep(2^((x+1)/2), 3), rep(2^((20-x)/2), 3)),
                    simplify=FALSE)
samp_cov2 <- sapply(c(s=0:19),
                    function(x) rnbinom(6, size=1.2, pro=.005) + 1,
                    simplify=FALSE)

## construct concomps
simple_df1 <- data.frame(chr=ichr, seqlengths=iseqlengths,
                         gStart=igStart, gStop=igStop,
                         kind=ikind, strand=istr, samp_cov1)
simple_cc1 <- concomp(simple_df1)
simple_df2 <- data.frame(chr=ichr, seqlengths=iseqlengths,
                         gStart=igStart, gStop=igStop,
                         kind=ikind, strand=istr, samp_cov2)
simple_cc2 <- concomp(simple_df2)

## construct smaller (less samples) concomp
simple_df3 <- data.frame(chr=ichr, seqlengths=iseqlengths,
                         gStart=igStart, gStop=igStop,
                         kind=ikind, strand=istr, samp_cov2[1:10])
simple_cc3 <- concomp(simple_df3)

## fix corresponding reference genomes
seqinfo(eRanges(simple_cc1)) <- seqinfo(exbytx)[seqlevels(eRanges(simple_cc1))]
seqinfo(jRanges(simple_cc1)) <- seqinfo(exbytx)[seqlevels(jRanges(simple_cc1))]

seqinfo(eRanges(simple_cc2)) <- seqinfo(exbytx)[seqlevels(eRanges(simple_cc2))]
seqinfo(jRanges(simple_cc2)) <- seqinfo(exbytx)[seqlevels(jRanges(simple_cc2))]

seqinfo(eRanges(simple_cc3)) <- seqinfo(exbytx)[seqlevels(eRanges(simple_cc3))]
seqinfo(jRanges(simple_cc3)) <- seqinfo(exbytx)[seqlevels(jRanges(simple_cc3))]


## ##############################################################################
## actual test cases

test_that("splicegrahm2 default call works", {
    ## default plot with simple dataset
    expect_silent(plt <- splicegrahm2(simple_cc1, simple_cc1))

    ## default splicegrahm plot
    plt_1 <- splicegrahm(simple_cc1)

    ## generate plot
    bld_p1 <- ggplot_build(plt@plot[[1]])
    bld_p2 <- ggplot_build(plt@plot[[2]])
    bld_1 <- ggplot_build(plt_1)

    layers_p1 <- sapply(bld_p1$plot$layer, function(x) class(x$geom)[1])
    layers_p2 <- sapply(bld_p2$plot$layer, function(x) class(x$geom)[1])
    layers_1 <- sapply(bld_1$plot$layer, function(x) class(x$geom)[1])
    
    ## check that top panel is standard splicegrahm
    expect_equal(bld_p1$data, bld_1$data)
    expect_equal(layers_p1, layers_1)
    
    ## check that bottom panel is same geoms but flipped
    expect_equal(layers_p2, layers_1)
    ## check exon frame layers are flipped
    expect_equal(c(bld_p1$data[[2]]$ymin, bld_p1$data[[2]]$ymax),
                 (-1) * c(bld_p2$data[[2]]$ymax, bld_p2$data[[2]]$ymin))
    ## check that heatmap layers are flipped
    expect_equal(c(bld_p1$data[[1]]$ymin, bld_p1$data[[1]]$ymax),
                 (-1) * c(bld_p2$data[[1]]$ymax, bld_p2$data[[1]]$ymin))
    ## check that arcs are flipped
    expect_equal(do.call(rbind, bld_p1$data[which(layers_p1 == "GeomPath")])$y,
                 (-1) * do.call(rbind, bld_p2$data[which(layers_p2 == "GeomPath")])$y)
})


test_that("splicegrahm2 accepts j_incl parameter", {
    ## default plot with simple dataset
    expect_silent(plt <- splicegrahm2(simple_cc1, simple_cc2, j_incl=TRUE))

    ## default splicegrahm plot
    plt_1 <- splicegrahm(simple_cc1, j_incl=TRUE)
    plt_2 <- splicegrahm(simple_cc2, j_incl=TRUE)

    ## generate plot
    bld_p1 <- ggplot_build(plt@plot[[1]])
    bld_p2 <- ggplot_build(plt@plot[[2]])
    bld_1 <- ggplot_build(plt_1)
    bld_2 <- ggplot_build(plt_2)

    layers_p1 <- sapply(bld_p1$plot$layer, function(x) class(x$geom)[1])
    layers_p2 <- sapply(bld_p2$plot$layer, function(x) class(x$geom)[1])
    layers_1 <- sapply(bld_1$plot$layer, function(x) class(x$geom)[1])
    layers_2 <- sapply(bld_2$plot$layer, function(x) class(x$geom)[1])
    
    ## check that top panel is splicegrahm w/ j_incl=TRUE
    expect_equal(bld_p1$data, bld_1$data)
    expect_equal(layers_p1, layers_1)
    
    ## check that bottom panel is flipped splicegrahm w/ j_incl=TRUE
    expect_equal(layers_p2, layers_2)
    ## check that arcs are flipped
    expect_equal(do.call(rbind, bld_2$data[which(layers_2 == "GeomPath")])$y,
                 (-1) * do.call(rbind, bld_p2$data[which(layers_p2 == "GeomPath")])$y)
    ## check that all rects, heatmap, frames are flipped
    expect_equal(do.call(rbind, bld_2$data[which(layers_2 == "GeomRect")])$ymin,
                 (-1) * do.call(rbind, bld_p2$data[which(layers_p2 == "GeomRect")])$ymax)
    expect_equal(do.call(rbind, bld_2$data[which(layers_2 == "GeomRect")])$ymax,
                 (-1) * do.call(rbind, bld_p2$data[which(layers_p2 == "GeomRect")])$ymin)
    ## check that all labels  are flipped
    expect_equal(do.call(rbind, bld_2$data[which(layers_2 == "GeomText")])$y,
                 (-1) * do.call(rbind, bld_p2$data[which(layers_p2 == "GeomText")])$y)    
})


test_that("splicegrahm2 accepts sort_sep parameter", {
    ## default splicegrahm2 plot with sort_sep = TRUE (default is sort_sep = FALSE)
    expect_silent(plt_sep <- splicegrahm2(simple_cc1, simple_cc2, j_incl=TRUE, sort_sep=TRUE))
    ## expect_silent(plt_same <- splicegrahm2(simple_cc1, simple_cc2, j_incl=TRUE, sort_sep=FALSE))

    ## create splicegrahm plot with sort_sep = TRUE (default is sort_sep = FALSE)
    plt_sep_1 <- splicegrahm(simple_cc1, j_incl=TRUE, sort_sep=TRUE)
    plt_sep_2 <- splicegrahm(simple_cc2, j_incl=TRUE, sort_sep=TRUE)
    ## plt_same_1 <- splicegrahm(simple_cc1, j_incl=TRUE, sort_sep=FALSE)
    ## plt_same_2 <- splicegrahm(simple_cc2, j_incl=TRUE, sort_sep=FALSE)

    ## generate plot summaries
    bld_sep_p1 <- ggplot_build(plt_sep@plot[[1]])
    bld_sep_p2 <- ggplot_build(plt_sep@plot[[2]])
    bld_sep_1 <- ggplot_build(plt_sep_1)
    bld_sep_2 <- ggplot_build(plt_sep_2)
    ## bld_same_p1 <- ggplot_build(plt_same@plot[[1]])
    ## bld_same_p2 <- ggplot_build(plt_same@plot[[2]])
    ## bld_same_1 <- ggplot_build(plt_same_1)
    ## bld_same_2 <- ggplot_build(plt_same_2)

    layers_sep_p1 <- sapply(bld_sep_p1$plot$layer, function(x) class(x$geom)[1])
    layers_sep_p2 <- sapply(bld_sep_p2$plot$layer, function(x) class(x$geom)[1])
    layers_sep_1 <- sapply(bld_sep_1$plot$layer, function(x) class(x$geom)[1])
    layers_sep_2 <- sapply(bld_sep_2$plot$layer, function(x) class(x$geom)[1])
    ## layers_same_p1 <- sapply(bld_same_p1$plot$layer, function(x) class(x$geom)[1])
    ## layers_same_p2 <- sapply(bld_same_p2$plot$layer, function(x) class(x$geom)[1])
    ## layers_same_1 <- sapply(bld_same_1$plot$layer, function(x) class(x$geom)[1])
    ## layers_same_2 <- sapply(bld_same_2$plot$layer, function(x) class(x$geom)[1])

    ## check that top panel is splicegrahm w/ sort_sep=TRUE
    expect_equal(bld_sep_p1$data, bld_sep_1$data)
    expect_equal(layers_sep_p1, layers_sep_1)
    
    ## check that bottom panel is flipped splicegrahm w/ sort_sep=TRUE
    expect_equal(layers_sep_p2, layers_sep_2)
    ## check that arcs are flipped
    expect_equal(do.call(rbind, bld_sep_2$data[which(layers_sep_2 == "GeomPath")])$y,
                 (-1) * do.call(rbind, bld_sep_p2$data[which(layers_sep_p2 == "GeomPath")])$y)
    ## check that all rects, heatmap, frames are flipped
    expect_equal(do.call(rbind, bld_sep_2$data[which(layers_sep_2 == "GeomRect")])$ymin,
                 (-1) * do.call(rbind, bld_sep_p2$data[which(layers_sep_p2 == "GeomRect")])$ymax)
    expect_equal(do.call(rbind, bld_sep_2$data[which(layers_sep_2 == "GeomRect")])$ymax,
                 (-1) * do.call(rbind, bld_sep_p2$data[which(layers_sep_p2 == "GeomRect")])$ymin)
    ## check that all labels  are flipped
    expect_equal(do.call(rbind, bld_sep_2$data[which(layers_sep_2 == "GeomText")])$y,
                 (-1) * do.call(rbind, bld_sep_p2$data[which(layers_sep_p2 == "GeomText")])$y)    
})


test_that("splicegrahm2 accepts sort_idx to maually specify sort order", {
    ## plot splicegrahm2 with sort_idx specified
    sidx1 <- c(11:20, 1:10)
    sidx2 <- c(16:20, 15:1)
    expect_silent(plt <- splicegrahm2(simple_cc1, simple_cc2, j_incl=TRUE,
                                      sort_idx1=sidx1, sort_idx2=sidx2))

    plt_1 <- splicegrahm(simple_cc1, j_incl=TRUE, sort_idx = sidx1)
    plt_2 <- splicegrahm(simple_cc2, j_incl=TRUE, sort_idx = sidx2)

    ## generate plot summaries
    bld_p1 <- ggplot_build(plt@plot[[1]])
    bld_p2 <- ggplot_build(plt@plot[[2]])
    bld_1 <- ggplot_build(plt_1)
    bld_2 <- ggplot_build(plt_2)
    
    ## check that data used for plotting is same in splicegrahm2, splicegrahm plots 
    expect_equivalent(bld_p1$plot$data, bld_1$plot$data)
    expect_equivalent(bld_p2$plot$data[, -which(grepl("^y", names(bld_p2$plot$data)))],
                      bld_2$plot$data[, -which(grepl("^y", names(bld_2$plot$data)))])
})


test_that("splicegrahm2 accepts log_base, log_shift to scale heatmap colorscale", {
    ## plot splicegrahm2 with log_base, log_shift parameters 
    expect_silent(plt_b2s0 <- splicegrahm2(simple_cc1, simple_cc2, j_incl=TRUE,
                                           log_base = 2, log_shift = 0))
    expect_silent(plt_b0s1 <- splicegrahm2(simple_cc1, simple_cc2, j_incl=TRUE,
                                           log_base = 0, log_shift = 1))

    plt_b2s0_1 <- splicegrahm(simple_cc1, j_incl=TRUE, log_base = 2, log_shift = 0)
    plt_b2s0_2 <- splicegrahm(simple_cc2, j_incl=TRUE, log_base = 2, log_shift = 0)
    plt_b0s1_1 <- splicegrahm(simple_cc1, j_incl=TRUE, log_base = 0, log_shift = 1)
    plt_b0s1_2 <- splicegrahm(simple_cc2, j_incl=TRUE, log_base = 0, log_shift = 1)

    ## generate plot summaries
    bld_b2s0_p1 <- ggplot_build(plt_b2s0@plot[[1]])
    bld_b2s0_p2 <- ggplot_build(plt_b2s0@plot[[2]])
    bld_b2s0_1 <- ggplot_build(plt_b2s0_1)
    bld_b2s0_2 <- ggplot_build(plt_b2s0_2)
    bld_b0s1_p1 <- ggplot_build(plt_b0s1@plot[[1]])
    bld_b0s1_p2 <- ggplot_build(plt_b0s1@plot[[2]])
    bld_b0s1_1 <- ggplot_build(plt_b0s1_1)
    bld_b0s1_2 <- ggplot_build(plt_b0s1_2)
    
    ## check that data used for plotting is same in splicegrahm2, splicegrahm plots 
    expect_equivalent(bld_b2s0_p1$plot$data, bld_b2s0_1$plot$data)
    expect_equivalent(bld_b2s0_p2$plot$data[, -which(grepl("^y", names(bld_b2s0_p2$plot$data)))],
                      bld_b2s0_2$plot$data[, -which(grepl("^y", names(bld_b2s0_2$plot$data)))])
    
    ## check that data used for plotting is same in splicegrahm2, splicegrahm plots
    ## - fails with expect_equal
    expect_equivalent(bld_b0s1_p1$plot$data, bld_b0s1_1$plot$data)
    expect_equivalent(bld_b0s1_p2$plot$data[, -which(grepl("^y", names(bld_b0s1_p2$plot$data)))],
                      bld_b0s1_2$plot$data[, -which(grepl("^y", names(bld_b0s1_2$plot$data)))])
})


test_that("splicegrahm2 accepts bin FALSE to plot on continuous colorscale", {
    ## plot splicegrahm2 with bin = FALSE (default is bin = TRUE)
    expect_silent(plt_nobin <- splicegrahm2(simple_cc1, simple_cc2, j_incl=TRUE, bin=FALSE))
    plt_nobin_p1 <- plt_nobin@plot[[1]]
    plt_nobin_p2 <- plt_nobin@plot[[2]]
    
    ## check that color scales is continuous when not binning
    expect_is(plt_nobin_p1$scales$scales[[4]], "ScaleContinuous")
    expect_is(plt_nobin_p2$scales$scales[[4]], "ScaleContinuous")
    
    ## check that value for color is not factor when not binning
    expect_false(is.factor(plt_nobin_p1$data$value))
    expect_false(is.factor(plt_nobin_p2$data$value))
})


test_that("splicegrahm2 accepts genomic, ex_use input to adjust x-axis scaling", {
    ## plot with simple dataset
    expect_silent(plt_gen <- splicegrahm2(simple_cc1, simple_cc2, genomic = TRUE))
    expect_silent(plt_non <- splicegrahm2(simple_cc1, simple_cc2, genomic = FALSE))
    expect_silent(plt_non_10 <- splicegrahm2(simple_cc1, simple_cc2, genomic = FALSE, ex_use = 1.0))

    ## plot with ex_use less than exon proportion in genomic space
    expect_warning(plt_non_01 <- splicegrahm2(simple_cc1, simple_cc2, genomic = FALSE, ex_use = 0.1),
                   paste0("Using 'genomic = TRUE' since exons occupy more than specified 'ex_use' ",
                          "proportion of plot in genomic coordinates. No need to squish."))

    bld_gen_p1 <- ggplot_build(plt_gen@plot[[1]])
    bld_non_p1 <- ggplot_build(plt_non@plot[[1]])
    bld_non_10_p1 <- ggplot_build(plt_non_10@plot[[1]])
    
    ## check that coordinate range shrinks in expected proportions
    width_gen <- diff(bld_gen_p1$panel$ranges[[1]]$x.range)
    width_non <- diff(bld_non_p1$panel$ranges[[1]]$x.range)
    width_non_10 <- diff(bld_non_10_p1$panel$ranges[[1]]$x.range)
    expect_lt(width_non, width_gen)
    expect_lt(width_non_10, width_gen)
    expect_equal(width_non_10 / width_non, 2/3, tolerance=0.001)

    ## add tests with concomp1, concomp2 having different exon definitions
})


test_that("splicegrahm2 accepts flip_neg to adjust whether stranded-ness matters", {
    ## only including basic tests - flip_neg should be refactored out
    expect_silent(plt_neg_noflip <- splicegrahm2(simple_cc1, simple_cc2, flip_neg=FALSE))
    expect_silent(plt_neg_flip <- splicegrahm2(simple_cc1, simple_cc2, flip_neg=TRUE))
    expect_silent(plt_neg_annot <- splicegrahm2(simple_cc1, simple_cc2, flip_neg=FALSE, txlist=exbytx))

    bld_neg_noflip <- ggplot_build(plt_neg_noflip@plot[[1]])
    bld_neg_flip <- ggplot_build(plt_neg_flip@plot[[1]])
    bld_neg_annot <- ggplot_build(plt_neg_annot@plot[[1]])

    ## check that coord flipped w/ neg strand
    expect_equal(bld_neg_noflip$panel$ranges[[1]]$x.range,
                 (-1) * rev(bld_neg_flip$panel$ranges[[1]]$x.range))
    
    ## check that coord not flipped w/ neg strand and neg annot if FALSE
    expect_equal(bld_neg_noflip$panel$ranges[[1]]$x.range,
                 bld_neg_annot$panel$ranges[[1]]$x.range)
})


test_that("splicegrahm2 accepts use_blk input to invert background color", {
    ## check that parameter is accepted without error/warnings
    expect_silent(plt_blk <- splicegrahm2(simple_cc1, simple_cc2, j_incl = TRUE, use_blk = TRUE))
    bld_blk_p1 <- ggplot_build(plt_blk@plot[[1]])
    bld_blk_p2 <- ggplot_build(plt_blk@plot[[2]])

    ## check that parameter is accepted without error/warnings
    plt_base <- splicegrahm2(simple_cc1, simple_cc2, j_incl = TRUE)
    bld_base_p1 <- ggplot_build(plt_base@plot[[1]])
    bld_base_p2 <- ggplot_build(plt_base@plot[[2]])

    ## check that heatmap frames are not drawn in darker plot (2 less 'GeomRect's)
    layers_base_p1 <- sapply(bld_base_p1$plot$layer, function(x) class(x$geom)[1])
    layers_blk_p1 <- sapply(bld_blk_p1$plot$layer, function(x) class(x$geom)[1])
    expect_equal(sum(layers_base_p1 == 'GeomRect') - sum(layers_blk_p1 == 'GeomRect'), 2)
    
    ## check that junction label colors have changed, everything else the same
    text_base <- do.call(rbind, bld_base_p1$data[layers_base_p1 == 'GeomText'])
    text_blk <- do.call(rbind, bld_blk_p1$data[layers_blk_p1 == 'GeomText'])
    expect_equal(text_base[, names(text_base) != 'colour'],
                 text_blk[, names(text_blk) != 'colour'])
    expect_true(all(text_base$colour != text_blk$colour))
})


test_that("splicegrahm2 accepts txlist, txdb, orgdb input to add gene annotations", {
    ## only including basic tests - parsing tested more w/ splicegrahm cases
    plt_base <- splicegrahm2(simple_cc1, simple_cc2, j_incl=TRUE)
    bld_base_p1 <- ggplot_build(plt_base@plot[[1]])
    bld_base_p2 <- ggplot_build(plt_base@plot[[2]])

    ## create a plot with a transcript annotations track
    expect_silent(plt_tx <- splicegrahm2(simple_cc1, simple_cc2, txlist=exbytx,
                                         txdb=txdb, orgdb=org.Hs.eg.db, j_incl=TRUE))
    expect_is(plt_tx, "Tracks")
    expect_true(length(plt_tx@plot) == 3)
    expect_is(plt_tx@plot[[1]], "ggplot")
    expect_is(plt_tx@plot[[2]], "GGbio")
    expect_is(plt_tx@plot[[3]], "ggplot")
    
    bld_tx_p1 <- ggplot_build(plt_tx@plot[[1]])
    bld_tx_p2 <- ggplot_build(plt_tx@plot[[2]])
    bld_tx_p3 <- ggplot_build(plt_tx@plot[[3]])
    
    ## check that top, bottom tracks are same as basic plot (at least data and layers)
    expect_equivalent(bld_tx_p1$data, bld_base_p1$data)
    expect_equivalent(bld_tx_p1$plot$layers, bld_base_p1$plot$layers)
    expect_equivalent(bld_tx_p3$data, bld_base_p2$data)
    expect_equivalent(bld_tx_p3$plot$layers, bld_base_p2$plot$layers)

    ## create a plot with only transcripts fully contained in concomp range (eps=0)
    expect_silent(plt_tx_e0 <- splicegrahm2(simple_cc1, simple_cc2, txlist=exbytx, eps=0))
    expect_is(plt_tx_e0, "Tracks")
    bld_tx_e0_p2 <- ggplot_build(plt_tx_e0@plot[[2]]) 

    ## check that eps=0 only keeps fully contained transcript models (just one here)
    expect_equal(bld_tx_e0_p2$panel$ranges[[1]]$y.labels, "70043")
})


test_that("splicegrahm2 accepts title parameter to add title to plot", {
    ## default plot with simple dataset
    plt_base <- splicegrahm2(simple_cc1, simple_cc2)

    ## plot with title
    expect_silent(plt_ttl <- splicegrahm2(simple_cc1, simple_cc2, title = "Simple Title"))

    ## check that title only exists when explicitly specified
    expect_equal(plt_base@main, "")
    expect_equal(plt_ttl@main, "Simple Title")
})


test_that("splicegrahm2 accepts mirror parameter for flipping bottom plot", {
    ## plot with bottom plot not mirrored (default mirror = TRUE)
    expect_silent(plt_nom <- splicegrahm2(simple_cc1, simple_cc1, mirror = FALSE))

    ## check that top and bottom panels are the same
    expect_equal(plt_nom@plot[[1]], plt_nom@plot[[2]])
})


test_that("splicegrahm2 accepts same_scale parameter for keeping plots on same scale", {
    ## plot with bottom plot not mirrored (default same_scale = TRUE), use mirror to make easier
    expect_silent(plt_base <- splicegrahm2(obj1=simple_cc1, obj2=simple_cc3,
                                           mirror = FALSE, same_scale = TRUE))
    expect_silent(plt_dscal <- splicegrahm2(obj1=simple_cc1, obj2=simple_cc3,
                                            mirror = FALSE, same_scale = FALSE))

    plt_1 <- splicegrahm(simple_cc1)
    plt_2 <- splicegrahm(simple_cc3)
    
    bld_base_p1 <- ggplot_build(plt_base@plot[[1]])
    bld_base_p2 <- ggplot_build(plt_base@plot[[2]])
    bld_dscal_p1 <- ggplot_build(plt_dscal@plot[[1]])
    bld_dscal_p2 <- ggplot_build(plt_dscal@plot[[2]])
    bld_1 <- ggplot_build(plt_1)
    bld_2 <- ggplot_build(plt_2)

    ## check that y ranges are same when same_scale = TRUE
    expect_equal(bld_base_p1$panel$ranges[[1]]$y.range,
                 bld_base_p2$panel$ranges[[1]]$y.range)

    ## check that y ranges are same as indep splicegrahm when same_scale = FALSE
    expect_equal(bld_dscal_p1$panel$ranges[[1]]$y.range,
                 bld_1$panel$ranges[[1]]$y.range)
    expect_equal(bld_dscal_p2$panel$ranges[[1]]$y.range,
                 bld_2$panel$ranges[[1]]$y.range)
})


