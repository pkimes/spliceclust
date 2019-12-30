context("splicegrahm plotting method")
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

## fix corresponding reference genomes
seqinfo(eRanges(simple_cc)) <- seqinfo(exbytx)[seqlevels(eRanges(simple_cc))]
seqinfo(jRanges(simple_cc)) <- seqinfo(exbytx)[seqlevels(jRanges(simple_cc))]


## ##############################################################################
## actual test cases

test_that("splicegrahm default call works", {
    ## default plot with simple dataset
    expect_silent(plt <- splicegrahm(simple_cc))
    
    ## check that ggplot is produced
    expect_is(plt, "ggplot")
    bld <- ggplot_build(plt)

    ## check ranges of plot 
    expect_equal((-1) * bld$layout$panel_scales_x[[1]]$range$range[2], 51532000)
    expect_equal((-1) * bld$layout$panel_scales_x[[1]]$range$range[1], 51538000)

    ## determine layer geoms
    layers_base <- sapply(plt$layer, function(x) class(x$geom)[1])

    ## check one GeomPaths per junction
    expect_equal(sum(layers_base == "GeomPath"), length(jRanges(simple_cc)))
    
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
    expect_equal(sum(layers_base == "GeomPath"), length(jRanges(simple_cc)))
    
    ## check GeomRects for exon/junc heatmap, exon outlines, junc outlines
    expect_equal(sum(layers_base == "GeomRect"), 3)

    ## check two GeomText for junction labels
    expect_equal(sum(layers_base == "GeomText"), 2)

    ## check ranges of splice junctions (flip for negative strand)
    gp_idx <- which(layers_base == "GeomPath")
    gp_starts <- (-1) * sapply(gp_idx, function(idx) max(bld$data[[idx]]$x))
    gp_ends <- (-1) * sapply(gp_idx, function(idx) min(bld$data[[idx]]$x))
    expect_equal(start(jRanges(simple_cc)), gp_starts)
    expect_equal(end(jRanges(simple_cc)), gp_ends)

    ## check number of rects matchs #exons, #junctions, #heatmap boxes
    gr_idx <- which(layers_base == "GeomRect")
    gr_rows <- sort(sapply(gr_idx, function(idx) nrow(bld$data[[idx]])))
    expect_equal(gr_rows, sort(c(length(eRanges(simple_cc)), length(jRanges(simple_cc)),
                                 prod(dim(eCoverage(simple_cc))) +
                                     prod(dim(jCoverage(simple_cc))))))

    ## check heatmap boxes placed in correct positions
    exon_frame_color <- "#3C3C3C"
    gr_dat <- do.call(rbind, bld$data[gr_idx])
    gr_dat <- gr_dat[gr_dat$colour == exon_frame_color, ]
    gr_starts <- (-1) * gr_dat$xmin
    gr_ends <- (-1) * gr_dat$xmax
    expect_equal(start(eRanges(simple_cc)), gr_starts)
    expect_equal(end(eRanges(simple_cc)), gr_ends)

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
    width_gen <- diff(bld_gen$layout$panel_scales_x[[1]]$range$range)
    width_non <- diff(bld_non$layout$panel_scales_x[[1]]$range$range)
    width_non_10 <- diff(bld_non_10$layout$panel_scales_x[[1]]$range$range)
    expect_lt(width_non, width_gen)
    expect_lt(width_non_10, width_gen)
    expect_equal(width_non_10 / width_non, 2/3, tolerance=0.001)
    
    ## check that bad (small) ex_use creates same plot as genomic = TRUE
    expect_equal(bld_gen$data, bld_non_01$data)
    expect_equal(bld_gen$layout, bld_non_01$layout)
    expect_equal(plt_gen[names(plt_gen) != "plot_env"],
                 plt_non_01[names(plt_non_01) != "plot_env"])

    ## check that txlist creates Tracks plot
    expect_silent(plt_tx_10 <- splicegrahm(simple_cc, genomic = FALSE, ex_use = 1.0,
                                           txlist = exbytx))
    expect_is(plt_tx_10, "Tracks")

    ## check that plot is wider since adjustment accounts for gene models from txlist
    bld_tx_10_p1 <- ggplot_build(plt_tx_10@plot[[1]])
    width_tx_10 <- diff(bld_tx_10_p1$layout$panel_scales_x[[1]]$range$range)
    expect_gt(width_tx_10, width_non_10)
})


test_that("splicegrahm accepts flip_neg to adjust whether stranded-ness matters", {
    ## datasets with sequence on '+' or '*' strands
    pos_strand_cc <- simple_cc
    strand(eRanges(pos_strand_cc)) <- '+'
    strand(jRanges(pos_strand_cc)) <- '+'

    ## dataset with sequences on '*' strand
    unk_strand_cc <- simple_cc
    strand(eRanges(unk_strand_cc)) <- '*'
    strand(jRanges(unk_strand_cc)) <- '*'

    ## simple_cc uses negative strand annotations
    neg_strand_cc <- simple_cc

    ## default plot with simple dataset
    expect_silent(plt_pos_noflip <- splicegrahm(pos_strand_cc, j_incl=TRUE, flip_neg=FALSE))
    expect_silent(plt_pos_flip <- splicegrahm(pos_strand_cc, j_incl=TRUE, flip_neg=TRUE))
    expect_silent(plt_pos_annot <- splicegrahm(pos_strand_cc, j_incl=TRUE, flip_neg=TRUE, txlist=exbytx))
    bld_pos_noflip <- ggplot_build(plt_pos_noflip)
    bld_pos_flip <- ggplot_build(plt_pos_flip)
    bld_pos_annot <- ggplot_build(plt_pos_annot@plot[[1]])
    
    expect_silent(plt_neg_noflip <- splicegrahm(neg_strand_cc, j_incl=TRUE, flip_neg=FALSE))
    expect_silent(plt_neg_flip <- splicegrahm(neg_strand_cc, j_incl=TRUE, flip_neg=TRUE))
    expect_silent(plt_neg_annot <- splicegrahm(neg_strand_cc, j_incl=TRUE, flip_neg=FALSE, txlist=exbytx))
    bld_neg_noflip <- ggplot_build(plt_neg_noflip)
    bld_neg_flip <- ggplot_build(plt_neg_flip)
    bld_neg_annot <- ggplot_build(plt_neg_annot@plot[[1]])

    expect_silent(plt_unk_noflip <- splicegrahm(unk_strand_cc, j_incl=TRUE, flip_neg=FALSE))
    expect_silent(plt_unk_flip <- splicegrahm(unk_strand_cc, j_incl=TRUE, flip_neg=TRUE))
    expect_silent(plt_unk_annot <- splicegrahm(unk_strand_cc, j_incl=TRUE, flip_neg=TRUE, txlist=exbytx))
    bld_unk_noflip <- ggplot_build(plt_unk_noflip)
    bld_unk_flip <- ggplot_build(plt_unk_flip)
    bld_unk_annot <- ggplot_build(plt_unk_annot@plot[[1]])
    
    ## check that coord flipped w/ neg strand
    expect_equal(bld_neg_noflip$layout$panel_scales_x[[1]]$range$range,
                 (-1) * rev(bld_neg_flip$layout$panel_scales_x[[1]]$range$range))
    
    ## check that coord flipped w/ unknown strand and neg annot
    expect_equal(bld_unk_flip$layout$panel_scales_x[[1]]$range$range,
                 (-1) * rev(bld_unk_annot$layout$panel_scales_x[[1]]$range$range))

    ## check that coord not flipped w/ pos strand
    expect_equal(bld_pos_flip$layout$panel_scales_x[[1]]$range$range,
                 bld_pos_noflip$layout$panel_scales_x[[1]]$range$range)
    
    ## check that coord not flipped w/ unknown strand
    expect_equal(bld_unk_flip$layout$panel_scales_x[[1]]$range$range,
                 bld_unk_noflip$layout$panel_scales_x[[1]]$range$range)

    ## check that coord not flipped w/ pos strand and neg annot
    expect_equal(bld_pos_flip$layout$panel_scales_x[[1]]$range$range,
                 bld_pos_annot$layout$panel_scales_x[[1]]$range$range)

    ## check that coord not flipped w/ neg strand and neg annot if FALSE
    expect_equal(bld_neg_noflip$layout$panel_scales_x[[1]]$range$range,
                 bld_neg_annot$layout$panel_scales_x[[1]]$range$range)
})


test_that("splicegrahm accepts highlight input to add group annotations", {
    ## check that string labels
    hl_s1 <- rep("A", each=20)
    hl_s2 <- rep(c("A", "B"), each=10)
    hl_s3 <- c(rep(c("A", "B"), each=6), rep("C", 8))
    hl_s4 <- rep(c("A", "B", "C", "D"), each=5)
    hl_f2 <- factor(hl_s2)
    hl_i2 <- rep(1:2, each=10)
    hl_short <- rep("A", 12)
    hl_long <- rep("A", 35)
    hl_reord <- rep(c("A", "B", "C", "D"), 5)

    ## simple plot w/ no highlighting
    expect_silent(plt_base <- splicegrahm(simple_cc, j_inc=TRUE))

    ## check that reasonable highlight labels are accepted
    expect_silent(plt_s1 <- splicegrahm(simple_cc, j_incl=TRUE, highlight=hl_s1))
    expect_silent(plt_s2 <- splicegrahm(simple_cc, j_incl=TRUE, highlight=hl_s2))
    expect_silent(plt_s3 <- splicegrahm(simple_cc, j_incl=TRUE, highlight=hl_s3))

    ## check that too many levels still accepted (>3)
    expect_silent(plt_s4 <- splicegrahm(simple_cc, j_incl=TRUE, highlight=hl_s4))

    ## check that factors can be specified
    expect_silent(plt_f2 <- splicegrahm(simple_cc, j_incl=TRUE, highlight=hl_f2))

    ## check that integers can be specified
    expect_silent(plt_i2 <- splicegrahm(simple_cc, j_incl=TRUE, highlight=hl_i2))
    
    ## chech that incorrect length vector still accepted
    expect_silent(plt_short <- splicegrahm(simple_cc, j_incl=TRUE, highlight=hl_short))
    expect_silent(plt_long <- splicegrahm(simple_cc, j_incl=TRUE, highlight=hl_long))

    ## check that samples are reordered by highlight
    expect_silent(plt_reord <- splicegrahm(simple_cc, j_incl=TRUE, highlight=hl_reord))

    ## layers in ggplot
    layers_base <- sapply(plt_base$layer, function(x) class(x$geom)[1])
    layers_s1 <- sapply(plt_s1$layer, function(x) class(x$geom)[1])
    layers_s2 <- sapply(plt_s2$layer, function(x) class(x$geom)[1])
    layers_s3 <- sapply(plt_s3$layer, function(x) class(x$geom)[1])
    layers_s4 <- sapply(plt_s4$layer, function(x) class(x$geom)[1])
    layers_f2 <- sapply(plt_f2$layer, function(x) class(x$geom)[1])
    layers_i2 <- sapply(plt_i2$layer, function(x) class(x$geom)[1])
    layers_short <- sapply(plt_short$layer, function(x) class(x$geom)[1])
    layers_long <- sapply(plt_long$layer, function(x) class(x$geom)[1])
    layers_reord <- sapply(plt_reord$layer, function(x) class(x$geom)[1])

    ## check that all have more GeomRects (highlight boxes) than base plot
    expect_equal(sum(layers_s1 == "GeomRect") - sum(layers_base == "GeomRect"),
                 2 * length(unique(hl_s1)))
    expect_equal(sum(layers_s2 == "GeomRect") - sum(layers_base == "GeomRect"),
                 2 * length(unique(hl_s2)))
    expect_equal(sum(layers_s3 == "GeomRect") - sum(layers_base == "GeomRect"),
                 2 * length(unique(hl_s3)))
    expect_equal(sum(layers_s4 == "GeomRect") - sum(layers_base == "GeomRect"),
                 2 * length(unique(hl_s4)))
    expect_equal(sum(layers_short == "GeomRect") - sum(layers_base == "GeomRect"),
                 2 * length(unique(hl_short)))
    expect_equal(sum(layers_long == "GeomRect") - sum(layers_base == "GeomRect"),
                 2 * length(unique(hl_long)))
    expect_equal(sum(layers_reord == "GeomRect") - sum(layers_base == "GeomRect"),
                 2 * length(unique(hl_reord)))

    ## check that plots are same using character, factor, integer vector
    expect_equivalent(plt_s2[names(plt_s2) != "plot_env"],
                      plt_f2[names(plt_f2) != "plot_env"])
    expect_equivalent(plt_s2[names(plt_s2) != "plot_env"],
                      plt_i2[names(plt_i2) != "plot_env"])

    ## check that heatmap is reordered to group by highlight parameter
    dat_e <- unique(plt_reord$data[plt_reord$data$kind == "e", c("variable", "ymin")])
    dat_j <- unique(plt_reord$data[plt_reord$data$kind == "j", c("variable", "ymin")])
    dat_e <- as.character(dat_e$variable[order(dat_e$ymin)])
    dat_j <- as.character(dat_j$variable[order(dat_j$ymin)])
    expect_equal(dat_e, paste0("s", order(hl_reord)))
    expect_equal(dat_j, paste0("s", order(hl_reord)))
})


test_that("splicegrahm accepts use_blk input to invert background color", {
    ## default plot with simple dataset
    plt_base <- splicegrahm(simple_cc, j_incl = TRUE)
    bld_base <- ggplot_build(plt_base)

    ## check that parameter is accepted
    expect_silent(plt_blk <- splicegrahm(simple_cc, j_incl = TRUE, use_blk = TRUE))
    bld_blk <- ggplot_build(plt_blk)

    ## check that heatmap frames are not drawn in darker plot (2 less 'GeomRect's)
    layers_base <- sapply(plt_base$layer, function(x) class(x$geom)[1])
    layers_blk <- sapply(plt_blk$layer, function(x) class(x$geom)[1])
    expect_equal(sum(layers_base == 'GeomRect') - sum(layers_blk == 'GeomRect'), 2)
    
    ## check that junction label colors have changed, everything else the same
    text_base <- do.call(rbind, bld_base$data[layers_base == 'GeomText'])
    text_blk <- do.call(rbind, bld_blk$data[layers_blk == 'GeomText'])
    expect_equal(text_base[, names(text_base) != 'colour'],
                 text_blk[, names(text_blk) != 'colour'])
    expect_true(all(text_base$colour != text_blk$colour))
})


test_that("splicegrahm accepts eps, txlist, txdb, orgdb input to add gene annotations", {
    ## default plot with simple dataset
    plt_base <- splicegrahm(simple_cc, j_incl=TRUE)
    bld_base <- ggplot_build(plt_base)
    

    ## create a plot with a transcript annotations track
    expect_silent(plt_tx <- splicegrahm(simple_cc, j_incl=TRUE, txlist=exbytx))
    expect_is(plt_tx, "Tracks")
    expect_is(plt_tx@plot[[1]], "ggplot")
    expect_is(plt_tx@plot[[2]], "GGbio")
    bld_tx_p1 <- ggplot_build(plt_tx@plot[[1]])
    bld_tx_p2 <- ggplot_build(plt_tx@plot[[2]]@ggplot)
    
    ## check that top track is same as basic plot
    expect_equivalent(bld_tx_p1$data, bld_base$data)
    expect_equivalent(bld_tx_p1$plot$layers, plt_base$layers)

    ## check that gene models are included by checking y labels
    expect_equal(bld_tx_p2$layout$panel_scales_y[[1]]$labels,
                 c("70043", "70044", "70045", "70046",
                   "70047", "70048", "70049", "70050"))

    
    ## create a plot with only transcripts fully contained in concomp range (eps=0)
    expect_silent(plt_tx_e0 <- splicegrahm(simple_cc, j_incl=TRUE, txlist=exbytx, eps=0))
    expect_is(plt_tx_e0, "Tracks")
    bld_tx_e0_p2 <- ggplot_build(plt_tx_e0@plot[[2]]@ggplot) 

    ## check that eps=0 only keeps fully contained transcript models
    expect_equal(bld_tx_e0_p2$layout$panel_scales_y[[1]]$labels, "70043")

    
    ## create a plot with named transcript track (w/ txdb)
    expect_silent(plt_txdb <- splicegrahm(simple_cc, j_incl=TRUE, txlist=exbytx, txdb=txdb))
    expect_is(plt_txdb, "Tracks")
    bld_txdb_p1 <- ggplot_build(plt_txdb@plot[[1]])
    bld_txdb_p2 <- ggplot_build(plt_txdb@plot[[2]]@ggplot)

    ## check that top track is same as basic plot
    expect_equivalent(bld_txdb_p1$data, bld_base$data)
    expect_equivalent(bld_txdb_p1$plot$layers, plt_base$layers)

    ## check that txdb gives transcript IDs as y-axis labels
    expect_equal(bld_txdb_p2$layout$panel_scales_y[[1]]$labels,
                 c("uc002pvg.1", "uc010ycp.1", "uc010ycq.1",
                   "uc010ycr.1", "uc010ycs.1", "uc002pvh.1",
                   "uc002pvi.1", "uc002pvj.1"))

    
    ## create a plot with named transcript track (w/ txdb)
    expect_silent(plt_orgdb <-
        splicegrahm(simple_cc, j_incl=TRUE, txlist=exbytx, txdb=txdb, orgdb=org.Hs.eg.db))
    expect_is(plt_orgdb, "Tracks")
    bld_orgdb_p1 <- ggplot_build(plt_orgdb@plot[[1]])
    bld_orgdb_p2 <- ggplot_build(plt_orgdb@plot[[2]]@ggplot)

    ## check that top track is same as basic plot
    expect_equivalent(bld_orgdb_p1$data, bld_base$data)
    expect_equivalent(bld_orgdb_p1$plot$layers, plt_base$layers)

    ## check that orgdb gives gene symbols as y-axis labels
    expect_equal(bld_orgdb_p2$layout$panel_scales_y[[1]]$labels,
                 c("KLK12", paste0("KLK12_", 1:7)))    
})


test_that("splicegrahm accepts title parameter to add title to plot", {
    ## default plot with simple dataset
    plt_base <- splicegrahm(simple_cc)
    bld_base <- ggplot_build(plt_base)

    ## plot with title
    plt_ttl <- splicegrahm(simple_cc, title = "Simple Title")
    bld_ttl <- ggplot_build(plt_ttl)

    ## check that title only exists when explicitly specified
    expect_equal(plt_base$labels$title, "")
    expect_equal(plt_ttl$labels$title, "Simple Title")

    ## check that after 're-setting' title, plots are same
    plt_ttl$labels$title <- ""
    expect_equal(plt_base, plt_ttl)
})

