context("splicepca plotting method")
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
seqinfo(exons(simple_cc)) <- seqinfo(exbytx)[seqlevels(exons(simple_cc))]
seqinfo(juncs(simple_cc)) <- seqinfo(exbytx)[seqlevels(juncs(simple_cc))]


## ##############################################################################
## actual test cases

test_that("splicepca default call works", {
    ## default plot with simple dataset
    expect_silent(plt <- splicepca(simple_cc))

    ## default number of PCs
    npc <- 3
    
    ## check that ggplot is produced
    expect_is(plt, "ggplot")
    bld <- ggplot_build(plt)

    ## check ranges of plot - allow +/-1 wiggle room
    expect_equal((-1) * bld$panel$ranges[[1]]$x.range[2], 51532000 - 300, tolerance=1)
    expect_equal((-1) * bld$panel$ranges[[1]]$x.range[1], 51538000 + 300, tolerance=1)

    ## determine layer geoms
    layers_base <- sapply(plt$layer, function(x) class(x$geom)[1])

    ## check that one GeomPaths per junction, per PC (3 by default)
    expect_equal(sum(layers_base == "GeomPath"),
                 npc * length(juncs(simple_cc)))
    
    ## check that 2 GeomRects drawn: 1 for exons, 1 for panels
    expect_equal(sum(layers_base == "GeomRect"), 2)

    ## check that exons GeomRect draws correct # exons in each panel
    epanels <- bld$data[[1]]$PANEL
    expect_is(epanels, "factor")
    expect_equal(levels(epanels), as.character(1:npc))
    expect_true(all(table(epanels) == nrow(exons(simple_cc))))

    ## check that text labels are included as single geom, 1 text per panel
    lyr_text <- which(layers_base == "GeomText")
    dat_text <- bld$data[[lyr_text]]
    expect_length(lyr_text, 1)
    expect_true(all(grepl("EXON var expl = ", dat_text$label)))
    expect_true(all(grepl("JUNC var expl = ", dat_text$label)))
    expect_true(all(table(dat_text$PANEL) == 1))    
})


test_that("splicepca PCA decomposition can be modified with npc, pc_sep, ej_w parameters", {
    ## default plot with simple dataset
    plt_base <- splicepca(simple_cc)

    ## create plots with pc_sep = TRUE (default)
    expect_silent(plt_npc2 <- splicepca(simple_cc, npc = 2))
    expect_silent(plt_ej10b <- splicepca(simple_cc, ej_w = c(10, 0)))
    ## create plots with pc_sep = FALSE 
    expect_silent(plt_joint <- splicepca(simple_cc, pc_sep = FALSE))
    expect_silent(plt_npc6 <- splicepca(simple_cc, npc = 6, pc_sep = FALSE))
    expect_silent(plt_ej01 <- splicepca(simple_cc, pc_sep = FALSE, ej_w = c(0, 1)))
    expect_silent(plt_ej10 <- splicepca(simple_cc, pc_sep = FALSE, ej_w = c(1, 0)))
    expect_silent(plt_ej55 <- splicepca(simple_cc, pc_sep = FALSE, ej_w = c(5, 5)))

    ## parse plots for analysis
    expect_is(plt_base, "ggplot")
    expect_is(plt_npc2, "ggplot")
    expect_is(plt_joint, "ggplot")
    expect_is(plt_npc6, "ggplot")
    expect_is(plt_ej01, "ggplot")
    expect_is(plt_ej10, "ggplot")
    expect_is(plt_ej55, "ggplot")
    expect_is(plt_ej10b, "ggplot")
    bld_base <- ggplot_build(plt_base)
    bld_npc2 <- ggplot_build(plt_npc2)
    bld_joint <- ggplot_build(plt_joint)
    bld_npc6 <- ggplot_build(plt_npc6)
    bld_ej01 <- ggplot_build(plt_ej01)
    bld_ej10 <- ggplot_build(plt_ej10)
    bld_ej55 <- ggplot_build(plt_ej55)
    bld_ej10b <- ggplot_build(plt_ej10b)

    ## determine layer geoms
    layers_base <- sapply(plt_base$layer, function(x) class(x$geom)[1])
    layers_npc2 <- sapply(plt_npc2$layer, function(x) class(x$geom)[1])
    layers_joint <- sapply(plt_joint$layer, function(x) class(x$geom)[1])

    layers_npc6 <- sapply(plt_npc6$layer, function(x) class(x$geom)[1])
    layers_ej01 <- sapply(plt_ej01$layer, function(x) class(x$geom)[1])
    layers_ej10 <- sapply(plt_ej10$layer, function(x) class(x$geom)[1])
    layers_ej55 <- sapply(plt_ej55$layer, function(x) class(x$geom)[1])
    layers_ej10b <- sapply(plt_ej10b$layer, function(x) class(x$geom)[1])

    ## check that npc parameter adjusts number of panels shown
    ## check that one GeomPaths per junction, per PC
    expect_equal(sum(layers_npc2 == "GeomPath"), 2 * length(juncs(simple_cc)))
    expect_equal(sum(layers_npc6 == "GeomPath"), 6 * length(juncs(simple_cc)))
    ## check that exons GeomRect draws exons in npc panel
    expect_equal(levels(bld_npc2$data[[1]]$PANEL), as.character(1:2))
    expect_equal(levels(bld_npc6$data[[1]]$PANEL), as.character(1:6))
    
    ## check that PCs 1,2 same when npc = 2, 3 (expect PC assoc. levels)
    expect_equivalent(bld_base$data[[which(layers_base == "GeomText")]][1:2, ],
                      bld_npc2$data[[which(layers_npc2 == "GeomText")]])
    
    ## check that PCs 1,2,3 also same when pc_sep = FALSE
    expect_equivalent(bld_joint$data[[which(layers_joint == "GeomText")]],
                      bld_npc6$data[[which(layers_npc6 == "GeomText")]][1:3, ])

    ## check that specifying pc_sep = FALSE generates different PC decomp coloring
    expect_false(all(bld_base$data[[1]]$fill == bld_joint$data[[1]]$fill))
    
    ## check that specifying pc_sep = FALSE generates labels with only total var
    expect_true(all(grepl("TOTAL var expl = ",
                          bld_joint$data[[which(layers_joint == "GeomText")]]$label)))

    ## check that placing all weight on exons plots juncs w/ null color
    arcs_ej10 <- do.call(rbind, bld_ej10$data[which(layers_ej10 == "GeomPath")])
    expect_true(all(arcs_ej10$fill == "#F7F7F7"))

    ## check that placing all weight on juncs plots exons w/ null color
    expect_true(all(bld_ej01$data[[1]]$fill == "#F7F7F7"))
    
    ## check that specified ej_w weights are scaled such that c(5, 5) same as default
    expect_equivalent(bld_joint, bld_ej55)
    
    ## check that specified ej_w weights w/ pc_sep = TRUE is ignored
    expect_equivalent(bld_base, bld_ej10b)
    
    ## check that warning returned if too many PCs specified (more than rank)
    expect_warning(plt_npc4 <- splicepca(simple_cc, npc = 4, pc_sep = TRUE),
                   "npc larger than dim of data, using npc = 3")
    expect_warning(plt_npc7 <- splicepca(simple_cc, npc = 7, pc_sep = FALSE),
                   "npc larger than dim of data, using npc = 6")

    ## check that after warning
    expect_equivalent(plt_npc4, plt_base)
    expect_equivalent(plt_npc7, plt_npc6)
})
    

test_that("splicepca accepts log_base, log_shift to scale heatmap colorscale", {
    ## create default plot with log_base = 10, log_shift = 1
    plt_base <- splicepca(simple_cc)

    ## create plots based on log scaling and shift - loadings plots
    expect_silent(plt_b2s1 <- splicepca(simple_cc, log_base=10, log_shift=1))
    expect_silent(plt_b2s0 <- splicepca(simple_cc, log_base=10, log_shift=0))
    expect_silent(plt_b0s1 <- splicepca(simple_cc, log_base=0, log_shift=1))
    expect_silent(plt_b0s0 <- splicepca(simple_cc, log_base=0, log_shift=0))
    
    expect_is(plt_b2s1, "ggplot")
    expect_is(plt_b2s0, "ggplot")
    expect_is(plt_b0s1, "ggplot")
    expect_is(plt_b0s0, "ggplot")
    bld_base <- ggplot_build(plt_base)
    bld_b2s1 <- ggplot_build(plt_b2s1)
    bld_b2s0 <- ggplot_build(plt_b2s0)
    bld_b0s1 <- ggplot_build(plt_b0s1)
    bld_b0s0 <- ggplot_build(plt_b0s0)

    layers_base <- sapply(plt_base$layer, function(x) class(x$geom)[1])
    layers_b2s1 <- sapply(plt_b2s1$layer, function(x) class(x$geom)[1])
    layers_b2s0 <- sapply(plt_b2s0$layer, function(x) class(x$geom)[1])
    layers_b0s1 <- sapply(plt_b0s1$layer, function(x) class(x$geom)[1])
    layers_b0s0 <- sapply(plt_b0s0$layer, function(x) class(x$geom)[1])
    
    ## check that PC loadings changed (i.e. decomposition changed) when log_shift changed
    expect_true(all(as.character(bld_base$data[[which(layers_base == "GeomText")]]$label) != 
                        as.character(bld_b2s0$data[[which(layers_b2s1 == "GeomText")]]$label)))

    ## check that PC loadings changed when no log transformation is used
    expect_true(all(as.character(bld_base$data[[which(layers_base == "GeomText")]]$label) != 
                        as.character(bld_b0s0$data[[which(layers_b2s1 == "GeomText")]]$label)))

    ## check that PC loadings are the same when only log base is changed
    expect_true(all(as.character(bld_base$data[[which(layers_base == "GeomText")]]$label) == 
                        as.character(bld_b2s1$data[[which(layers_b2s1 == "GeomText")]]$label)))

    ## check that log_shift is unused when log_base = 0
    expect_equivalent(bld_b0s1, bld_b0s0)
    

    ## create plots based on log scaling and shft - scores plots 
    expect_silent(plt_bases <- splicepca(simple_cc, scores=TRUE))
    expect_silent(plt_b2s1s <- splicepca(simple_cc, log_base=2, log_shift=1, scores=TRUE))

    ## check that matrix plots are returned when scores = TRUE
    expect_is(plt_bases, "ggmatrix")
    expect_is(plt_b2s1s, "ggmatrix")

    ## check that PC scores are different when log_base is changed
    expect_true(all(plt_bases$data != plt_b2s1s$data))
})


test_that("splicepca accepts genomic, ex_use input to adjust x-axis scaling", {
    ## create default plot with genomic = TRUE 
    plt_base <- splicepca(simple_cc)

    ## create plots on non-genomic scale
    expect_silent(plt_non <- splicepca(simple_cc, genomic = FALSE))
    expect_silent(plt_non_10 <- splicepca(simple_cc, genomic = FALSE, ex_use = 1.0))

    ## plot with ex_use less than exon proportion in genomic space
    expect_warning(plt_non_01 <- splicepca(simple_cc, genomic = FALSE, ex_use = 0.1),
                   paste0("Using 'genomic = TRUE' since exons occupy more than specified 'ex_use' ",
                          "proportion of plot in genomic coordinates. No need to squish."))

    ## check that ggplot is produced
    expect_is(plt_base, "ggplot")
    expect_is(plt_non, "ggplot")
    expect_is(plt_non_10, "ggplot")
    expect_is(plt_non_01, "ggplot")
    bld_base <- ggplot_build(plt_base)
    bld_non <- ggplot_build(plt_non)
    bld_non_10 <- ggplot_build(plt_non_10)
    bld_non_01 <- ggplot_build(plt_non_01)

    ## check that coordinate range shrinks in expected proportions
    width_base <- diff(bld_base$panel$ranges[[1]]$x.range)
    width_non <- diff(bld_non$panel$ranges[[1]]$x.range)
    width_non_10 <- diff(bld_non_10$panel$ranges[[1]]$x.range)
    expect_lt(width_non, width_base)
    expect_lt(width_non_10, width_base)
    expect_equal(width_non_10 / width_non, 2/3, tolerance=0.001)
    
    ## check that bad (small) ex_use creates same plot as genomic = TRUE
    expect_equal(bld_base$data, bld_non_01$data)
    expect_equal(bld_base$panel, bld_non_01$panel)
    expect_equal(plt_base[names(plt_base) != "plot_env"],
                 plt_non_01[names(plt_non_01) != "plot_env"])
})


test_that("splicepca accepts flip_neg to adjust whether stranded-ness matters", {
    ## only including basic tests - flip_neg should be refactored out
    expect_silent(plt_neg_noflip <- splicepca(simple_cc, flip_neg=FALSE))
    expect_silent(plt_neg_flip <- splicepca(simple_cc, flip_neg=TRUE))
    expect_silent(plt_neg_annot <- splicepca(simple_cc, flip_neg=FALSE, txlist=exbytx))

    expect_is(plt_neg_noflip, "ggplot")
    expect_is(plt_neg_flip, "ggplot")
    expect_is(plt_neg_annot, "Tracks")
    bld_neg_noflip <- ggplot_build(plt_neg_noflip)
    bld_neg_flip <- ggplot_build(plt_neg_flip)
    bld_neg_annot <- ggplot_build(plt_neg_annot@plot[[1]])

    ## check that coord flipped w/ neg strand
    expect_equal(bld_neg_noflip$panel$ranges[[1]]$x.range,
                 (-1) * rev(bld_neg_flip$panel$ranges[[1]]$x.range))
    
    ## check that coord not flipped w/ neg strand and neg annot if FALSE
    expect_equal(bld_neg_noflip$panel$ranges[[1]]$x.range,
                 bld_neg_annot$panel$ranges[[1]]$x.range)
})


test_that("splicepca accepts use_blk input to invert background color", {
    ## default plot with simple dataset
    plt_base <- splicepca(simple_cc, j_incl = TRUE)
    bld_base <- ggplot_build(plt_base)

    ## check that parameter is accepted
    expect_silent(plt_blk <- splicepca(simple_cc, use_blk = TRUE))
    bld_blk <- ggplot_build(plt_blk)

    ## check that exon frames are not drawn in darker plot (1 less 'GeomRect's)
    layers_base <- sapply(plt_base$layer, function(x) class(x$geom)[1])
    layers_blk <- sapply(plt_blk$layer, function(x) class(x$geom)[1])
    expect_equal(sum(layers_base == 'GeomRect') - sum(layers_blk == 'GeomRect'), 1)
    
    ## check that junction label colors have changed, everything else the same
    text_base <- do.call(rbind, bld_base$data[layers_base == 'GeomText'])
    text_blk <- do.call(rbind, bld_blk$data[layers_blk == 'GeomText'])
    expect_equal(text_base[, names(text_base) != 'colour'],
                 text_blk[, names(text_blk) != 'colour'])
    expect_true(all(text_base$colour != text_blk$colour))
})


test_that("splicepca accepts txlist, txdb, orgdb input to add gene annotations", {
    ## default plot with simple dataset
    plt_base <- splicepca(simple_cc)
    bld_base <- ggplot_build(plt_base)

    ## create a plot with a transcript annotations track
    expect_silent(plt_tx <- splicepca(simple_cc, txlist=exbytx,
                                      txdb=txdb, orgdb=org.Hs.eg.db))
    expect_is(plt_tx, "Tracks")
    expect_is(plt_tx@plot[[1]], "ggplot")
    expect_is(plt_tx@plot[[2]], "ggplot")
    bld_tx_p1 <- ggplot_build(plt_tx@plot[[1]])
    bld_tx_p2 <- ggplot_build(plt_tx@plot[[2]])
    
    ## check that top track is same as basic plot
    expect_equivalent(bld_tx_p1$data, bld_base$data)
    expect_equivalent(bld_tx_p1$plot$layers, plt_base$layers)
})


test_that("splicepca accepts scores parameter to plot scores, not loadings", {
    ## default plot with simple dataset
    plt_base <- splicepca(simple_cc)
    bld_base <- ggplot_build(plt_base)

    ## check that scores plot can be created
    expect_silent(plt_scores <- splicepca(simple_cc, scores=TRUE))
    expect_is(plt_scores, "ggmatrix")

    ## check that data used for scores plot is as expected
    expect_equal(names(plt_scores$data),
                 c(paste0("E_PC", 1:3), paste0("J_PC", 1:3)))
    expect_equal(nrow(plt_scores$data), 20)
    expect_equal(length(plt_scores$plots), (2 * 3)^2)
})


test_that("splicepca accepts plot parameter to determine whether to just return values", {
    ## default plot with simple dataset
    plt_base <- splicepca(simple_cc)
    bld_base <- ggplot_build(plt_base)

    ## check that values can be returned instead of plot
    expect_silent(plt_vals <- splicepca(simple_cc, plot = FALSE))
    expect_is(plt_vals, "list")
    expect_equal(names(plt_vals), c("pca_e", "pca_j"))
    expect_is(plt_vals$pca_e, "prcomp")
    expect_is(plt_vals$pca_j, "prcomp")

    ## check that PCA decompositions are different for exons, juncs (no accidental copying)
    expect_true(all(dim(plt_vals$pca_e$x) == c(20, length(exons(simple_cc)))))
    expect_true(all(dim(plt_vals$pca_j$x) == c(20, length(juncs(simple_cc)))))
    expect_true(any(plt_vals$pca_e$x != plt_vals$pca_j$x))
})



