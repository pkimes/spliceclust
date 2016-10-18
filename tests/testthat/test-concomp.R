context("concomp constructor methods")
library("GenomicRanges")


test_that("concomp constructor accepts included GRanges data object", {
    ## load included dataset
    data("klk12_lusc_gr")

    ## split GRanges object to GrangesList by 'kind' (e, j)
    klk12_gl <- split(klk12_lusc_gr, klk12_lusc_gr$kind)

    ## check that example dataset generates a valid concomp
    expect_silent(klk12_cc <- concomp(klk12_gl))
    expect_is(klk12_cc, "concomp")
})


test_that("concomp constructor accepts empty GRangesList with 'e' and 'j' entries", {
    ## simple empty GRangesList with "e" and "j" entries
    empty_grl <- GRangesList(e = GRanges(),
                             j = GRanges())

    ## check that a valid concomp is constructed
    expect_silent(empty_cc <- concomp(empty_grl))
    expect_is(empty_cc, "concomp")
})


test_that("concomp constructor does not accept GRangesList with too few entries", {
    ## simple emty GRangesList with missing "j" entry 
    short_grl <- GRangesList(e = GRanges())

    ## check that concomp returns an error 
    expect_error(concomp(short_grl),
                 "GRangesList must only contain two objects")
})


test_that("concomp constructor does not accept GRangesList with too many entries", {
    ## simple emty GRangesList with extra "k" entry 
    excess_grl <- GRangesList(e = GRanges(),
                              j = GRanges(),
                              k = GRanges())

    ## check that concomp returns an error 
    expect_error(concomp(excess_grl),
                 "GRangesList must only contain two objects")
})


test_that("concomp constructor accepts data.frame", {
    ## simple data.frame with 2 samples
    simple_df <- data.frame(chr = rep("chr9", 3),
                            seqlengths = rep(1e4, 3),
                            gStart = c(100, 200, 300),
                            gStop = c(199, 299, 399),
                            kind = c("e", "j", "e"),
                            s1 = c(20, 30, 20),
                            s2 = c(100, 120, 30))

    ## check that a valid concomp is constructed
    expect_silent(simple_cc <- concomp(simple_df))
    expect_is(simple_cc, "concomp")

    ## check that concomp slots match input values
    expect_is(exons(simple_cc), "GRanges")
    expect_equal(ranges(exons(simple_cc)),
                 IRanges(start=c(100, 300),
                         end=c(199, 399)))

    expect_is(juncs(simple_cc), "GRanges")
    expect_equal(ranges(juncs(simple_cc)),
                 IRanges(start=200, end=299))

    expect_is(exonValues(simple_cc), "matrix")
    expect_equal(exonValues(simple_cc),
                 as.matrix(data.frame(s1=c(20, 20),
                                      s2=c(100, 30))))
    
    expect_is(juncValues(simple_cc), "matrix")
    expect_equal(juncValues(simple_cc),
                 as.matrix(data.frame(s1=30,
                                      s2=120)))
})


test_that("concomp constructor does not accept malformed data.frame", {
    ## dataset missing 'chr' column
    nochr_df <- data.frame(gStart=100, gStop=199, kind="e",
                           s1=20, s2=100)
    
    ## dataset missing 'kind' column
    nokind_df <- data.frame(chr="chr9", seqlengths=1e4,
                            gStart=100, gStop=199,
                            s1=20, s2=100)

    ## dataset missing 'seqlengths column
    noseqlen_df <- data.frame(chr="chr9",
                              gStart=100, gStop=199, kind="e",
                              s1=20, s2=100)

    ## check that correct error messages are returned
    expect_error(concomp(nochr_df),
                 "'chr' column must be specified with e.g. 'chr9'.")
    expect_error(concomp(nokind_df),
                 paste0("'kind' column must be specified with 'e' and 'j'",
                        " to differentiate exons and juncs."))
    expect_error(concomp(noseqlen_df),
                 paste0("'seqlengths' column must be specified, see manual for",
                        " directions on how to obtain appropriate values,",
                        " note that only the first value is used."))
})









