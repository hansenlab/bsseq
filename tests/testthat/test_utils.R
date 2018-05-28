context("utils")

test_that("plogis()", {
    set.seed(666)
    x <- seq(-100, 100, 0.01)
    x <- x + runif(length(x), -0.001, 0.001)
    expect_equal(plogis(x), bsseq:::.oldTrans(x))
    # NOTE: There was a bug in .oldTrans() when x = 0 that mean it returned 0
    #       instead of 0.5
    expect_equal(bsseq:::.oldTrans(0), 0)
    expect_equal(plogis(0), 0.5)
})
