test_plogis <- function() {
    set.seed(666)
    x <- seq(-100, 100, 0.01)
    x <- x + runif(length(x), -0.001, 0.001)
    checkEqualsNumeric(plogis(x), bsseq:::.oldTrans(x))
    # NOTE: There was a bug in .oldTrans() when x = 0 that mean it returned 0
    #       instead of 0.5
    checkEqualsNumeric(bsseq:::.oldTrans(0), 0)
    checkEqualsNumeric(plogis(0), 0.5)
}
