# Functions/methods that would be good to have in DelayedArray

.rowVars <- function(x, rows = NULL, cols = NULL, ...) {
    if (is(x, "DelayedArray")) {
        if (!is.null(rows)) {
            x <- x[rows, ]
        }
        if (!is.null(cols)) {
            x <- x[, cols]
        }
        row_vars <- rowVars(as.array(x), ...)
    } else {
        row_vars <- rowVars(x, rows = rows, cols = cols, ...)
    }
    row_vars
}

.rowSds <- function(x, rows = NULL, cols = NULL, ...) {
    row_vars <- .rowVars(x, rows = rows, cols = cols, ...)
    sqrt(row_vars)
}

.quantile <- function(x, ...) {
    if (is(x, "DelayedArray")) {
        x <- as.array(x)
    }
    quantile(x, ...)
}

.DelayedMatrix <- function(x) {
    x_name <- deparse(substitute(x))
    X <- try(DelayedArray(x), silent = TRUE)
    if (is(X, "try-error")) {
        stop("Could not construct DelayedMatrix from '", x_name, "'",
             call. = FALSE)
    }
    if (!is(X, "DelayedMatrix")) {
        stop("'", x_name, "' must be matrix-like", call. = FALSE)
    }
    X
}

.isSimpleDelayedMatrix <- function(x) {
    is(x@seed, "matrix")
}

# NOTE: Equivalent to rowSums2(x[, j, drop = FALSE]) but does it using a
#       delayed operation and always returns a nrow(x) x 1 DelayedMatrix
.delayed_rowSums2 <- function(x, j) {
    Reduce(`+`, lapply(j, function(jj) x[, jj, drop = FALSE]))
}

# NOTE: Equivalent to colSums2(x[i, , drop = FALSE]) but does it using a
#       delayed operation and always returns a 1 x ncol(x) DelayedMatrix
.delayed_colSums2 <- function(x, i) {
    Reduce(`+`, lapply(i, function(ii) x[ii, , drop = FALSE]))
}

# MARGIN = 1: collapse using rowSums
# MARGIN = 2: collapse using colSums
.collapseDelayedMatrix <- function(x, sp, MARGIN, BACKEND = NULL) {
    stopifnot(is(x, "DelayedMatrix"))
    if (MARGIN == 1) {
        if (is.null(BACKEND)) {
            collapsed_x <- do.call(cbind, lapply(sp, function(j) {
                rowSums2(x[, j, drop = FALSE])
            }))
        } else {
            collapsed_x <- do.call(cbind, lapply(sp, function(j) {
                .delayed_rowSums2(x, j)
            }))
            # NOTE: Need to manually add colnames when using this method
            colnames(collapsed_x) <- names(sp)
        }
    } else if (MARGIN == 2) {
        if (is.null(BACKEND)) {
            collapsed_x <- do.call(rbind, lapply(sp, function(i) {
                colSums2(x[i, , drop = FALSE])
            }))
        } else {
            collapsed_x <- do.call(rbind, lapply(sp, function(i) {
                .delayed_colSums2(x, i)
            }))
            # NOTE: Need to manually add rownames when using this method
            rownames(collapsed_x) <- names(sp)
        }
    } else {
        stop("'MARGIN' must be 1 or 2")
    }
    realize(collapsed_x, BACKEND = BACKEND)
}
