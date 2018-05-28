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

.zero_type <- function(type) {
    if (identical(type, "integer")) {
        fill <- 0L
    } else if (identical(type, "double")) {
        fill <- 0
    } else {
        stop("'type' = ", type, " is not supported")
    }
}
