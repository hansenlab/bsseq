# NOTE: DelayedMatrix objects don't suport subassignment ([<- or [[<-). This
#       function supports a special case of subassignment for DelayedMatrix
#       objects namely:
#       - Create a DelayedMatrix, value, with `nrow` and `ncol`
#       - Fill in value[i, ] with `x`
#       - Fill in remaining rows with `fill`
.subassignRowsDelayedMatrix <- function(x, i, nrow, fill = NA_integer_,
                                        BACKEND = NULL, by_row = FALSE) {
    if (is.null(BACKEND)) {
        if (by_row) {
            warning("'by_row' ignored when 'BACKEND' is NULL")
        }
        z <- matrix(fill, nrow = nrow, ncol = ncol(x))
        z[i, ] <- as.array(x)
        return(realize(z, BACKEND = BACKEND))
    } else {
        # NOTE: Believe by_row == FALSE is better for the situations it's being
        #       called in by bsseq but not done extensive benchmarking
        if (by_row) {
            # Construct the 'missing' rows of the returned matrix in memory,
            # then realize(), then rbind these rows, then reorder the rows,
            # then realize() the returned matrix
            missing_i <- setdiff(seq_len(nrow), i)
            x2 <- .DelayedMatrix(matrix(fill,
                                        nrow = nrow - nrow(x),
                                        ncol = ncol(x)))
            x2 <- realize(x2, BACKEND = BACKEND)
            z <- rbind(x, x2)
            z <- z[order(c(i, missing_i)), ]
            return(realize(z, BACKEND = BACKEND))
        } else {
            # Construct each column of the returned matrix in memory, then
            # realize(), then cbind these columns
            zjs <- lapply(seq_len(ncol(x)), function(j) {
                zj <- matrix(fill, nrow = nrow, ncol = 1)
                zj[i] <- as.array(x[, j])
                realize(zj, BACKEND = BACKEND)
            })
            do.call(cbind, zjs)
        }
    }
}

# X is a list of DelayedMatrix objects
# I is a list of row indices where I[[k]] is the vector of row indices for
# X[[k]]
.combineListOfDelayedMatrixObjects <- function(X, I, nrow, ncol,
                                               dimnames = NULL,
                                               fill = NA_integer_,
                                               BACKEND = NULL) {
    stopifnot(length(X) == length(I))
    is_delayed_matrix <- vapply(X, function(XX) is(XX, "DelayedMatrix"),
                                logical(1L))
    stopifnot(all(is_delayed_matrix))
    is_simple_delayed_matrix <- vapply(X, .isSimpleDelayedMatrix, logical(1L))
    # NOTE: If all DelayedMatrix objects are 'simple' (i.e. the @seed of each
    #       is a base::matrix), then can use this much faster and more
    #       memory-efficient method to combine the DelayedMatrix objects into
    #       a single DelayedMatrix
    if (all(is_simple_delayed_matrix) && is.null(BACKEND)) {
        X <- endoapply(X, as.array)
        z <- matrix(fill, nrow, ncol, dimnames = dimnames)
        j0 <- 0
        for (j in seq_along(X)) {
            i <- I[[j]]
            jj <- j0 + seq_len(ncol(X[[j]]))
            z[i, jj] <- X[[j]]
            j0 <- j0 + ncol(X[[j]])
        }
        return(realize(z, BACKEND = BACKEND))
    } else {
        z <- mapply(.subassignRowsDelayedMatrix, x = X, i = I,
                    MoreArgs = list(nrow = nrow, fill = fill,
                                    BACKEND = BACKEND))
        z <- do.call("cbind", z)
        dimnames(z) <- dimnames
        z
    }
}
