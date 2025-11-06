#Likelihood functions for calling methylation status from MethylCounts objects generated from read.bedMethyl

#estimateErrorRate
estimateErrorRate <- function(
        MethylCounts,
        minCov = 10,
        maxCov = 100,
        minRatio = 0.8,
        plotErrorProfile = FALSE) {

    # Input validation of MethylCounts
    if (!inherits(MethylCounts, "MethylCounts")) {
        stop("Error: The input object must be of class 'MethylCounts'.")
    }
    if (ncol(MethylCounts) != 1) {
        stop("Error: The input 'MethylCounts' must have exactly one column.")
    }
    if (!all(c("M", "U", "D") %in% names(assays(MethylCounts)))) {
        stop("Error: The 'MethylCounts' must contain 'M', 'U', and 'D' assays.")
    }

    # Input validation for thresholds
    if (minCov < 0 || maxCov < 0) {
        stop("Error: Both 'minCov' and 'maxCov' must be non-negative.")
    }
    if (maxCov <= minCov) {
        stop("Error: 'maxCov' must be greater than 'minCov'.")
    }
    if (minRatio < 0 || minRatio > 1) {
        stop("Error: 'minRatio' must be between 0 and 1.")
    }

    # Access assay data safely
    cov <- getMethylCounts(MethylCounts, type="Cov")
    d <- getMethylCounts(MethylCounts, type="D")

    # Calculate CpG-to-non-CpG ratio and filter sites
    ratio <- cov / (cov + d)
    sites <- cov >= minCov &
        cov <= maxCov &
        ratio >= minRatio

    # Plot the error profile if requested
    if (plotErrorProfile) {
        if (any(sites)) {
            p_all <- hist(ratio, breaks = 100, plot = FALSE)
            p_used <- hist(ratio[sites], breaks = (1 - minRatio) * 100, plot = FALSE)
            plot(p_all, xlim = c(0, 1), xlab = "CpG to non-CpG ratio", main = "Error Profile")
            plot(p_used, col = "black", xlim = c(0, 1), add = TRUE)
        } else {
            warning("No sites satisfy the filtering criteria. No error profile will be plotted.")
        }
    }

    # Calculate and return the error rate
    error_rate <- sum(d[sites]) / (sum(d[sites]) + sum(cov[sites]))
    return(error_rate)
}

#genotype_likelihood
#internal function to calculate the likelihood of a genotype given the observed data
.genotype_likelihood <- function(n = 2, g, e, cov, D) {
  if (g == 0) {
    p_mismatch <- n * e
    p_match <- n * (1 - e)
  } else if (g == 1) {
    p_mismatch <- n / 2
    p_match <- n / 2
  } else if (g == 2) {
    p_mismatch <- n * (1 - e)
    p_match <- n * e
  } else {
    stop("Invalid genotype value. Must be 0, 1, or 2.")
  }
  likelihood <- (p_mismatch^D) * (p_match^cov)
  return(likelihood)
}


#.getScaledLikelihoods
#internal function to calculate the scaled likelihoods for each genotype
.getScaledLikelihoods <- function(
        MethylCounts,
        e = NULL) {
    # Input validation of BSseq
    if (!inherits(MethylCounts, "MethylCounts")) {
        stop("Error: The input object must be of class 'MethylCounts'.")
    }
    if (ncol(MethylCounts) != 1) {
        stop("Error: The input 'MethylCounts' must have exactly one column.")
    }
    if (!all(c("M", "U", "D") %in% names(assays(MethylCounts)))) {
        stop("Error: The 'MethylCounts' must contain 'M', 'U', and 'D' assays.")
    }
    # Input validation of e
    if (!is.null(e) && (e < 0 || e > 1)) {
        stop("Error: 'e' must be between 0 and 1.")
    }
    if (is.null(e)) {
        e <- estimateErrorRate(MethylCounts)
    }
    cov <- getMethylCounts(MethylCounts, type="Cov")
    d <- getMethylCounts(MethylCounts, type="D")

    l <- matrix(0, nrow = length(cov), ncol = 3)
    l[, 1] <- .genotype_likelihood(g = 0, e = e, cov = cov, D = d)
    l[, 2] <- .genotype_likelihood(g = 1, e = e, cov = cov, D = d)
    l[, 3] <- .genotype_likelihood(g = 2, e = e, cov = cov, D = d)

    SL <- l / rowSums(l)
    return(SL)
}

getCpGs <- function(
    MethylCounts,
    type = c("homozygous", "heterozygous", "allCpG", "nonCpG"),
    threshold = 0.99,
    e = NULL) {
  # Input validation of MethylCounts
  if (!inherits(MethylCounts, "MethylCounts")) {
    stop("Error: The input object must be of class 'MethylCounts'.")
  }
  if (ncol(MethylCounts) != 1) {
    stop("Error: The input 'MethylCounts' must have exactly one column.")
  }
  if (!all(c("M", "U", "D") %in% names(assays(MethylCounts)))) {
    stop("Error: The 'MethylCounts' must contain 'M', 'U', and 'D' assays.")
  }
  # Input validation of type
  if (!type %in% c("homozygous", "heterozygous", "allCpG", "nonCpG")) {
    stop("Error: 'type' must be one of 'homozygous', 'heterozygous', 'allCpG', or 'nonCpG'.")
  }
  # Input validation of threshold
  if (threshold < 0 || threshold > 1) {
    stop("Error: 'threshold' must be between 0 and 1.")
  }
  # Input validation of e
  if (!is.null(e) && (e < 0 || e > 1)) {
    stop("Error: 'e' must be between 0 and 1.")
  }
  # Get scaled likelihoods
  SL <- .getScaledLikelihoods(MethylCounts, e)

  # Determine the CpG sites based on the type
  loci <- switch(type,
                   "homozygous" = which(SL[, 1] > threshold),
                   "heterozygous" = which(SL[, 2] > threshold),
                   "allCpG" = which(rowSums(SL[, 1:2]) > threshold),
                   "nonCpG" = which(SL[, 3] > threshold))

  return(loci)
}

getCpGMatrix <- function(
    MethylCounts,
    e = NULL,
    allCpG = FALSE) {

    # Input validation of MethylCounts
    if (!inherits(MethylCounts, "MethylCounts")) {
        stop("Error: The input object must be of class 'MethylCounts'.")
    }
    if (!all(c("M", "U", "D") %in% names(assays(MethylCounts)))) {
        stop("Error: The 'MethylCounts' must contain 'M', 'U', and 'D' assays.")
    }

  # Input validation of e
  if (!is.null(e)) {
    if (!is.numeric(e) || (length(e) != ncol(MethylCounts))) {
      stop("Error: 'e' must be NULL or a numeric vector with the same length as the number of columns in 'MethylCounts'.")
    }
    if (any(e < 0 | e > 1)) {
      stop("Error: All elements of 'e' must be between 0 and 1.")
    }
  }

  # Initialize CpG matrix
  G <- matrix(nrow = nrow(MethylCounts), ncol = ncol(MethylCounts), byrow = FALSE)

  for (i in 1:ncol(MethylCounts)) {
    # Use the provided error rate or estimate it
    e_i <- if (is.null(e)) estimateErrorRate(MethylCounts[, i]) else e[i]

    # Get scaled likelihoods
    SL_i <- .getScaledLikelihoods(MethylCounts[, i], e_i)

    # Calculate CpG matrix
    if (allCpG) {
      G[, i] <- ifelse(SL_i[, 3] >= 1/3, 2, 0)
    } else {
      G[, i] <- max.col(SL_i, ties.method = "last") - 1
    }
  }

  return(G)
}

getMaxLikelihoodMatrix <- function(
    MethylCounts,
    e = NULL,
    allCpG = FALSE) {

    # Input validation of MethylCounts
    if (!inherits(MethylCounts, "MethylCounts")) {
        stop("Error: The input object must be of class 'MethylCounts'.")
    }
    if (!all(c("M", "U", "D") %in% names(assays(MethylCounts)))) {
        stop("Error: The 'MethylCounts' must contain 'M', 'U', and 'D' assays.")
    }

  # Input validation of e
  if (!is.null(e)) {
    if (!is.numeric(e) || (length(e) != ncol(MethylCounts))) {
      stop("Error: 'e' must be NULL or a numeric vector with the same length as the number of columns in 'MethylCounts'.")
    }
    if (any(e < 0 | e > 1)) {
      stop("Error: All elements of 'e' must be between 0 and 1.")
    }
  }

  # Initialize quality matrix
  Q <- matrix(nrow = nrow(MethylCounts), ncol = ncol(MethylCounts), byrow = FALSE)

  for (i in 1:ncol(MethylCounts)) {
    # Use the provided error rate or estimate it
    e_i <- if (is.null(e)) estimateErrorRate(MethylCounts[, i]) else e[i]

    # Get scaled likelihoods
    SL_i <- .getScaledLikelihoods(MethylCounts[, i], e_i)

    # Calculate quality matrix
    if (allCpG) {
      Q[, i] <- ifelse(SL_i[, 3] >= 1/3, SL_i[, 3], (1 - SL_i[, 3]))
    } else {
      Q[, i] <- DelayedMatrixStats::rowMaxs(SL_i)
    }
  }

  return(Q)
}






