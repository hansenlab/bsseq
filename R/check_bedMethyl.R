#' Check the input for `read.bedMethyl` compatibility
#'
#' @param file Character. Path to the input bedMethyl file.
#' @param output Character. Type of object to output ("MethylCounts" or "BSseq").
#' @param n_max Integer. Maximum number of rows to read for checking (default: 100000).
#' @return None. Prints messages or stops on errors if input is invalid.
#'
.check_bedMethyl <- function(file, output, n_max = 100000) {
    if (!file.exists(file)) {
        stop("Error: File does not exist: ", file)
    }

    if (!output %in% c("MethylCounts", "BSseq")) {
        stop("Error: Invalid output specified. Use 'MethylCounts' or 'BSseq'.")
    }

    if (!is.numeric(n_max) || n_max <= 0) {
        stop("Error: `n_max` must be a positive integer.")
    }

    message("Checking '", file, "' as suitable input for '", output, "' object...")

    # Read initial rows for validation
    f <- data.table::fread(file, nrows = n_max)

    ## General file structure checks
    if (ncol(f) != 18) {
        stop("Error: The bedMethyl file must have exactly 18 columns.")
    }

    # Check unique modifications in column 4
    modification_types <- unique(f$V4)
    .assess_modifications(modification_types)

    # Stop if essential modifications are missing
    if (!any(modification_types %in% c("m", "h"))) {
        stop("Error: No 5mC ('m') or 5hmC ('h') to import; cannot proceed.")
    }

    # Filter relevant rows and calculate non-CpG proportion
    relevant_rows <- f[f$V4 %in% c("m", "h"), ]
    non_cpg_proportion <- .calculate_non_cpg_proportion(relevant_rows)

    # Check strandedness
    stranded <- all(f$V6 %in% c("+", "-"))
    collapsed <- all(f$V6 == ".")

    if (output == "MethylCounts") {
        .validate_methylcounts_input(f, non_cpg_proportion, stranded, collapsed, file)
    } else {
        .validate_bsseq_input(f, non_cpg_proportion, stranded, collapsed, file)
    }
}

# Helper function to assess modifications
.assess_modifications <- function(modification_types) {
    if ("m" %in% modification_types && "h" %in% modification_types) {
        message("5mC and 5hmC are present and will be imported.")
    } else if ("m" %in% modification_types) {
        message("5mC is present and will be imported.")
    } else if ("h" %in% modification_types) {
        message("5hmC is present and will be imported.")
    }

    if ("a" %in% modification_types) {
        message("6mA is present but will not be imported.")
    }
    if ("21839" %in% modification_types) {
        message("4mC is present but will not be imported.")
    }
}

# Helper function to calculate non-CpG proportion
.calculate_non_cpg_proportion <- function(relevant_rows) {
    sum_low_fp <- sum(relevant_rows$V5 / (relevant_rows$V5 + relevant_rows$V15 +
                                              relevant_rows$V17 + relevant_rows$V18) < 0.2, na.rm = TRUE)
    total_rows <- max(nrow(relevant_rows), 1)
    return(sum_low_fp / total_rows)
}

# Validation for MethylCounts output
.validate_methylcounts_input <- function(f, non_cpg_proportion, stranded, collapsed, file) {
    message("Validating for output = 'MethylCounts'...")

    if (any(f$V18 != 0)) {
        message("Modification calls use CG-context model.")
    } else {
        stop("Error: All-context model detected. Use 'BSseq' instead.")
    }

    if (non_cpg_proportion > 0.05 && stranded) {
        message("Input contains ~", round(non_cpg_proportion * 100, 2),
                "% non-CpG loci. Filtering recommended.")
    } else {
        message("Non-CpG proportion is acceptable.")
    }
}

# Validation for BSseq output
.validate_bsseq_input <- function(f, non_cpg_proportion, stranded, collapsed, file) {
    message("Validating for output = 'BSseq'...")

    if (any(f$V18 != 0)) {
        message("CG-context model detected. Non-CpG loci proportion: ",
                round(non_cpg_proportion * 100, 2), "%.")
    } else {
        message("All-context model detected. Please use coverage filtering.")
    }

    if (stranded) {
        message("Input is stranded; verify if non-CpGs are included.")
    } else if (collapsed) {
        message("Input is strand-collapsed, likely containing fewer non-CpGs.")
    }
}


