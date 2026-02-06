# Helper function to calculate the non-CpG proportion from relevant rows
.calculate_non_cpg_proportion <- function(relevant_rows) {
    # Ensure that the necessary columns are available and calculate proportions
    if (!all(c("V5", "V15", "V17", "V18") %in% colnames(relevant_rows))) {
        stop("Error: Missing necessary columns (V5, V15, V17, V18) for non-CpG proportion calculation.")
    }

    # Calculate the low-fraction CpG scores and sum them
    sum_low_fp <- sum(
        relevant_rows$V5 / (relevant_rows$V5 + relevant_rows$V15 +
                                relevant_rows$V17 + relevant_rows$V18) < 0.2,
        na.rm = TRUE
    )

    # Calculate the total count of rows (to handle edge cases where nrow = 0)
    total_rows <- max(nrow(relevant_rows), 1)

    # Return the non-CpG proportion as the ratio of low-fraction counts to total rows
    proportion <- sum_low_fp / total_rows
    return(proportion)
}

.assess_modifications <- function(modification_types) {
    # Logical checks for the presence of modification types
    detected_mods <- c(
        `5mC` = "m" %in% modification_types,
        `5hmC` = "h" %in% modification_types,
        `4mC` = "21839" %in% modification_types,
        `6mA` = "a" %in% modification_types
    )
    return(detected_mods)
}

.check_bedMethyls <- function(files, output, n_max = 100000) {
    # Initialize an empty data.frame to store results
    info_df <- data.frame(
        File = basename(files),
        Mod_5mC = logical(length(files)),
        Mod_5hmC = logical(length(files)),
        Mod_4mC = logical(length(files)),
        Mod_6mA = logical(length(files)),
        CGContext = logical(length(files)),
        Stranded = logical(length(files)),
        NonCpGs = numeric(length(files)),
        stringsAsFactors = FALSE
    )

    # Iterate over all files and collect the information
    for (i in seq_along(files)) {
        file <- files[i]
        if (!file.exists(file)) {
            stop("File does not exist: ", file)
        }

        # Read a subset of the file for validation
        f <- data.table::fread(file, nrows = n_max)

        # General structure checks
        if (ncol(f) != 18) {
            stop("Error: The bedMethyl file must have exactly 18 columns.")
        }

        # Assess modifications (logical vector)
        mod_df <- .assess_modifications(unique(f$V4))

        # Calculate the non-CpG proportion
        relevant_rows <- f[f$V4 %in% c("m", "h"), ]
        non_cpg_proportion <- .calculate_non_cpg_proportion(relevant_rows)

        # Check strandedness and CG-context details
        stranded <- all(f$V6 %in% c("+", "-"))
        cg_context <- any(f$V18 != 0)

        # Write results into the data.frame
        info_df[i, "Mod_5mC"] <- mod_df["5mC"]
        info_df[i, "Mod_5hmC"] <- mod_df["5hmC"]
        info_df[i, "Mod_4mC"] <- mod_df["4mC"]
        info_df[i, "Mod_6mA"] <- mod_df["6mA"]
        info_df[i, "CGContext"] <- cg_context
        info_df[i, "Stranded"] <- stranded
        info_df[i, "NonCpGs"] <- round(non_cpg_proportion, 4)  # Limit to 4 decimals
    }

    return(info_df)
}
