bedMethylCheck
###test
files<-list.files("~/Documents/Likelihood_ms/HG002_PAW70337/6mA",pattern="bed",full.names=TRUE)
f<-data.table::fread(files[1], nrows=100000)
table(f$V6)
sum(f$V18)
sum(f$V5)

##basecaller
files<-list.files("~/Documents/Likelihood_ms/HG002_PAW70337/5mCG_5hmCG",pattern="bed",full.names=TRUE)

bedMethylCheck <- function(file_path) {
    message("Checking input file from the first 100,000 lines of ", file_path)
    f <- data.table::fread(file_path, nrows = 100000)

    # Check that the file has 18 columns
    if (ncol(f) != 18) {
        stop("Error: The bedMethyl file should have 18 columns.")
    }

    #Check the mods in column 4
    mods <- unique(f$V4)

    # STOP if neither "m" nor "h" is present
    if (!("m" %in% mods || "h" %in% mods)) {
        stop("Error: Neither 5mC ('m') nor 5hmC ('h') is present in column 4. Import will not proceed.")
    }

    # Print messages based on the presence of modifications
    if ("m" %in% mods && "h" %in% mods) {
        message("5mC and 5hmC are present and will be imported.")
    } else if ("m" %in% mods) {
        message("5mC is present and will be imported.")
    } else if ("h" %in% mods) {
        message("5hmC is present and will be imported.")
    }

    if ("a" %in% mods && "21839" %in% mods) {
        message("6mA and 4mC are present but will not be imported.")
    } else if ("a" %in% mods) {
        message("6mA is present but will not be imported.")
    } else if ("21839" %in% mods) {
        message("Only 21839 is present but will not be imported.")
    }

    #Check if a CG context modification model was used
    if (any(f$V18 != 0))
        message("CG context modification model detected.")
    if (all(f$V18 == 0))
        message("all context modification model detected.")

    ##check pileup settings
    #check if strand-specific or strand-collapsed pileup was used
    if (all((f$V6 == "+" | f$V6 == "-"))) {
        message("input is stranded.")
    } else if (all(f$V6 == ".")) {
        message("input is strand-collapsed")
    }
    #check if --cpg flag was used
    ff<-f[f$V4 %in% c("m","h"),]
    FP_prop<-sum(ff$V5/(ff$V5+ff$V15+ff$V17+ff$V18)<0.2)/nrow(ff)
    message("Proportion of highly likely non-CpGs: ", round(FP_prop, 4))
}

check_bedMethyl <- function(file, output, n_max = 100000) {
    message("Checking ", file, " as suitable input for ", output, "-object")

    f <- data.table::fread(file, nrows = n_max)

    ## general checks to ensure valid bedMethyl file
    # Check that the file has 18 columns
    if (ncol(f) != 18) {
        message("Error: The bedMethyl file should have 18 columns.")
    }

    # Check the mods in column 4
    mods <- unique(f$V4)

    # Print messages based on the presence of modifications
    if ("m" %in% mods && "h" %in% mods) {
        message("5mC and 5hmC are present and will be imported to the object.")
    } else if ("m" %in% mods) {
        message("5mC is present and will be imported to the object.")
    } else if ("h" %in% mods) {
        message("5hmC is present and will be imported to the object.")
    }

    if ("a" %in% mods && "21839" %in% mods) {
        message("6mA and 4mC are present but will not be imported to the object.")
    } else if ("a" %in% mods) {
        message("6mA is present but will not be imported to the object.")
    } else if ("21839" %in% mods) {
        message("4mC is present but will not be imported to the object.")
    }

    # STOP if neither "m" nor "h" is present
    if (!("m" %in% mods || "h" %in% mods)) {
        stop("Error: No 5mC ('m') or 5hmC ('h') to import; cannot proceed.")
    }

    ff <- f[f$V4 %in% c("m", "h"), ]
    FP_prop <- sum(ff$V5 / (ff$V5 + ff$V15 + ff$V17 + ff$V18) < 0.2, na.rm = TRUE) / max(nrow(ff), 1)

    stranded <- all(f$V6 %in% c("+", "-"))
    collapsed <- all(f$V6 == ".")

    ## Output-specific checks for methylCounts
    if (output == "MethylCounts") {
        message("Checking if suitable for output='MethylCounts'")

        # CG context modification model check
        message("checking modification model..")
        if (any(f$V18 != 0)) {
            message("OK")
        } else {
            message("Fail: The modification calls in ", file, " are called using an all-context model and should be imported as a BSseq object.")
        }

        # Pileup strategy check
        message("checking pileup strategy..")
        if (FP_prop > 0.05 && stranded) {
            message("OK. Input contains ~", round(FP_prop * 100, 2),"% non-CpG which should be filtered, but it seems to includes non-reference CpGs.")
        } else if ((FP_prop <= 0.05 && stranded) || (FP_prop > 0.05 && collapsed) || (FP_prop <= 0.05 && collapsed)) {
            message("Warning: Input contains ~", round(FP_prop * 100, 2),"% non-CpG which can be filtered, however, it seems to be lacking non-reference CpGs.")
        }
    }
    else if (output == "BSseq") {
        message("Checking if suitable for output='BSseq'")

        # Modification model check
        if (any(f$V18 != 0)) {
            message("Modification calls from CG-specific modification model.
                    Input contains ~", round(FP_prop * 100, 2)," non-CpG loci. Remove these using coverage filtering or import file as MethylCounts for likelihood filtering.")

        } else if (all(f$V18 == 0)) {
            message("Modification calls from all-context modification model.
                    checking pileup strategy..")
        }
        message("checking pileup strategy..")
        if (stranded) {
            message("Warning: Input is stranded indicating non-CpGs may be included.")
        } else if (collapsed) {
            message("Input is strand-collapsed indicating it does not contain many non-CpGs.")
        }
        if (FP_prop > 0.02 && stranded) {
            message("Input contains ~", round(FP_prop * 100, 2),"% non-CpG which should be removed using refernce-guided filtering")
        } else if ((FP_prop <= 0.02 && stranded) || (FP_prop > 0.05 && collapsed) || (FP_prop <= 0.05 && collapsed)) {
            message("Input contains ~", round(FP_prop * 100, 2),"% non-CpG. Remove these using coverage filtering. Heterozygous CpGs may remain.")
        }
    } else {
        stop("Error: Invalid output object specified. Use 'MethylCounts' or 'BSseq'.")
    }
}

files<-list.files("~/Documents/Likelihood_ms/HG002_PAW70337",pattern="bed",full.names=TRUE,recursive = T)
files<-list.files("~/Documents/Likelihood_ms/HG002_PAW70337/5mC_5hmC/",pattern="bed",full.names=TRUE,recursive = T)

f<-data.table::fread(files[3], nrows=100000)
ff <- f[f$V4 %in% c("m", "h"), ]
FP_prop <- sum(ff$V5 / (ff$V5 + ff$V15 + ff$V17 + ff$V18) < 0.2) / nrow(ff)


check_bedMethyl(files[1], output = "MethylCounts")
check_bedMethyl(files[2], output = "MethylCounts")
check_bedMethyl(files[3], output = "MethylCounts")
check_bedMethyl(files[4], output = "MethylCounts")
check_bedMethyl(files[5], output = "MethylCounts")
check_bedMethyl(files[6], output = "MethylCounts")
check_bedMethyl(files[7], output = "MethylCounts")
check_bedMethyl(files[8], output = "MethylCounts")
check_bedMethyl(files[9], output = "MethylCounts")
check_bedMethyl(files[10], output = "MethylCounts")
check_bedMethyl(files[11], output = "MethylCounts")
check_bedMethyl(files[12], output = "MethylCounts")
check_bedMethyl(files[13], output = "MethylCounts")
check_bedMethyl(files[14], output = "MethylCounts")
check_bedMethyl(files[15], output = "MethylCounts")
check_bedMethyl(files[16], output = "MethylCounts")
check_bedMethyl(files[17], output = "MethylCounts")
check_bedMethyl(files[18], output = "MethylCounts")
check_bedMethyl(files[19], output = "MethylCounts")
check_bedMethyl(files[20], output = "MethylCounts")

check_bedMethyl(files[1], output = "BSseq")
check_bedMethyl(files[2], output = "BSseq")
check_bedMethyl(files[3], output = "BSseq")
check_bedMethyl(files[4], output = "BSseq")
check_bedMethyl(files[5], output = "BSseq")
check_bedMethyl(files[6], output = "BSseq")
check_bedMethyl(files[7], output = "BSseq")
check_bedMethyl(files[8], output = "BSseq")
check_bedMethyl(files[9], output = "BSseq")
check_bedMethyl(files[10], output = "BSseq")
check_bedMethyl(files[11], output = "BSseq")
check_bedMethyl(files[12], output = "BSseq")
check_bedMethyl(files[13], output = "BSseq")
check_bedMethyl(files[14], output = "BSseq")
check_bedMethyl(files[15], output = "BSseq")
check_bedMethyl(files[16], output = "BSseq")
check_bedMethyl(files[17], output = "BSseq")
check_bedMethyl(files[18], output = "BSseq")
check_bedMethyl(files[19], output = "BSseq")
check_bedMethyl(files[20], output = "BSseq")
#manual checks



f<-data.table::fread(files[4], nrows=100000)
hist((f$V5)/(f$V5+f$V15+f$V17+f$V18))

  print("both 5mC and 5hmC are present")
elif m in mods:
  print("only m present")
elif h in mods:
  print("only h present")
else:
  print("no mods present")

#check which intergers are in column 6
f<-data.table::fread(files[1], nrows=100000)
sum(f$V18)



#check strands in column 6
table(f$V6)
table(f$18)


files<-list.files("~/Documents/Likelihood_ms/HG002_PAW70337/5mCG_5hmCG",pattern="bed",full.names=TRUE)
all<-bsseq:::.readbedMethylAsFWGRanges(files[1], strandCollapse = F, rmZeroCov=T, verbose = T)
all_5mCG_5hmCG<-bsseq:::.readbedMethylAsDT(files[1], verbose = T)
f<-data.table::fread(files[1], nrows=100000)
sum(all_5mCG_5hmCG$NO)

files<-list.files("~/Documents/Likelihood_ms/HG002_PAW70337/5mC_5hmC",pattern="bed",full.names=TRUE)
all_5mC_5hmC<-bsseq:::.readbedMethylAsDT(files[1], verbose = T)
sum(all_5mC_5hmC$NO)


