context("MethylCounts")

test_that("Likelihood functions works for bedMethyl files", {

infiles <- c(system.file("extdata/HG002_nanopore_test.bedMethyl.gz",
                                   package = "bsseq"),
             system.file("extdata/HG002_pacbio_test.bedMethyl.gz",
                                   package = "bsseq"))

mc <- read.bedMethyl(files = infiles,
                        colData = DataFrame(row.names = c("test_nanopore","test_pacbio")),
                        rmZeroCov = TRUE,
                        strandCollapse = TRUE,
                        output = "MethylCounts",
                        verbose = FALSE)

# Filter CpG sites for the Nanopore dataset and test number of loci
mc_nano <- mc[, 1]
mc_nano_99All_filtered <- mc_nano[getCpGs(mc_nano,
                                           type = "allCpG", threshold = 0.99)]
expect_equal(dim(mc_nano_99All_filtered), c(828L, 1L))

# Filter CpG sites for the PacBio dataset and test number of loci
mc_pacbio <- mc[, 2]
mc_pacbio_99All_filtered <- mc_pacbio[getCpGs(mc_pacbio,
                                             type = "allCpG", threshold = 0.99)]
expect_equal(dim(mc_pacbio_99All_filtered), c(827L, 1L))

CpGMatrix <- getCpGMatrix(mc, allCpG = TRUE)
MaxLikelihoodMatrix <- getMaxLikelihoodMatrix(mc, allCpG = TRUE)

# Filter for allCpG loci with a likelihood > 0.99 in both samples
mc_combined_99All_filtered <- mc[which(rowAlls(CpGMatrix == 0)
                                       & rowMins(MaxLikelihoodMatrix) > 0.99)]
expect_equal(dim(mc_combined_99All_filtered), c(827L, 2L))
})

