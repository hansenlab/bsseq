## The following script downloads and constructs the BS.chr22 dataset,
## included in the bsseq package

library(bsseq)

## First we download.  Each file is slightly less than 200 MB

download.file(
    url = "ftp://ftpuser3:s3qu3nc3@neomorph.salk.edu/mc/mc_imr90_r1.tar.gz",
    destfile = "mc_imr90_r1.tar.gz")
untar("mc_imr90_r1.tar.gz",  "mc_imr90_r1/mc_imr90_r1_22", compressed = TRUE)
download.file(
    url = "ftp://ftpuser3:s3qu3nc3@neomorph.salk.edu/mc/mc_imr90_r2.tar.gz",
    destfile = "mc_imr90_r2.tar.gz")
untar("mc_imr90_r2.tar.gz",  "mc_imr90_r2/mc_imr90_r2_22", compressed = TRUE)

## Now the workhorse function

read.lister <- function(file) {
    dat <- read.table(
        file,
        skip = 1,
        row.names = NULL,
        col.names = c("chr", "pos", "strand", "context", "M", "Cov"),
        colClasses = c("character", "integer", "character", "character",
                       "integer", "integer"))
    ## we remove all non-CpG calls.  This includes SNPs
    dat <- dat[dat$context == "CG", ]
    dat$context <- NULL
    dat$chr <- paste("chr", dat$chr, sep = "")
    ## Now we need to handle that the data has separate lines for each strand
    ## We join these
    tmp <- dat[dat$strand == "+",]
    BS.forward <- BSseq(
        pos = tmp$pos,
        chr = tmp$chr,
        M = as.matrix(tmp$M, ncol = 1),
        Cov = as.matrix(tmp$Cov, ncol = 1),
        sampleNames = "forward")
    tmp <- dat[dat$strand == "-",]
    BS.reverse <- BSseq(
        pos = tmp$pos - 1L,
        chr = tmp$chr,
        M = as.matrix(tmp$M, ncol = 1),
        Cov = as.matrix(tmp$Cov, ncol = 1),
        sampleNames = "reverse")
    BS <- combine(BS.forward, BS.reverse)
    BS <- collapseBSseq(BS, columns = c("a", "a"))
    BS
}

BS.r1 <- read.lister("mc_imr90_r1/mc_imr90_r1_22")
sampleNames(BS.r1) <- "r1"
BS.r2 <- read.lister("mc_imr90_r2/mc_imr90_r2_22")
sampleNames(BS.r2) <- "r2"

BS.chr22 <- combine(BS.r1, BS.r2)
pData(BS.chr22)$Rep <- c("replicate1", "replicate2")
validObject(BS.chr22)
pData(BS.chr22)

# NOTE: To reduce size of object, set the storage.mode() of assays to integer.
#       (at this point they use doubles because collapseBSseq() returns
#       doubles).
storage.mode(assay(BS.chr22, 1)) <- "integer"
storage.mode(assay(BS.chr22, 2)) <- "integer"

save(BS.chr22, file = "BS.chr22.rda")
library(tools)
resaveRdaFiles("BS.chr22.rda")
