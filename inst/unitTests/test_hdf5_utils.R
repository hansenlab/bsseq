test_.isHDF5ArrayBacked <- function() {
    matrix <- matrix(1:10, ncol = 2)
    da <- bsseq:::.DelayedMatrix(matrix)
    ha <- realize(da, BACKEND = "HDF5Array")
    ha_ha <- rbind(ha, ha)
    da_da <- rbind(da, da)
    ha_ha_ha_ha <- rbind(ha_ha, ha_ha)
    da_da_da_da <- rbind(da_da, da_da)
    da_da_ha_ha <- rbind(da_da, ha_ha)
    da_with_da_seed <- bsseq:::.collapseDelayedMatrix(x = da,
                                                     sp = list(1:2, 3:5),
                                                     MARGIN = 2)
    ha_with_ha_seed <- bsseq:::.collapseDelayedMatrix(
        x = ha,
        sp = list(1:2, 3:5),
        MARGIN = 2,
        BACKEND = "HDF5Array")

    checkTrue(!bsseq:::.isHDF5ArrayBacked(matrix))
    checkTrue(bsseq:::.isHDF5ArrayBacked(ha))
    checkTrue(!bsseq:::.isHDF5ArrayBacked(da))
    checkTrue(bsseq:::.isHDF5ArrayBacked(ha_ha))
    checkTrue(!bsseq:::.isHDF5ArrayBacked(da_da))
    checkTrue(bsseq:::.isHDF5ArrayBacked(ha_ha_ha_ha))
    checkTrue(!bsseq:::.isHDF5ArrayBacked(da_da_da_da))
    checkTrue(bsseq:::.isHDF5ArrayBacked(da_da_ha_ha))

    # A complicated DelayedArray with another DelayedArray as a seed
    checkTrue(!bsseq:::.isHDF5ArrayBacked(da_with_da_seed))
    checkTrue(bsseq:::.isHDF5ArrayBacked(ha_with_ha_seed))
}
