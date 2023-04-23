context("BSmooth")

test_that("Errors on bad input", {
    expect_error(
        object = BSmooth(
            BSseq = as(bsseq_test, "SummarizedExperiment"),
            BPPARAM = SerialParam()),
        regexp = "'BSseq' must be a BSseq object")
    expect_error(
        object =  BSmooth(rev(bsseq_test), BPPARAM = SerialParam()),
        regexp = "'BSseq' must be sorted")
    expect_error(
        object = BSmooth(resize(bsseq_test, 2), BPPARAM = SerialParam()),
        regexp = "All loci in 'BSseq' must have width == 1")
})

test_that("BSmooth properly inherits 'dir'", {
    infile <- system.file(
        "extdata", "test_data.fastq_bismark.bismark.cov.gz",
        package = "bsseq")
    bsseq <- read.bismark(
        files = infile,
        BACKEND = "HDF5Array",
        BPPARAM = SerialParam())
    bsseq <- BSmooth(bsseq, BPPARAM = SerialParam())
    expect_is(bsseq, "BSseq")
    expect_true(hasBeenSmoothed(bsseq))
    expect_error(
        BSmooth(bsseq),
        "already contains a dataset named 'coef' or 'se.coef'")
})

test_that("Checks compatability of 'BPPARAM' with the realization backend", {
        expect_error(
            BSmooth(
                BSseq = realize(bsseq_test, "HDF5Array"),
                BPPARAM = SnowParam(c("node1", "node2"))),
            "The parallelisation strategy must use a single machine when using an on-disk realization backend")
    }
)

test_that("Backends work when using SerialParam parallelisation backend", {
    bsseq_NULL <- BSmooth(
        BSseq = realize(bsseq_test, NULL), BPPARAM = SerialParam())
    bsseq_HDF5Array <- BSmooth(bsseq_test, BPPARAM = SerialParam())
    expect_equivalent_SE(bsseq_NULL, bsseq_HDF5Array)
})

test_that(
    "Expected parallelisation backends work with in-memory realization backend", {
        setAutoRealizationBackend(NULL)

        bsseq_serial_param <- BSmooth(bsseq_test, BPPARAM = SerialParam())
        bsseq_multicore_param <- BSmooth(
            bsseq_test,
            BPPARAM = MulticoreParam(2))
        expect_equivalent_SE(bsseq_serial_param, bsseq_multicore_param)

        bsseq_snow_param_single_machine <- BSmooth(
            bsseq_serial_param,
            BPPARAM = SnowParam(2))
        expect_equivalent_SE(bsseq_test, bsseq_snow_param_single_machine)

        bsseq_batch_jobs_param <- BSmooth(
            bsseq_test,
            BPPARAM = BatchtoolsParam(2, progressbar = FALSE))
        expect_equivalent_SE(bsseq_serial_param, bsseq_batch_jobs_param)

        skip_if_not_installed("doParallel")
        doParallel::registerDoParallel(2)
        bsseq_dopar_param <- BSmooth(bsseq_test, BPPARAM = DoparParam())
        expect_equivalent_SE(bsseq_serial_param, bsseq_dopar_param)

        skip("Skipping tests that require SnowParam with multiple machines")
        bsseq_snow_param_multiple_machines <- BSmooth(
            bsseq_test,
            BPPARAM = SnowParam(c("fake1", "fake2")))
        expect_equivalent_SE(
            bsseq_serial_param,
            bsseq_snow_param_multiple_machines)

        skip("Skipping tests that require SnowParam with MPI")
        bsseq_snow_param_mpi <- BSmooth(
            bsseq_test,
            BPPARAM = SnowParam(2, type = "MPI"))
        expect_equivalent_SE(bsseq_serial_param, bsseq_snow_param_mpi)
    }
)

test_that(
    "Expected parallelisation backends work with on-disk realization backend", {
        bsseq_test <- realize(bsseq_test, "HDF5Array")

        bsseq_serial_param <- BSmooth(bsseq_test, BPPARAM = SerialParam())
        bsseq_multicore_param <- BSmooth(
            bsseq_test,
            BPPARAM = MulticoreParam(2))
        expect_equivalent_SE(bsseq_serial_param, bsseq_multicore_param)

        bsseq_snow_param_single_machine <- BSmooth(
            bsseq_test,
            BPPARAM = SnowParam(2))
        expect_equivalent_SE(
            bsseq_serial_param,
            bsseq_snow_param_single_machine)

        # Only single-machine backends support on-disk realization.
        expect_error(
            BSmooth(
                bsseq_serial_param,
                BPPARAM = BatchtoolsParam(2, progressbar = FALSE)))
        expect_error(
            BSmooth(bsseq_test,
                    BPPARAM = SnowParam(workers = c("fake1", "fake2"))))
        expect_error(
            BSmooth(bsseq_test, BPPARAM = SnowParam(2, type = "MPI")))
        skip_if_not_installed("doParallel")
        doParallel::registerDoParallel(2)
        # Only single-machine backends support on-disk realization.
        # TODO: Don't yet have good way to test for single-machine backend when
        #       using DoparaParam; see
        #       https://github.com/Bioconductor/BiocParallel/issues/76
        expect_error(BSmooth(bsseq_test, BPPARAM = DoparParam()))

        skip("Skipping tests that require SnowParam with multiple machines")
        bsseq_snow_param_multiple_machines <- BSmooth(
            bsseq_test,
            BPPARAM = SnowParam(c("fake1", "fake2")))
        expect_equivalent_SE(
            bsseq_serial_param,
            bsseq_snow_param_multiple_machines)

        skip("Skipping tests that require SnowParam with MPI")
        bsseq_snow_param_mpi <- BSmooth(
            bsseq_test,
            BPPARAM = SnowParam(2, type = "MPI"))
        expect_equivalent_SE(bsseq_serial_param, bsseq_snow_param_mpi)
    }
)

test_that("keep.se works", {
    expect_is(
        assay(BSmooth(bsseq_test, keep.se = TRUE, BPPARAM = SerialParam()),
              "se.coef"),
        "matrix")
    expect_error(
        assay(BSmooth(bsseq_test, keep.se = FALSE, BPPARAM = SerialParam()),
              "'se.coef' not in names"))

    expect_is(
        assay(
            BSmooth(
                BSseq = realize(bsseq_test, "HDF5Array"),
                keep.se = TRUE,
                BPPARAM = SerialParam()),
              "se.coef",
            withDimnames = FALSE),
        "HDF5Matrix")
    expect_error(
        assay(BSmooth(bsseq_test, keep.se = FALSE, BPPARAM = SerialParam()),
              "'se.coef' not in names"))
})
