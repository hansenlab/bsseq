context("BSmooth")

test_that("Errors on bad input", {
    expect_error(
        BSmooth(as(bsseq_test, "SummarizedExperiment")),
        "'BSseq' must be a BSseq object")
    expect_error(BSmooth(rev(bsseq_test)), "'BSseq' must be sorted")
    expect_error(
        BSmooth(resize(bsseq_test, 2)),
        "All loci in 'BSseq' must have width == 1")
})

test_that("Issues warning for deprecated arguments", {
    expect_warning(
        BSmooth(bsseq_test, parallelBy = "sample", BPPARAM = SerialParam()),
        "'parallelBy' is deprecated")
    expect_warning(
        BSmooth(bsseq_test, mc.preschedule = FALSE, BPPARAM = SerialParam()),
        "'mc.preschedule' is deprecated")
    expect_warning(
        BSmooth(bsseq_test, mc.cores = 1, BPPARAM = SerialParam()),
        "'mc.cores' is deprecated")
    expect_warning(
        BSmooth(bsseq_test, verbose = FALSE, BPPARAM = SerialParam()),
        "'verbose' is deprecated")
})

test_that(
    "Checks compatability of realization backend with backend(s) of BSseq object", {
        setRealizationBackend(NULL)
        expect_error(
            BSmooth(realize(bsseq_test, "HDF5Array"),
                    "Using an in-memory backend for a disk-backed BSseq object is not supported")
        )
        bsseq_test2 <- bsseq_test
        assay(bsseq_test2, withDimnames = FALSE) <- realize(
            assay(bsseq_test2, withDimnames = FALSE),
            "HDF5Array")
        expect_error(
            BSmooth(realize(bsseq_test2, "HDF5Array"),
                    "Using an in-memory backend for a disk-backed BSseq object is not supported")
        )

        setRealizationBackend("RleArray")
        expect_error(
            BSmooth(realize(bsseq_test, "HDF5Array"),
                    "Using an in-memory backend for a disk-backed BSseq object is not supported")
        )
        setRealizationBackend("HDF5Array")
        expect_silent(
            BSmooth(realize(bsseq_test, "HDF5Array"), BPPARAM = SerialParam()))
    })

test_that("Checks compatability of 'BPPARAM' with the realization backend", {
        setRealizationBackend("HDF5Array")
        expect_error(
            BSmooth(bsseq_test, BPPARAM = SnowParam(c("node1", "node2"))),
            "The parallelisation strategy must use a single machine when using an on-disk realization backend")
    }
)

test_that(
    "Expected realization backends work when using SerialParam parallelisation backend", {
        setRealizationBackend(NULL)
        bsseq_NULL <- BSmooth(bsseq_test, BPPARAM = SerialParam())

        setRealizationBackend("HDF5Array")
        bsseq_HDF5Array <- BSmooth(bsseq_test, BPPARAM = SerialParam())
        expect_equivalent_SE(bsseq_NULL, bsseq_HDF5Array)

        setRealizationBackend("RleArray")
        expect_error(
            BSmooth(bsseq_test, BPPARAM = SerialParam()),
            "The 'RleArray' realization backend is not supported")
})

test_that(
    "Expected parallelisation backends work with in-memory realization backend", {
        setRealizationBackend(NULL)

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
            BPPARAM = BatchJobsParam(2, progressbar = FALSE))
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
        setRealizationBackend("HDF5Array")

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
                BPPARAM = BatchJobsParam(2, progressbar = FALSE)))
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
    setRealizationBackend(NULL)
    expect_is(
        assay(BSmooth(bsseq_test, keep.se = TRUE, BPPARAM = SerialParam()),
              "se.coef"),
        "matrix")
    expect_error(
        assay(BSmooth(bsseq_test, keep.se = FALSE, BPPARAM = SerialParam()),
              "'se.coef' not in names"))

    setRealizationBackend("HDF5Array")
    expect_is(
        assay(BSmooth(bsseq_test, keep.se = TRUE, BPPARAM = SerialParam()),
              "se.coef", withDimnames = FALSE),
        "HDF5Matrix")
    expect_error(
        assay(BSmooth(bsseq_test, keep.se = FALSE, BPPARAM = SerialParam()),
              "'se.coef' not in names"))
})
