#include "BSseq.h"

#include "Rtatami.h"

#include "utils.h"

#include <string>
#include <vector>

// NOTE: Returning Rcpp::CharacterVector rather than throwing an error because
//       this function is used within a validity method.

SEXP check_M_and_Cov(SEXP M, SEXP Cov, SEXP nt) {
    BEGIN_RCPP

    Rtatami::BoundNumericPointer M_bound(M);
    const auto& M_bm = *(M_bound->ptr);
    Rtatami::BoundNumericPointer Cov_bound(Cov);
    const auto& Cov_bm = *(Cov_bound->ptr);

    // Get the dimensions of 'M' and 'Cov' and check these are compatible.
    const int M_nrow = M_bm.nrow();
    const int Cov_nrow = Cov_bm.nrow();
    if (M_nrow != Cov_nrow) {
        return Rcpp::CharacterVector(
            "'M' and 'Cov' must have the same number of rows.");
    }
    const int M_ncol = M_bm.ncol();
    const int Cov_ncol = Cov_bm.ncol();
    if (M_ncol != Cov_ncol) {
        return Rcpp::CharacterVector(
            "'M' and 'Cov' must have the same number of columns.");
    }

    Rcpp::IntegerVector raw_nt(nt);
    if (raw_nt.size() != 1 || raw_nt[0] <= 0) {
        return Rcpp::CharacterVector(
            "Number of threads should be a positive integer.");
    }
    int nthreads = raw_nt[0];

    // Simultaneously loop over columns of 'M' and 'Cov', checking that
    // `all(0 <= M <= Cov) && !anyNA(M) && !anyNA(Cov)` && all(is.finite(Cov)).
    std::vector<std::string> errors(nthreads);
    tatami::parallelize([&](int tid, int start, int length) {
        std::vector<double> M_buffer(M_nrow), Cov_buffer(Cov_nrow);
        auto M_ext = tatami::consecutive_extractor<false>(&M_bm, false, start, length);
        auto Cov_ext = tatami::consecutive_extractor<false>(&Cov_bm, false, start, length);

        for (int c = start, cend = start + length; c < cend; ++c) {
            auto M_ptr = M_ext->fetch(M_buffer.data());
            auto Cov_ptr = Cov_ext->fetch(Cov_buffer.data());

            for (int r = 0; r < M_nrow; ++r) {
                auto M_current = M_ptr[r];
                auto Cov_current = Cov_ptr[r];

                if (isNA(M_current)) {
                    errors[tid] = "'M' must not contain NAs.";
                    return;
                }
                if (isNA(Cov_current)) {
                    errors[tid] = "'Cov' must not contain NAs.";
                    return;
                }
                if (M_current < 0) {
                    errors[tid] = "'M' must not contain negative values.";
                    return;
                }
                if (M_current > Cov_current) {
                    errors[tid] = "All values of 'M' must be less than or equal to the corresponding value of 'Cov'.";
                    return;
                }
                if (!R_FINITE(Cov_current)) {
                    errors[tid] = "All values of 'Cov' must be finite.";
                    return;
                }
            }
        }
    }, M_ncol, nthreads);

    for (const auto& msg : errors) {
        if (!msg.empty()) {
            return Rcpp::CharacterVector(msg.c_str());
        }
    }

    return R_NilValue;
    END_RCPP
}

// TODOs -----------------------------------------------------------------------

// TODO: Add code path to process ordinary R vectors (for use within
//       read.bismark() funcionality)?
// TODO: Wasn't able to figure out how to use get_const_col. See
//       https://gist.github.com/PeteHaitch/cf6cffd40d4f8bd7082eec1b0f330082
//       and note that `check_M_and_Cov_using_const_columns` never fails when
//       matrix input contains NA
