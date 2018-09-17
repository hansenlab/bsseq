#include "BSseq.h"

// NOTE: Returning Rcpp::CharacterVector rather than throwing an error because
//       this function is used within a validity method.

template <class M_column_class, class Cov_column_class, class M_class,
          class Cov_class>
Rcpp::RObject check_M_and_Cov_internal(M_class M_bm, Cov_class Cov_bm) {
    BEGIN_RCPP

    // Get the dimensions of 'M' and 'Cov' and check these are compatible.
    const size_t M_nrow = M_bm->get_nrow();
    const size_t Cov_nrow = Cov_bm->get_nrow();
    if (M_nrow != Cov_nrow) {
        return Rcpp::CharacterVector(
            "'M' and 'Cov' must have the same number of rows.");
    }
    const size_t M_ncol = M_bm->get_ncol();
    const size_t Cov_ncol = Cov_bm->get_ncol();
    if (M_ncol != Cov_ncol) {
        return Rcpp::CharacterVector(
            "'M' and 'Cov' must have the same number of columns.");
    }

    // Simultaneously loop over columns of 'M' and 'Cov', checking that
    // `all(0 <= M <= Cov) && !anyNA(M) && !anyNA(Cov)` && all(is.finite(Cov)).
    M_column_class M_column(M_nrow);
    Cov_column_class Cov_column(Cov_nrow);
    for (size_t j = 0; j < M_ncol; ++j) {
        // Copy the j-th column of M to M_column and the j-th column of Cov to
        // Cov_column
        M_bm->get_col(j, M_column.begin());
        Cov_bm->get_col(j, Cov_column.begin());
        // Construct iterators
        // NOTE: Iterators constructed outside of loop because they may be of
        //       different type, which is not supported within a for loop
        //       constructor.
        auto M_column_it = M_column.begin();
        auto Cov_column_it = Cov_column.begin();
        for (M_column_it = M_column.begin(), Cov_column_it = Cov_column.begin();
             M_column_it != M_column.end();
             ++M_column_it, ++Cov_column_it) {
            if (isNA(*M_column_it)) {
                return Rcpp::CharacterVector("'M' must not contain NAs.");
            }
            if (isNA(*Cov_column_it)) {
                return Rcpp::CharacterVector("'Cov' must not contain NAs.");
            }
            if (*M_column_it < 0) {
                return Rcpp::CharacterVector(
                    "'M' must not contain negative values.");
            }
            if (*M_column_it > *Cov_column_it) {
                return Rcpp::CharacterVector(
                    "All values of 'M' must be less than or equal to the corresponding value of 'Cov'.");
            }
            if (!R_FINITE(*Cov_column_it)) {
                return Rcpp::CharacterVector("All values of 'Cov' must be finite.");
            }
        }
    }

    return R_NilValue;
    END_RCPP
}

SEXP check_M_and_Cov(SEXP M, SEXP Cov) {
    BEGIN_RCPP

    // Get the type of 'M' and 'Cov',
    int M_type = beachmat::find_sexp_type(M);
    int Cov_type = beachmat::find_sexp_type(Cov);
    if (M_type == INTSXP && Cov_type == INTSXP) {
        auto M_bm = beachmat::create_integer_matrix(M);
        auto Cov_bm = beachmat::create_integer_matrix(Cov);
        return check_M_and_Cov_internal<
            Rcpp::IntegerVector, Rcpp::IntegerVector>(M_bm.get(), Cov_bm.get());
    } else if (M_type == REALSXP && Cov_type == REALSXP) {
        auto M_bm = beachmat::create_numeric_matrix(M);
        auto Cov_bm = beachmat::create_numeric_matrix(Cov);
        return check_M_and_Cov_internal<
            Rcpp::NumericVector, Rcpp::NumericVector>(M_bm.get(), Cov_bm.get());
    } else if (M_type == INTSXP && Cov_type == REALSXP) {
        auto M_bm = beachmat::create_integer_matrix(M);
        auto Cov_bm = beachmat::create_numeric_matrix(Cov);
        return check_M_and_Cov_internal<
            Rcpp::IntegerVector, Rcpp::NumericVector>(M_bm.get(), Cov_bm.get());
    } else if (M_type == REALSXP && Cov_type == INTSXP) {
        auto M_bm = beachmat::create_numeric_matrix(M);
        auto Cov_bm = beachmat::create_integer_matrix(Cov);
        return check_M_and_Cov_internal<
            Rcpp::NumericVector, Rcpp::IntegerVector>(M_bm.get(), Cov_bm.get());
    }
    else {
        return Rcpp::CharacterVector(
            "'M' and 'Cov' must contain integer or numeric values.");
    }
    END_RCPP
}
// TODOs -----------------------------------------------------------------------

// TODO: Add code path to process ordinary R vectors (for use within
//       read.bismark() funcionality)?
// TODO: Wasn't able to figure out how to use get_const_col. See
//       https://gist.github.com/PeteHaitch/cf6cffd40d4f8bd7082eec1b0f330082
//       and note that `check_M_and_Cov_using_const_columns` never fails when
//       matrix input contains NA
