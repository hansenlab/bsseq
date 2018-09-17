#ifndef BSSEQ_H
#define BSSEQ_H

#include "Rcpp.h"
#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"

// TODO: Remove if not using std::runtime_error()
#include <stdexcept>

// Functions to be called from R.
extern "C" {

    // Validity checking.

    SEXP check_M_and_Cov(SEXP, SEXP);
}

#include "utils.h"

#endif
