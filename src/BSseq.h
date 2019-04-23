#ifndef BSSEQ_H
#define BSSEQ_H

#include "Rcpp.h"

// Functions to be called from R.
extern "C" {

    // Validity checking.

    SEXP check_M_and_Cov(SEXP, SEXP);
}

#endif
