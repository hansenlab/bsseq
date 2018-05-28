#include "BSseq.h"

// Special function to check for NA'ness.

bool isNA(int x) {
    return x==NA_INTEGER;
}

bool isNA(double x) {
    return ISNAN(x);
}
