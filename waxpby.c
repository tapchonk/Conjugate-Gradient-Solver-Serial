#include "waxpby.h"

/**
 * @brief Compute the update of a vector with the sum of two scaled vectors
 * 
 * @param n Number of vector elements
 * @param alpha Scalars applied to x
 * @param x Input vector
 * @param beta Scalars applied to y
 * @param y Input vector
 * @param w Output vector
 * @return int 0 if no error
 */
int waxpby (const int n, const double alpha, const double * const x, const double beta, const double * const y, double * const w) {  
  if (alpha==1.0) {
    for (int i=0; i<n; i++) {
      w[i] = x[i] + beta * y[i];
    }
  } else if(beta==1.0) {
    for (int i=0; i<n; i++) {
      w[i] = alpha * x[i] + y[i];
    }
  } else {
    for (int i=0; i<n; i++) {
      w[i] = alpha * x[i] + beta * y[i];
    }
  }

  return 0;
}
