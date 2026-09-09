/**
 * @file bsxfun.cpp
 * @brief BNN neural-network clock: element-wise broadcast helper.
 */

// Include Files
#include "include/rt_nonfinite.h"
#include "include/bnn.h"
#include "include/bsxfun.h"
#include <Rcpp.h>
// Function Definitions

//
// Arguments    : const double a[1059]
//                const double b[353]
//                double c[1059]
// Return Type  : void
//
void b_bsxfun( double *a,  double *b, double *c, int cpgs, int samples)
{
  int ak = 0;
  int ck;
  int k;
  
  for (ck = 0; ck <= (cpgs*(samples-1))+1; ck += cpgs) {
    for (k = 0; k < cpgs; k++) {
      c[ck + k] = a[ak + k] * b[k];
    }
    
    ak += cpgs;
  }
}

//
// Arguments    : const double a[1059]
//                const double b[353]
//                double c[1059]
// Return Type  : void
//
void bsxfun( double *a,  double *b, double *c, int cpgs, int samples)
{
  int ak = 0;
  int ck;
  int k;
  
  for (ck = 0; ck <= (cpgs*(samples-1))+1; ck += cpgs) {
    for (k = 0; k < cpgs; k++) {
      c[ck + k] = a[ak + k] - b[k];
    }
    
    ak += cpgs;
  }
}

//
// Arguments    : const double a[1059]
//                double c[1059]
// Return Type  : void
//
void c_bsxfun( double *a, double *c, int cpgs, int samples)
{
  int ak = 0;
  int ck;
  int k;
  
  for (ck = 0; ck <= (cpgs*(samples-1))+1; ck += cpgs) {
    for (k = 0; k < cpgs; k++) {
      c[ck + k] = a[ak + k] - 1.0;
    }
    
    ak += cpgs;
  }
}
