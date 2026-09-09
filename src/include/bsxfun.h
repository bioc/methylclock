/**
 * @file bsxfun.h
 * @brief BNN neural-network clock: element-wise broadcast helper.
 */

#ifndef BSXFUN_H
#define BSXFUN_H

// Include Files
#include <cmath>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "bnn_types.h"

// Function Declarations
void b_bsxfun( double *a,  double *b, double *c, int cpgs, int samples);
void bsxfun( double *a,  double *b, double *c, int cpgs, int samples);
void c_bsxfun( double *a, double *c, int cpgs, int samples);

#endif
