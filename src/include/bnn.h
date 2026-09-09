/**
 * @file bnn.h
 * @brief BNN neural-network clock: core network.
 */

#ifndef BNN_H
#define BNN_H

// Include Files
#include <Rcpp.h>
#include <cmath>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "bnn_types.h"

// Function Declarations
extern void bnn(Rcpp::NumericMatrix x1, double b_y1[], int cpgs, int samples);

#endif
