/**
 * @file bnn_emxutil.h
 * @brief BNN neural-network clock: dynamic-array (emx) utilities.
 */

#ifndef BNN_EMXUTIL_H
#define BNN_EMXUTIL_H

// Include Files
#include <cmath>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "bnn_types.h"

// Function Declarations
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);

#endif
