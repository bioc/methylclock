/**
 * @file bnn_initialize.cpp
 * @brief BNN neural-network clock: module initialization.
 */

// Include Files
#include "include/rt_nonfinite.h"
#include "include/bnn.h"
#include "include/bnn_initialize.h"

// Function Definitions

//
// Arguments    : void
// Return Type  : void
//
void bnn_initialize()
{
  rt_InitInfAndNaN(8U);
}
