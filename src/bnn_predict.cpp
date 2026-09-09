/**
 * @file bnn_predict.cpp
 * @brief BNN neural-network clock: forward-pass entry point (adapted to R/Rcpp).
 */




// Include Files
#include<Rcpp.h>
#include<vector>

#include "include/rt_nonfinite.h"
#include "include/bnn.h"
#include "include/bnn_predict.h"
#include "include/bnn_terminate.h"
#include "include/bnn_emxAPI.h"
#include "include/bnn_initialize.h"


// Function Declarations
Rcpp::NumericVector bnn_predict(Rcpp::RObject odata );

// Function Definitions

//
// Arguments    :
//    idat : Number of observed data x column
//    odata : data
// Return Type  : void
//
// [[Rcpp::export]]
Rcpp::NumericVector bnn_predict(Rcpp::RObject odata )
{
  
  Rcpp::NumericMatrix data = Rcpp::as<Rcpp::NumericMatrix>(odata);
  
  int icpgs = data.nrow(); // Common CpGs
  int isamples= data.ncol(); // Samples
  std::vector<double> b_y1(isamples);


  // Call the entry-point 'bnn'.
  bnn(data, b_y1.data(), icpgs, isamples);

  return Rcpp::wrap(b_y1);


}
