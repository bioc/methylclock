/**
 * @file altumage.cpp
 * @brief AltumAge neural-network clock forward pass in C++/Eigen.
 *
 * The network weights live in HDF5 and are read through the BigDataStatMeth C++
 * API (mc::read_full), which reconciles the R/HDF5 row-column-major layout so the
 * matrices arrive in R-view order. The input is the beta matrix already aligned
 * to the clock's reference CpGs (features x samples) and passed from R; the
 * scaler and every layer weight come from the HDF5 file. Architecture (eval
 * mode): scaler -> [bn_k -> linear_k -> SELU]_{k=1..5} -> bn6 -> linear6.
 */

// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <RcppEigen.h>
#include <BigDataStatMeth.hpp>
#include "methylclock_io.h"
using namespace Rcpp;

/**
 * @brief SELU activation, applied element-wise.
 * @param x Pre-activation matrix.
 * @return SELU(x), same shape as x.
 */
static inline Eigen::MatrixXd mc_selu(const Eigen::MatrixXd& x) {
    const double a = 1.6732632423543772848170429916717;
    const double s = 1.0507009873554804934193349852946;
    return x.unaryExpr([a, s](double v) {
        return v > 0.0 ? s * v : s * a * (std::exp(v) - 1.0);
    });
}

//' AltumAge forward pass (C++/Eigen, weights from HDF5)
//'
//' @param weights_file HDF5 file with the scaler and layer weights.
//' @param Xin Beta matrix aligned to the reference CpGs (features x samples).
//' @return A numeric vector, one predicted age per sample.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector altumage_forward(std::string weights_file,
                                     Rcpp::NumericMatrix Xin) {
    try {
        H5::Exception::dontPrint();
        const int F = Xin.nrow(), N = Xin.ncol();
        Eigen::Map<Eigen::MatrixXd> X(Xin.begin(), F, N);   // features x samples

        // scaler: (x - center) / (scale + 1e-8), per feature (row).
        Eigen::VectorXd center = mc::read_vec(weights_file, "scaler", "center");
        Eigen::VectorXd scale  = mc::read_vec(weights_file, "scaler", "scale");
        Eigen::MatrixXd h = X.colwise() - center;
        Eigen::VectorXd inv = (scale.array() + 1e-8).inverse();
        h = h.array().colwise() * inv.array();

        for (int k = 1; k <= 6; ++k) {
            const std::string b = "bn" + std::to_string(k);
            Eigen::VectorXd g  = mc::read_vec(weights_file, "weights", b + ".weight");
            Eigen::VectorXd be = mc::read_vec(weights_file, "weights", b + ".bias");
            Eigen::VectorXd m  = mc::read_vec(weights_file, "weights", b + ".running_mean");
            Eigen::VectorXd v  = mc::read_vec(weights_file, "weights", b + ".running_var");
            // batchnorm (eval): (h - m)/sqrt(v + eps) * g + be, per feature.
            Eigen::VectorXd sd = (v.array() + 0.001).sqrt();
            h = h.colwise() - m;
            h = h.array().colwise() * (g.array() / sd.array());
            h = h.colwise() + be;

            const std::string l = "linear" + std::to_string(k);
            Eigen::MatrixXd W  = mc::read_full(weights_file, "weights", l + ".weight");
            Eigen::VectorXd lb = mc::read_vec(weights_file, "weights", l + ".bias");
            h = (W * h).colwise() + lb;                     // out x samples
            if (k < 6) h = mc_selu(h);
        }
        // h is 1 x N: one age per sample.
        Rcpp::NumericVector out(N);
        for (int j = 0; j < N; ++j) out[j] = h(0, j);
        return out;
    } catch (H5::Exception& e) {
        Rf_error("altumage_forward HDF5 error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("altumage_forward: %s", e.what());
    }
    return R_NilValue;
}
