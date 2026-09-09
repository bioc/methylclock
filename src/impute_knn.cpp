/**
 * @file impute_knn.cpp
 * @brief Out-of-core K-nearest-neighbour imputation of a methylation matrix,
 *        built on the BigDataStatMeth C++ algebra so the heavy linear algebra
 *        (block SVD, cross-products) and the HDF5 I/O run out of core and
 *        reconcile the R/HDF5 row/column-major layout through the library, never
 *        by hand.
 *
 * The matrix is CpGs (rows) by samples (columns). A missing value in row i,
 * column j is filled from the k CpGs most similar to i -- by averaged Euclidean
 * distance over their co-observed samples -- among those observed at column j,
 * as the unweighted mean of their values there. This matches the reference
 * gene-expression KNN imputer (Troyanskaya 2001), verified cell-by-cell in the
 * package tests. The file provides the in-memory exact kernel first; the block
 * split and the out-of-core HDF5 path are added on top of it.
 */

// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <RcppEigen.h>
#include <BigDataStatMeth.hpp>
#include "methylclock_io.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
using namespace Rcpp;

namespace mcknn {

/**
 * @brief Averaged Euclidean distance between two rows over their co-observed
 *        columns: sqrt(mean over co-observed (a-b)^2). Returns +inf when they
 *        share no observed column, so such a pair is never chosen as a
 *        neighbour. Missing entries are NaN. Matching the reference imputer, the
 *        k nearest rows are chosen once per target row by this distance
 *        (globally, not per missing column).
 * @param B CpGs x samples matrix (missing entries are NaN).
 * @param i Index of the target row.
 * @param g Index of the candidate neighbour row.
 * @return The averaged co-observed distance, or +inf if no column is shared.
 */
static inline double row_dist(const Eigen::MatrixXd& B, int i, int g) {
    double ss = 0.0;
    int nobs = 0;
    const int n = static_cast<int>(B.cols());
    for (int c = 0; c < n; ++c) {
        const double a = B(i, c), b = B(g, c);
        if (!std::isnan(a) && !std::isnan(b)) {
            const double d = a - b;
            ss += d * d;
            ++nobs;
        }
    }
    if (nobs == 0) return std::numeric_limits<double>::infinity();
    return std::sqrt(ss / static_cast<double>(nobs));
}

/**
 * @brief Per-column mean over observed entries, used as the fallback fill for a
 *        cell with no usable neighbour, and for rows too sparse to trust for KNN.
 * @param B CpGs x samples matrix (missing entries are NaN).
 * @return A length-ncols vector of per-column means over observed entries.
 */
static inline Eigen::VectorXd col_means(const Eigen::MatrixXd& B) {
    const int m = static_cast<int>(B.rows()), n = static_cast<int>(B.cols());
    Eigen::VectorXd mu = Eigen::VectorXd::Zero(n);
    for (int c = 0; c < n; ++c) {
        double s = 0.0;
        int cnt = 0;
        for (int r = 0; r < m; ++r) {
            const double v = B(r, c);
            if (!std::isnan(v)) { s += v; ++cnt; }
        }
        mu(c) = cnt ? s / cnt : 0.0;
    }
    return mu;
}

/**
 * @brief Exact in-memory KNN imputation of a block, in place. Mirrors the
 *        reference imputer cell-for-cell: rows with more than trunc(rowmax*n)
 *        missing entries are set aside and filled with the per-column mean (over
 *        the kept rows) and are not used as neighbours; every other row with gaps
 *        takes the k globally nearest kept rows (by averaged Euclidean distance,
 *        chosen once), and each missing cell is the mean of those neighbours
 *        observed at that column (falling back to the column mean when none are).
 * @param B CpGs x samples block; imputed in place (missing entries are NaN).
 * @param k Number of neighbours averaged for each imputed value.
 * @param rowmax A row with a greater missing fraction is filled with the column
 *        mean instead of by KNN, and excluded as a neighbour.
 * @param colmax A column with a greater missing fraction is counted for the
 *        caller (a warning in R); imputation still proceeds.
 * @return The count of columns exceeding colmax (0 = clean).
 */
static int knn_impute_block(Eigen::MatrixXd& B, int k, double rowmax,
                            double colmax) {
    const int m = static_cast<int>(B.rows()), n = static_cast<int>(B.cols());
    if (m == 0 || n == 0) return 0;

    // Read everything (distances and neighbour values) from the original data;
    // write imputed values into B. Keeping the source read-only stops a just-
    // imputed cell from feeding into another row's neighbours -- every distance
    // and average is over the observed input, as in the reference imputer.
    const Eigen::MatrixXd A = B;

    Eigen::VectorXi row_na = Eigen::VectorXi::Zero(m);
    Eigen::VectorXi col_na = Eigen::VectorXi::Zero(n);
    for (int r = 0; r < m; ++r)
        for (int c = 0; c < n; ++c)
            if (std::isnan(A(r, c))) { ++row_na(r); ++col_na(c); }
    int overcol = 0;
    for (int c = 0; c < n; ++c)
        if (static_cast<double>(col_na(c)) / m > colmax) ++overcol;

    // Kept rows (usable as neighbours): missing count within rowmax; the rest are
    // too sparse to trust for KNN (matches knnimp's row split).
    const int imax = static_cast<int>(std::trunc(rowmax * n));
    std::vector<char> keep(m);
    for (int r = 0; r < m; ++r) keep[r] = (row_na(r) <= imax);

    // Per-column mean over the kept rows' observed entries (fallback fill).
    Eigen::VectorXd mu = Eigen::VectorXd::Zero(n);
    {
        Eigen::VectorXi cnt = Eigen::VectorXi::Zero(n);
        for (int r = 0; r < m; ++r) {
            if (!keep[r]) continue;
            for (int c = 0; c < n; ++c)
                if (!std::isnan(A(r, c))) { mu(c) += A(r, c); ++cnt(c); }
        }
        for (int c = 0; c < n; ++c) if (cnt(c)) mu(c) /= cnt(c);
    }

    for (int r = 0; r < m; ++r) {
        if (row_na(r) == 0) continue;
        if (!keep[r]) {                          // too sparse: per-column mean
            for (int c = 0; c < n; ++c)
                if (std::isnan(A(r, c))) B(r, c) = mu(c);
            continue;
        }
        // The k nearest kept rows, chosen once for this row (globally, not per
        // column). Distance to self is +inf so it is never selected.
        std::vector<std::pair<double, int>> nbr;
        nbr.reserve(m);
        for (int g = 0; g < m; ++g) {
            if (g == r || !keep[g]) continue;
            const double dd = row_dist(A, r, g);
            if (std::isfinite(dd)) nbr.emplace_back(dd, g);
        }
        const int kk = std::min<int>(k, static_cast<int>(nbr.size()));
        std::partial_sort(nbr.begin(), nbr.begin() + kk, nbr.end());
        for (int c = 0; c < n; ++c) {
            if (!std::isnan(A(r, c))) continue;
            double s = 0.0;
            int cnt = 0;
            for (int t = 0; t < kk; ++t) {         // average observed neighbours
                const double v = A(nbr[t].second, c);
                if (!std::isnan(v)) { s += v; ++cnt; }
            }
            B(r, c) = cnt ? s / cnt : mu(c);
        }
    }
    return overcol;
}

// --- block split by SVD embedding (for row sets larger than maxp) ------------

// Number of embedding dimensions used to cluster rows before block-wise KNN.
static const int MC_EMBED_DIM = 20;

/**
 * @brief Low-dimensional row embedding through BigDataStatMeth's in-memory SVD.
 *        The matrix is zero-filled first (missing -> 0) purely so the
 *        decomposition runs; this fill only steers which rows cluster together,
 *        never the imputed values, which come from the observed data in the
 *        exact kernel. For a tall matrix (more CpGs than samples) BigDataStatMeth
 *        solves the small samples x samples Gram, so this stays cheap.
 * @param A CpGs x samples matrix (missing entries are NaN).
 * @param q Number of embedding dimensions requested.
 * @return A rows x q matrix of row coordinates (empty if the SVD failed).
 */
static Eigen::MatrixXd embed_rows_mem(const Eigen::MatrixXd& A, int q) {
    const int m = static_cast<int>(A.rows()), n = static_cast<int>(A.cols());
    Eigen::MatrixXd X = A;
    for (int c = 0; c < n; ++c)
        for (int r = 0; r < m; ++r)
            if (std::isnan(X(r, c))) X(r, c) = 0.0;
    int qq = std::min(q, std::min(m, n) - 1);
    if (qq < 1) qq = 1;
    BigDataStatMeth::svdeig s =
        BigDataStatMeth::RcppbdSVD(X, qq, 0, false, false);
    if (!s.bokuv || s.u.rows() != m || s.u.cols() < 1)
        return Eigen::MatrixXd();                  // caller falls back
    // Principal row coordinates: U scaled by the singular values.
    return s.u.array().rowwise() * s.d.transpose().array();
}

/**
 * @brief One 2-means pass (Lloyd) over the rows of an embedding. Deterministic
 *        init: the two rows at the extremes of the first coordinate, so the same
 *        data always yields the same split (no RNG).
 * @param E rows x q embedding to cluster.
 * @param iters Maximum Lloyd iterations.
 * @return A length-rows 0/1 cluster assignment.
 */
static std::vector<int> two_means(const Eigen::MatrixXd& E, int iters = 10) {
    const int m = static_cast<int>(E.rows());
    std::vector<int> cl(m, 0);
    if (m < 2) return cl;
    int lo = 0, hi = 0;
    for (int r = 1; r < m; ++r) {
        if (E(r, 0) < E(lo, 0)) lo = r;
        if (E(r, 0) > E(hi, 0)) hi = r;
    }
    if (lo == hi) return cl;
    Eigen::RowVectorXd c0 = E.row(lo), c1 = E.row(hi);
    for (int it = 0; it < iters; ++it) {
        bool changed = false;
        for (int r = 0; r < m; ++r) {
            const double d0 = (E.row(r) - c0).squaredNorm();
            const double d1 = (E.row(r) - c1).squaredNorm();
            const int a = (d0 <= d1) ? 0 : 1;
            if (a != cl[r]) { cl[r] = a; changed = true; }
        }
        if (it && !changed) break;
        Eigen::RowVectorXd s0 = Eigen::RowVectorXd::Zero(E.cols());
        Eigen::RowVectorXd s1 = Eigen::RowVectorXd::Zero(E.cols());
        int n0 = 0, n1 = 0;
        for (int r = 0; r < m; ++r)
            if (cl[r] == 0) { s0 += E.row(r); ++n0; }
            else            { s1 += E.row(r); ++n1; }
        if (n0 == 0 || n1 == 0) break;
        c0 = s0 / n0;
        c1 = s1 / n1;
    }
    return cl;
}

/**
 * @brief Recursively split rows (given by global indices into the embedding E)
 *        by 2-means until every block has at most maxp rows.
 * @param E rows x q embedding.
 * @param idx Global row indices of the current subset to split.
 * @param maxp Maximum block size; splitting stops once a block is this small.
 * @param out Output list; each finished block (a vector of global row indices)
 *        is appended.
 */
static void two_means_split(const Eigen::MatrixXd& E,
                            const std::vector<int>& idx, int maxp,
                            std::vector<std::vector<int>>& out) {
    if (static_cast<int>(idx.size()) <= maxp) { out.push_back(idx); return; }
    Eigen::MatrixXd sub(idx.size(), E.cols());
    for (size_t r = 0; r < idx.size(); ++r) sub.row(r) = E.row(idx[r]);
    std::vector<int> cl = two_means(sub);
    std::vector<int> a, b;
    for (size_t r = 0; r < idx.size(); ++r)
        (cl[r] == 0 ? a : b).push_back(idx[r]);
    if (a.empty() || b.empty()) {                  // could not separate: halve
        const size_t h = idx.size() / 2;
        a.assign(idx.begin(), idx.begin() + h);
        b.assign(idx.begin() + h, idx.end());
    }
    two_means_split(E, a, maxp, out);
    two_means_split(E, b, maxp, out);
}

/**
 * @brief KNN-impute a matrix in memory: exact when it fits in one block,
 *        otherwise split into blocks of similar rows (SVD embedding + 2-means)
 *        and impute each block exactly, so a large row set never forms a
 *        rows x rows distance matrix.
 * @param B CpGs x samples matrix; imputed in place (missing entries are NaN).
 * @param k Number of neighbours averaged for each imputed value.
 * @param maxp Largest block imputed exactly; larger row sets are split first.
 * @param rowmax Row missing-fraction threshold (see knn_impute_block).
 * @param colmax Column missing-fraction threshold (see knn_impute_block).
 * @return The number of columns exceeding colmax.
 */
static int knn_impute_inmem(Eigen::MatrixXd& B, int k, int maxp, double rowmax,
                            double colmax) {
    const int m = static_cast<int>(B.rows());
    if (m <= maxp) return knn_impute_block(B, k, rowmax, colmax);

    Eigen::MatrixXd E = embed_rows_mem(B, MC_EMBED_DIM);
    if (E.rows() != m)                              // SVD failed: one block
        return knn_impute_block(B, k, rowmax, colmax);

    std::vector<int> all(m);
    for (int r = 0; r < m; ++r) all[r] = r;
    std::vector<std::vector<int>> blocks;
    two_means_split(E, all, maxp, blocks);

    int overcol = 0;
    const int n = static_cast<int>(B.cols());
    for (const auto& idx : blocks) {
        Eigen::MatrixXd sub(idx.size(), n);
        for (size_t r = 0; r < idx.size(); ++r) sub.row(r) = B.row(idx[r]);
        overcol += knn_impute_block(sub, k, rowmax, colmax);
        for (size_t r = 0; r < idx.size(); ++r) B.row(idx[r]) = sub.row(r);
    }
    return overcol;
}

}  // namespace mcknn

//' KNN-impute a matrix in memory (C++/Eigen, block split via BigDataStatMeth)
//'
//' Fills the missing entries of a CpGs x samples matrix from each CpG's k
//' nearest co-varying CpGs. A row set within \code{maxp} is imputed exactly
//' (every CpG can be a neighbour of every other); a larger one is first split
//' into blocks of similar CpGs through a BigDataStatMeth SVD embedding and
//' 2-means, then each block is imputed exactly.
//'
//' @param Bin Numeric matrix (CpGs x samples); missing values are \code{NA}.
//' @param k Number of neighbours averaged for each imputed value.
//' @param maxp Largest block imputed exactly; larger row sets are split first.
//' @param rowmax A row with a greater missing fraction than this is filled with
//'   the per-column mean instead of by KNN.
//' @param colmax A column with a greater missing fraction than this is flagged
//'   (a warning is raised in R); imputation still proceeds.
//' @return The matrix with missing values filled.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix mc_knn_impute_mem(Rcpp::NumericMatrix Bin, int k = 10,
                                      int maxp = 1500, double rowmax = 0.5,
                                      double colmax = 0.8) {
    const int m = Bin.nrow(), n = Bin.ncol();
    Eigen::MatrixXd B(m, n);
    for (int c = 0; c < n; ++c)
        for (int r = 0; r < m; ++r)
            B(r, c) = Bin(r, c);
    const int overcol = mcknn::knn_impute_inmem(B, k, maxp, rowmax, colmax);
    if (overcol > 0)
        Rcpp::warning("KNN imputation: %d column(s) exceed colmax; imputed from "
                      "the few observed neighbours available.", overcol);
    Rcpp::NumericMatrix out(m, n);
    for (int c = 0; c < n; ++c)
        for (int r = 0; r < m; ++r)
            out(r, c) = B(r, c);
    out.attr("dimnames") = Bin.attr("dimnames");
    return out;
}

//' KNN-impute a matrix stored in HDF5 (reads through BigDataStatMeth)
//'
//' Reads a CpGs x samples matrix from an HDF5 dataset in R-view order through
//' the BigDataStatMeth C++ API (so the R/HDF5 layout is reconciled by the
//' library), imputes it with the same engine as \code{mc_knn_impute_mem}, and
//' returns the filled matrix. This validates the HDF5 read path against the
//' in-memory one; the memory-bounded block-wise variant builds on it.
//'
//' @param filename HDF5 file path.
//' @param group Group holding the dataset.
//' @param dataset Dataset name (CpGs x samples in R view).
//' @param k,maxp,rowmax,colmax See \code{mc_knn_impute_mem}.
//' @return The filled matrix (CpGs x samples), without dimnames.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix mc_knn_impute_hdf5(std::string filename, std::string group,
                                       std::string dataset, int k = 10,
                                       int maxp = 1500, double rowmax = 0.5,
                                       double colmax = 0.8) {
    Eigen::MatrixXd B = mc::read_full(filename, group, dataset);
    const int overcol = mcknn::knn_impute_inmem(B, k, maxp, rowmax, colmax);
    if (overcol > 0)
        Rcpp::warning("KNN imputation: %d column(s) exceed colmax.", overcol);
    return Rcpp::wrap(B);
}
