/**
 * @file methylclock_io.h
 * @brief Read HDF5 datasets into Eigen matrices in R-view order through the
 *        BigDataStatMeth C++ API, so the row/column-major reconciliation between
 *        R, HDF5 and C++ is handled by the library, never by hand. Include after
 *        BigDataStatMeth.hpp.
 */
#ifndef METHYLCLOCK_IO_HPP
#define METHYLCLOCK_IO_HPP

#include <BigDataStatMeth.hpp>

namespace mc {

/**
 * @brief Read a full dataset into an Eigen matrix in R-view order.
 * @param fn HDF5 file path.
 * @param grp Group holding the dataset.
 * @param ds Dataset name.
 * @return The dataset as an Eigen matrix (rows = nrows_r, cols = ncols_r).
 */
inline Eigen::MatrixXd read_full(const std::string& fn, const std::string& grp,
                                 const std::string& ds) {
    std::unique_ptr<BigDataStatMeth::hdf5Dataset> d(
        new BigDataStatMeth::hdf5Dataset(fn, grp, ds, false));
    d->openDataset();
    if (d->getDatasetptr() == nullptr)
        throw std::runtime_error("cannot open " + grp + "/" + ds);
    const std::size_t nr = d->nrows_r(), nc = d->ncols_r();
    Eigen::MatrixXd M(nr, nc);
    std::vector<hsize_t> off = {0, 0}, cnt = {(hsize_t)nc, (hsize_t)nr},
                         st = {1, 1}, bl = {1, 1};
    d->readDatasetBlock(off, cnt, st, bl, M.data());
    return M;
}

/**
 * @brief Read a dataset stored as 1 x N (or N x 1) into an Eigen vector.
 * @param fn HDF5 file path.
 * @param grp Group holding the dataset.
 * @param ds Dataset name.
 * @return The dataset flattened into an Eigen vector.
 */
inline Eigen::VectorXd read_vec(const std::string& fn, const std::string& grp,
                                const std::string& ds) {
    Eigen::MatrixXd M = read_full(fn, grp, ds);
    return Eigen::Map<Eigen::VectorXd>(M.data(), M.size());
}

/**
 * @brief Read a contiguous range of R-view rows [r0, r0 + rc) as an rc x ncols
 *        matrix from an already-open dataset. The row index is the HDF5 fast
 *        dimension, so a contiguous CpG range is a single hyperslab; the (nc, rc)
 *        count mirrors read_full.
 * @param d Open dataset to read from.
 * @param r0 First R-view row of the range.
 * @param rc Number of rows to read.
 * @return The rc x ncols block as an Eigen matrix.
 */
inline Eigen::MatrixXd read_rows(BigDataStatMeth::hdf5Dataset* d,
                                 std::size_t r0, std::size_t rc) {
    const std::size_t nc = d->ncols_r();
    Eigen::MatrixXd M(rc, nc);
    std::vector<hsize_t> off = {0, (hsize_t)r0}, cnt = {(hsize_t)nc, (hsize_t)rc},
                         st = {1, 1}, bl = {1, 1};
    d->readDatasetBlock(off, cnt, st, bl, M.data());
    return M;
}

/**
 * @brief Write an rc x ncols R-view block back to rows [r0, r0 + rc) of an open
 *        dataset, symmetric to read_rows (same offset and count, same
 *        column-major buffer), so a write then read round-trips. The dataset must
 *        already exist with the layout read_full expects (e.g. created from R
 *        with hdf5_create_matrix).
 * @param d Open dataset to write to.
 * @param r0 First R-view row of the destination range.
 * @param M The rc x ncols block to write (rc = M.rows()).
 */
inline void write_rows(BigDataStatMeth::hdf5Dataset* d, std::size_t r0,
                       const Eigen::MatrixXd& M) {
    const std::size_t rc = M.rows(), nc = M.cols();
    std::vector<double> buf(M.data(), M.data() + rc * nc);
    std::vector<hsize_t> off = {0, (hsize_t)r0}, cnt = {(hsize_t)nc, (hsize_t)rc},
                         st = {1, 1}, bl = {1, 1};
    d->writeDatasetBlock(buf, off, cnt, st, bl);
}

}  // namespace mc
#endif
