// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;

// This function expects a NumericVector with a "dim" attribute (a 3D array)
// with dimensions p x p x n, and returns a 3D array of the same size
// where each p x p slice is the inverse of the corresponding input matrix.
// [[Rcpp::export]]
NumericVector invert_3d_array(NumericVector arr) {
    // Get the dimensions of the array
    IntegerVector dims = arr.attr("dim");
    if (dims.size() != 3)
        stop("Input must be a 3D array");

    int p = dims[0];  // number of rows (and columns)
    int n = dims[2];  // number of matrices

    // Create an output vector of the same size and set its dimensions
    NumericVector out(arr.size());
    out.attr("dim") = dims;

    // Loop over each matrix slice (the matrices are stored contiguously)
    for (int k = 0; k < n; k++) {
        // Create an Eigen map for the k-th matrix in the input
        Eigen::Map<Eigen::MatrixXd> mat(arr.begin() + k * p * p, p, p);
        // Create an Eigen map for the k-th matrix in the output
        Eigen::Map<Eigen::MatrixXd> inv_mat(out.begin() + k * p * p, p, p);

        // Compute the inverse using Eigen's highly optimized routines
        inv_mat = mat.inverse();
    }

    return out;
}