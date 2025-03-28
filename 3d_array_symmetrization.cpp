// File: 3d_array_symmetrization.cpp

#include <Rcpp.h>
using namespace Rcpp;

//' Symmetrize a 3D array of matrices
//'
//' Given a 3D array \code{X} with dimensions \(d \times d \times n\), this function
//' symmetrizes each \(d \times d\) matrix (slice) by computing \((X + t(X))/2\) for each slice.
//'
//' @param X A numeric array of dimensions \(d \times d \times n\)
//' @return The symmetrized 3D array.
//' @export
// [[Rcpp::export]]
NumericVector symmetrize_3d_array(NumericVector X) {
    // Retrieve dimensions; ensure X is a 3D array.
    IntegerVector dims = X.attr("dim");
    if (dims.size() != 3) {
        stop("Input must be a 3D array with dimensions (d, d, n).");
    }

    int d1 = dims[0];  // Number of rows
    int d2 = dims[1];  // Number of columns
    int n_slices = dims[2]; // Number of matrices

    if (d1 != d2) {
        stop("Each matrix slice must be square (number of rows must equal number of columns).");
    }

    int d = d1;  // Dimension of the square matrix

    // Loop over each matrix slice.
    for (int k = 0; k < n_slices; k++) {
        // Loop over rows and columns. We loop j from i to d to avoid redundant calculations.
        for (int i = 0; i < d; i++) {
            for (int j = i; j < d; j++) {
                // Calculate the index for element (i, j) and (j, i) in the 3D array.
                int idx1 = i + j * d + k * d * d; // Element (i, j) of kth slice
                int idx2 = j + i * d + k * d * d; // Element (j, i) of kth slice

                // Compute the average of the two symmetric entries.
                double avg = (X[idx1] + X[idx2]) / 2.0;
                X[idx1] = avg;
                X[idx2] = avg;
            }
        }
    }

    return X;
}
