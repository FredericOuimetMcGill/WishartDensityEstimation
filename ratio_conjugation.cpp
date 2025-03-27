// ratio_conjugation.cpp
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;

//' Conjugated ratio of two SPD matrix arrays
//'
//' Given two 3D arrays `S1` and `S2` (each of dimension p x p x n), compute
//' for each index i in {1,...,n} the matrix
//' \eqn{ S2_i^{-1/2} * S1_i * S2_i^{-1/2} }.
//'
//' @param S1 A 3D array (p x p x n) in \code{NumericVector} form with dim attribute.
//' @param S2 A 3D array (p x p x n) in \code{NumericVector} form with dim attribute.
//' @return A 3D array (p x p x n) whose i-th slice is
//' \eqn{ S2_i^{-1/2} * S1_i * S2_i^{-1/2} }.
//' @export
// [[Rcpp::export]]
NumericVector conjugated_ratio(const NumericVector& S1, const NumericVector& S2)
{
    // Check dimensions
    if (!S1.hasAttribute("dim") || !S2.hasAttribute("dim")) {
        stop("Both S1 and S2 must have 'dim' attributes (be 3D arrays).");
    }

    IntegerVector dimS1 = S1.attr("dim");
    IntegerVector dimS2 = S2.attr("dim");

    if (dimS1.size() != 3 || dimS2.size() != 3) {
        stop("S1 and S2 must be 3D arrays.");
    }

    int p1 = dimS1[0];
    int p2 = dimS1[1];
    int n1 = dimS1[2];

    int p1b = dimS2[0];
    int p2b = dimS2[1];
    int n2 = dimS2[2];

    if (p1 != p2 || p1b != p2b || p1 != p1b) {
        stop("Non-conforming dimensions: S1 and S2 must be p x p x n with the same p.");
    }
    if (n1 != n2) {
        stop("S1 and S2 must have the same 'n' in the third dimension.");
    }

    int p = p1;  // dimension of each matrix
    int n = n1;  // number of slices

    // Prepare output: p x p x n
    NumericVector out(p * p * n);
    out.attr("dim") = Dimension(p, p, n);

    // Loop over each slice
    // We'll map each slice of S1 and S2 into Eigen matrices, compute S2^-1/2,
    // and then do the product: S2^-1/2 * S1 * S2^-1/2.

    for (int i = 0; i < n; i++) {
        // Map the i-th slice of S1
        Eigen::Map<const Eigen::MatrixXd> M1(&S1[i * p * p], p, p);
        // Map the i-th slice of S2
        Eigen::Map<const Eigen::MatrixXd> M2(&S2[i * p * p], p, p);

        // Compute the symmetric inverse-half of M2
        // 1) Eigen-decomposition: M2 = Q * D * Q^T (where D is diagonal of eigenvalues)
        // 2) Then M2^-1/2 = Q * D^(-1/2) * Q^T
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M2);
        if (es.info() != Eigen::Success) {
            stop("Eigen decomposition failed for slice %d of S2.", i + 1);
        }

        // Eigenvalues and eigenvectors
        Eigen::VectorXd evals = es.eigenvalues();
        Eigen::MatrixXd evecs = es.eigenvectors();

        // Take negative half power of the eigenvalues
        // (i.e., for each lambda, we want lambda^(-1/2))
        for (int k = 0; k < p; k++) {
            if (evals[k] <= 0.0) {
                stop("Matrix S2 slice %d is not positive definite (non-positive eigenvalue).", i + 1);
            }
        }
        Eigen::VectorXd evals_mhalf = evals.array().sqrt().inverse();  // 1 / sqrt(lambda)

        // Build diagonal matrix of lambda^(-1/2)
        Eigen::MatrixXd D_mhalf = evals_mhalf.asDiagonal();

        // Now compute S2^-1/2 = evecs * D_mhalf * evecs^T
        Eigen::MatrixXd S2_invhalf = evecs * D_mhalf * evecs.transpose();

        // Compute the conjugation: S2_invhalf * M1 * S2_invhalf
        Eigen::MatrixXd C = S2_invhalf * M1 * S2_invhalf;

        // Copy the result into the output
        Eigen::Map<Eigen::MatrixXd> outSlice(&out[i * p * p], p, p);
        outSlice = C;
    }

    return out;
}
