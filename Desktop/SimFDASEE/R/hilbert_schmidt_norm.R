#' Compute the Hilbert–Schmidt Norm of a Matrix
#'
#' This function computes the Hilbert–Schmidt norm (also known as the Frobenius norm)
#' of a numeric matrix \(A\), defined as the square root of the sum of the squared entries.
#'
#' @param A A numeric matrix.
#'
#' @return A single numeric value representing the Hilbert–Schmidt norm of the matrix.
#'
#' @examples
#' mat <- matrix(1:4, nrow = 2)
#' hilbert_schmidt_norm(mat)
#'
#' @export
hilbert_schmidt_norm <- function(A) {
  sqrt(sum(A^2))
}
