#' Estimate Causal DAG using Mendelian Randomization
#'
#' @description
#' \code{MR_DAG} infers a causal Directed Acyclic Graph (DAG) among a set of phenotypes (X) using
#' genetic instruments (Z). It first estimates the pairwise effects using Instrumental Variable (IV)
#' methods and then projects the resulting matrix onto the space of DAGs to resolve cycles and
#' enforce consistency.
#'
#' @param X A numeric matrix of dimensions \eqn{n \times p} representing \eqn{n} samples and
#'   \eqn{p} phenotypes (exposures).
#' @param Z A numeric matrix of dimensions \eqn{n \times q} representing \eqn{n} samples and
#'   \eqn{q} genetic variants (instruments).
#' @param S A list of length \eqn{p}. \code{S[[j]]} contains the integer indices of the columns in \code{Z}
#'   that are valid instruments for phenotype \code{X[, j]}.
#' @param lam Numeric. Sparsity penalty passed to \code{\link{notears_projection}}. Default is \code{0.0}.
#' @param max_iter Integer. Maximum iterations for projection. Default is \code{100}.
#' @param h_tol Numeric. Tolerance for acyclicity. Default is \code{1e-8}.
#' @param rho_max Numeric. Maximum penalty parameter. Default is \code{1e16}.
#' @param w_threshold Numeric. Threshold for zeroing out small coefficients. Default is \code{0.1}.
#'
#' @return A numeric matrix \eqn{B^*} of dimensions \eqn{p \times p}, representing the estimated causal
#'   adjacency matrix. Entry \eqn{B^*_{ij}} represents the causal effect of \eqn{X_j} on \eqn{X_i}.
#'
#' @details
#' The algorithm proceeds in two steps:
#' \enumerate{
#'   \item \strong{Initial Estimation:} It computes the reduced-form coefficients \eqn{\Gamma = X^T Z (Z^T Z)^{-1}}.
#'   Then, for each target phenotype \eqn{j}, it estimates the incoming edges from other phenotypes by
#'   regressing \eqn{\Gamma_{j, -\text{cis-SNPs}}} on \eqn{\Gamma_{-j, -\text{cis-SNPs}}}.
#'   \item \strong{Projection:} The matrix \eqn{B_{hat}} estimated in step 1 may contain cycles.
#'   The function calls \code{\link{notears_projection}} to find the nearest valid DAG.
#' }
#'
#' @importFrom MASS ginv
#' @importFrom stats lm.fit coef
#' @export
#'
MR_DAG <- function(X, Z, S,
                   lam = 0.0,
                   max_iter = 100,
                   h_tol = 1e-8,
                   rho_max = 1e16,
                   w_threshold = 0.1) {

  # 1. Compute Cross-products
  XTZ <- t(X) %*% Z
  ZTZ <- t(Z) %*% Z

  # 2. Invert ZTZ (Robustly)
  ZTZ_inv <- tryCatch(solve(ZTZ), error = function(e) MASS::ginv(ZTZ))

  # 3. Reduced form coefficients Gamma (p x q)
  Gamma_hat <- XTZ %*% ZTZ_inv

  p <- ncol(X)
  q <- ncol(Z)
  B_hat <- matrix(0, nrow = p, ncol = p)

  # 4. Estimate B_hat row by row
  for (j in 1:p) {
    Sj <- S[[j]]
    all_snps <- 1:q
    notSj <- setdiff(all_snps, Sj)

    # We regress Gamma of target j on Gamma of other phenotypes
    # using only instruments NOT associated with j (validity condition)
    if (length(notSj) > 0) {
      y <- Gamma_hat[j, notSj]
      other_prots <- setdiff(1:p, j)
      Xj <- Gamma_hat[other_prots, notSj, drop = FALSE]
      Xj_t <- t(Xj)

      # Fast linear regression
      fit <- stats::lm.fit(x = as.matrix(Xj_t), y = y)
      coefs <- stats::coef(fit)
      coefs[is.na(coefs)] <- 0

      B_hat[j, other_prots] <- coefs
    }
  }

  # 5. Project B_hat onto DAG space
  B_star <- notears_projection(B_init = B_hat,
                               lam = lam,
                               max_iter = max_iter,
                               h_tol = h_tol,
                               rho_max = rho_max,
                               w_threshold = w_threshold)
  return(B_star)
}
