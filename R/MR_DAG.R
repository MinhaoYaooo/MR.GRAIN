#' Estimate Causal DAG using Mendelian Randomization
#'
#' @description
#' `MR_DAG` infers a causal Directed Acyclic Graph (DAG) among a set of phenotypes (X) using
#' genetic instruments (Z). It first estimates the pairwise effects using two-stage least squares (2SLS)
#' methods and then projects the resulting matrix onto the space of DAGs to resolve cycles and
#' enforce consistency.
#'
#' @param X A numeric matrix of dimensions \eqn{n \times p} representing \eqn{n} samples and
#'   \eqn{p} phenotypes (exposures).
#' @param Z A numeric matrix of dimensions \eqn{n \times q} representing \eqn{n} samples and
#'   \eqn{q} genetic variants (instruments).
#' @param S A list of length \eqn{p}. \code{S[[j]]} contains the integer indices of the columns in \code{Z}
#'   that only have direct / cis-regulatory effects) for phenotype \code{X[, j]}.
#' @param lam Numeric. Sparsity penalty passed to \code{\link{notears_projection}}. Default is \code{0.0}.
#' @param max_iter Integer. Maximum iterations for projection. Default is \code{100}.
#' @param h_tol Numeric. Tolerance for acyclicity. Default is \code{1e-8}.
#' @param rho_max Numeric. Maximum penalty parameter. Default is \code{1e16}.
#' @param threshold_mode Character. Thresholding strategy passed to \code{\link{notears_projection}}.
#'   Either \code{"value"} (absolute cutoff) or \code{"quantile"} (data-driven cutoff based on the
#'   empirical distribution of edge magnitudes). Default is \code{"quantile"}.
#' @param w_threshold Numeric. Absolute threshold below which edge weights are set to zero in the
#'   final output when \code{threshold_mode = "value"}. Default is \code{0.1}.
#' @param q_threshold Numeric. Quantile threshold in \eqn{(0,1)} used when
#'   \code{threshold_mode = "quantile"}. Edges with magnitude below the
#'   \code{q_threshold}-quantile of the absolute off-diagonal weights are set to zero.
#'   Default is \code{0.95}.
#' @param quantile_ignore_zeros Logical. If \code{TRUE}, zero-valued weights are excluded when
#'   computing the quantile threshold. Default is \code{TRUE}.
#' @param ci Logical. If \code{TRUE}, computes standard errors, p-values, and confidence intervals
#'   for the non-zero entries of the estimated DAG. Default is \code{FALSE}.
#' @param alpha Numeric. Significance level if \code{ci=TRUE}. Default is \code{0.05}.
#'
#' @return
#' If \code{ci = FALSE}, returns a numeric matrix \eqn{B^*} of dimensions \eqn{p \times p}, representing
#' the estimated causal adjacency matrix.
#'
#' If \code{ci = TRUE}, returns a list containing:
#' \item{B_est}{The estimated weighted adjacency matrix \eqn{B^*}.}
#' \item{inference}{A data frame containing the Source, Target, Estimate, SE, P-value, and CI for each edge.}
#'
#' @details
#' The algorithm proceeds in two steps:
#' \enumerate{
#'   \item \strong{Initial Estimation:} It computes the reduced-form coefficients \eqn{\Gamma = X^T Z (Z^T Z)^{-1}}.
#'   Then, for each target phenotype \eqn{j}, it estimates the incoming edges from other phenotypes by
#'   regressing \eqn{\Gamma_{j, \text{valid}}} on \eqn{\Gamma_{-j, \text{valid}}}.
#'   \item \strong{Projection:} The matrix \eqn{B_{hat}} estimated in step 1 may contain cycles.
#'   The function calls \code{\link{notears_projection}} to find the nearest valid DAG.
#' }
#'
#' If \code{ci = TRUE}, post-selection inference is performed based on the asymptotic normality of the
#' MR-DAG estimator. This involves re-estimating residual variances from the reduced-form
#' equations and constructing a covariance matrix using the coefficients of non-target proteins.
#'
#' @importFrom MASS ginv
#' @importFrom stats lm.fit coef pnorm quantile
#' @export
MR_DAG <- function(X, Z, S,
                   lam = 0.0,
                   max_iter = 100,
                   h_tol = 1e-8,
                   rho_max = 1e16,
                   threshold_mode = c("quantile", "value"),
                   w_threshold = 0.1,
                   q_threshold = 0.95,
                   quantile_ignore_zeros = TRUE,
                   ci = FALSE,
                   alpha = 0.05) {

  # -------------------------------------------------------------
  # 1. Reduced-form Regression & Initialization
  # -------------------------------------------------------------
  # Center data for consistency with theory
  X_c <- scale(X, center = TRUE, scale = FALSE)
  Z_c <- scale(Z, center = TRUE, scale = FALSE)
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)

  # Gamma_hat = (X'Z)(Z'Z)^-1
  XTZ <- matrix(0, nrow = ncol(X_c), ncol = ncol(Z_c))
  for (j in 1:ncol(X_c)) {
    ok <- !is.na(X_c[, j])
    XTZ[j, ] <- crossprod(Z_c[ok, , drop = FALSE], X_c[ok, j])
  }
  ZTZ <- t(Z_c) %*% Z_c
  ZTZ_inv <- tryCatch(solve(ZTZ), error = function(e) MASS::ginv(ZTZ))
  Gamma_hat <- XTZ %*% ZTZ_inv

  # -------------------------------------------------------------
  # 2. Initial Matrix Estimate (B_hat)
  # -------------------------------------------------------------
  B_hat <- matrix(0, nrow = p, ncol = p)

  for (j in 1:p) {
    Sj <- S[[j]]
    all_snps <- 1:q
    notSj <- setdiff(all_snps, Sj)

    if (length(notSj) > 0) {
      y <- Gamma_hat[j, notSj]
      other_prots <- setdiff(1:p, j)
      Xj <- Gamma_hat[other_prots, notSj, drop = FALSE]
      Xj_t <- t(Xj)

      ok <- !is.na(y) & stats::complete.cases(Xj_t)

      if (sum(ok) > 0) {
        fit <- stats::lm.fit(x = Xj_t[ok, , drop = FALSE], y = y[ok])
        coefs <- stats::coef(fit)
        coefs[is.na(coefs)] <- 0
      } else {
        coefs <- rep(0, length(other_prots))
      }

      B_hat[j, other_prots] <- coefs
    }
  }

  # -------------------------------------------------------------
  # 3. NOTEARS Projection
  # -------------------------------------------------------------
  B_star <- notears_projection(B_init = B_hat,
                               lam = lam,
                               max_iter = max_iter,
                               h_tol = h_tol,
                               rho_max = rho_max,
                               threshold_mode = threshold_mode,
                               w_threshold = w_threshold,
                               q_threshold = q_threshold,
                               quantile_ignore_zeros = quantile_ignore_zeros)

  if (!is.null(colnames(X))) {
    rownames(B_star) <- colnames(X)
    colnames(B_star) <- colnames(X)
  } else {
    rownames(B_star) <- 1:p
    colnames(B_star) <- 1:p
  }

  if (!ci) {
    return(B_star)
  }

  # -------------------------------------------------------------
  # 4. Statistical Inference (Theorem 3)
  # -------------------------------------------------------------
  results_df <- data.frame(
    target = integer(),
    source = integer(),
    est = numeric(),
    se = numeric(),
    pval = numeric(),
    ci_lower = numeric(),
    ci_upper = numeric()
  )

  for (j in 1:p) {
    parents <- which(abs(B_star[j, ]) > 1e-5)
    if (length(parents) == 0) next

    Gamma_j_hat <- Gamma_hat[j, ]
    ok <- !is.na(X_c[, j])

    resid_rf <- X_c[ok, j] - (Z_c[ok, , drop = FALSE] %*% Gamma_j_hat)
    sigma_sq <- sum(resid_rf^2) / sum(ok)

    notSj <- setdiff(1:q, S[[j]])
    if (length(notSj) == 0) next

    other_nodes <- setdiff(1:p, j)
    Gamma_full <- Gamma_hat[other_nodes, notSj, drop = FALSE]
    Sigma_Gamma_full <- Gamma_full %*% t(Gamma_full)

    Sigma_inv_full <- tryCatch(solve(Sigma_Gamma_full),
                               error = function(e) MASS::ginv(Sigma_Gamma_full))

    parent_idx <- match(parents, other_nodes)
    Sigma_inv <- Sigma_inv_full[parent_idx, parent_idx, drop = FALSE]

    var_matrix <- (sigma_sq / n) * Sigma_inv
    var_diag <- diag(var_matrix)
    var_diag[var_diag < 0] <- 0
    se_vec <- sqrt(var_diag)

    for (k in seq_along(parents)) {
      parent_node <- parents[k]
      est_val <- B_star[j, parent_node]
      se_val <- se_vec[k]

      if (se_val < 1e-12) {
        t_stat <- Inf * sign(est_val)
        p_val <- 0.0
      } else {
        t_stat <- est_val / se_val
        p_val <- 2 * (1 - stats::pnorm(abs(t_stat)))
      }

      ci_lower <- est_val - stats::qnorm(1 - alpha/2) * se_val
      ci_upper <- est_val + stats::qnorm(1 - alpha/2) * se_val

      results_df <- rbind(results_df, data.frame(
        target = rownames(B_star)[j],
        source = rownames(B_star)[parent_node],
        est = est_val,
        se = se_val,
        pval = p_val,
        ci_lower = ci_lower,
        ci_upper = ci_upper
      ))
    }
  }

  list(B_est = B_star, inference = results_df)
}
