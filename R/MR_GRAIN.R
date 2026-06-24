#' Fit MR-GRAIN
#'
#' Reconstruct a directed acyclic biomolecular network using genetically anchored instrumental variables.
#'
#' @param X Matrix of biomolecular traits. Rows are samples and columns are biomolecules.
#' @param Z Matrix of genetic instruments. Rows are samples and columns are SNPs.
#' @param S List of SNP sets. `S[[j]]` contains the SNP indices or SNP names for biomolecule `j`.
#' @param lam L1 penalty parameter used in NOTEARS projection.
#' @param max_iter Maximum number of NOTEARS augmented-Lagrangian iterations.
#' @param h_tol Tolerance for the acyclicity constraint.
#' @param rho_max Maximum augmented-Lagrangian penalty.
#' @param w_threshold Threshold for pruning small edge weights after DAG projection.
#' @param center_tol Tolerance for detecting near-zero normalization constants.
#' @param inv_ridge Ridge parameter used in safe matrix inversion.
#' @param inference Logical. If `TRUE`, compute standard inference for active DAG edges. Default is `FALSE`.
#'
#' @return A list containing `B_init`, `B_DAG`, `sd`, and `pval`.
#' @export
MR_GRAIN <- function(X, Z, S,
                     lam = 0.0,
                     max_iter = 100,
                     h_tol = 1e-8,
                     rho_max = 1e16,
                     w_threshold = 0.1,
                     center_tol = 1e-8,
                     inv_ridge = 1e-6,
                     inference = FALSE) {

  ########## Input checks for inference ##########
  if (!is.logical(inference) || length(inference) != 1 || is.na(inference)) {
    stop("inference must be TRUE or FALSE.")
  }

  X <- as.matrix(X)
  Z <- as.matrix(Z)

  n <- nrow(X)
  p0 <- ncol(X)

  if (nrow(Z) != n) {
    stop("X and Z must have the same number of rows.")
  }
  if (length(S) != p0) {
    stop("S must be a list of length ncol(X).")
  }

  ########## Resolve S once ##########
  S_idx_full <- .resolve_S_indices(S, Z)

  ########## Exclude biomolecules with |S_j| = 0 ##########
  keep_idx <- which(lengths(S_idx_full) > 0)

  if (length(keep_idx) == 0) {
    stop("All biomolecules have empty S_j; nothing to estimate.")
  }

  X <- X[, keep_idx, drop = FALSE]
  S_idx <- S_idx_full[keep_idx]
  keep_names <- colnames(X)

  p <- ncol(X)

  ########## 1. Marginal reduced-form effect estimate ##########
  gamma_fit <- .estimate_gamma_hat(
    X = X,
    Z = Z,
    inv_ridge = inv_ridge
  )

  Gamma_hat <- gamma_fit$Gamma_hat
  ZTZ_list <- gamma_fit$ZTZ_list

  ########## 2. Estimate Lambda_hat ##########
  lambda_fit <- .estimate_lambda_hat(
    Gamma_hat = Gamma_hat,
    S_idx = S_idx,
    center_tol = center_tol
  )

  Lambda_hat <- lambda_fit$Lambda_hat
  G_hat <- lambda_fit$G_hat

  ########## 3. Initial unconstrained estimate ##########
  Lambda_inv <- .safe_inverse(Lambda_hat, ridge = inv_ridge)
  B_init <- diag(p) - Lambda_inv
  diag(B_init) <- 0.0
  rownames(B_init) <- keep_names
  colnames(B_init) <- keep_names

  ########## 4. DAG projection by NOTEARS ##########
  B_DAG <- .notears_projection(
    B_init = B_init,
    lam = lam,
    max_iter = max_iter,
    h_tol = h_tol,
    rho_max = rho_max,
    w_threshold = w_threshold
  )

  rownames(B_DAG) <- keep_names
  colnames(B_DAG) <- keep_names

  ########## 5. Standard inference if requested ##########
  sd_mat <- matrix(
    NA_real_,
    nrow = p,
    ncol = p,
    dimnames = list(keep_names, keep_names)
  )

  p_mat <- matrix(
    NA_real_,
    nrow = p,
    ncol = p,
    dimnames = list(keep_names, keep_names)
  )

  if (isTRUE(inference)) {
    cov_fit <- .compute_lambda_cov(
      X = X,
      Z = Z,
      Gamma_hat = Gamma_hat,
      S_idx = S_idx,
      ZTZ_list = ZTZ_list,
      Lambda_hat = Lambda_hat,
      G_hat = G_hat,
      n = n,
      inv_ridge = inv_ridge
    )

    V_Lambda_hat <- cov_fit$V_Lambda_hat

    H_B_hat <- kronecker(t(Lambda_inv), Lambda_inv)

    V_B_hat <- H_B_hat %*% V_Lambda_hat %*% t(H_B_hat)

    se_all <- matrix(
      sqrt(pmax(diag(V_B_hat) / n, 0)),
      nrow = p,
      ncol = p,
      dimnames = list(keep_names, keep_names)
    )

    active_idx <- which(B_DAG != 0)

    if (length(active_idx) > 0) {
      sd_mat[active_idx] <- se_all[active_idx]

      valid_idx <- active_idx[
        is.finite(se_all[active_idx]) &
          (se_all[active_idx] > 0)
      ]

      if (length(valid_idx) > 0) {
        z_stat <- B_DAG[valid_idx] / se_all[valid_idx]
        p_mat[valid_idx] <- 2 * stats::pnorm(-abs(z_stat))
      }
    }
  }

  ########## 6. Return ##########
  out <- list(
    B_init = B_init,
    B_DAG = B_DAG,
    sd = sd_mat,
    pval = p_mat
  )

  return(out)
}
