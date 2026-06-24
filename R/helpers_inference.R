.compute_lambda_cov <- function(X, Z, Gamma_hat, S_idx, ZTZ_list,
                                Lambda_hat, G_hat, n, inv_ridge) {
  p <- ncol(X)
  q <- ncol(Z)

  A_mat <- matrix(0, nrow = q, ncol = p)
  rownames(A_mat) <- colnames(Z)
  colnames(A_mat) <- colnames(X)

  for (j in 1:p) {
    Sj_idx <- S_idx[[j]]
    A_mat[Sj_idx, j] <- 1 / length(Sj_idx)
  }

  Eps_hat <- matrix(0, nrow = n, ncol = p)
  rownames(Eps_hat) <- rownames(X)
  colnames(Eps_hat) <- colnames(X)

  obs_mat <- !is.na(X)

  for (j in 1:p) {
    ok <- obs_mat[, j]
    if (any(ok)) {
      Eps_hat[ok, j] <- X[ok, j] - Z[ok, , drop = FALSE] %*% Gamma_hat[j, ]
    }
  }

  psi_G <- matrix(0, nrow = n, ncol = p * p)

  for (r in 1:p) {
    M_hat_r <- ZTZ_list[[r]] / n
    M_hat_r_inv <- .safe_inverse(M_hat_r, ridge = inv_ridge)

    score_r <- Z %*% M_hat_r_inv %*% A_mat

    eps_r <- Eps_hat[, r]
    obs_r <- obs_mat[, r]
    contrib_r <- as.numeric(obs_r) * eps_r

    for (c in 1:p) {
      idx_rc <- r + (c - 1) * p
      psi_G[, idx_rc] <- contrib_r * score_r[, c]
    }
  }

  V_G_hat <- crossprod(psi_G) / n

  D_list <- vector("list", p)

  for (j in 1:p) {
    e_j <- rep(0, p)
    e_j[j] <- 1
    D_list[[j]] <- (diag(p) - tcrossprod(Lambda_hat[, j], e_j)) / G_hat[j, j]
  }

  J_Lambda_hat <- as.matrix(Matrix::bdiag(D_list))
  V_Lambda_hat <- J_Lambda_hat %*% V_G_hat %*% t(J_Lambda_hat)

  list(V_Lambda_hat = V_Lambda_hat)
}
