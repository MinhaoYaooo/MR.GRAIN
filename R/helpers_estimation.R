.estimate_gamma_hat <- function(X, Z, inv_ridge) {
  p <- ncol(X)
  q <- ncol(Z)

  Gamma_hat <- matrix(NA_real_, nrow = p, ncol = q)
  rownames(Gamma_hat) <- colnames(X)
  colnames(Gamma_hat) <- colnames(Z)

  XTZ_list <- vector("list", p)
  ZTZ_list <- vector("list", p)
  ZTZ_inv_list <- vector("list", p)

  for (j in 1:p) {
    ok <- !is.na(X[, j])
    n_j <- sum(ok)

    if (n_j == 0) {
      stop(sprintf(
        "Column %d of X is entirely NA after excluding empty S_j biomolecules.",
        j
      ))
    }

    Z_j <- Z[ok, , drop = FALSE]
    x_j <- X[ok, j]

    XTZ_j <- as.numeric(crossprod(Z_j, x_j))
    ZTZ_j <- crossprod(Z_j)
    ZTZ_inv_j <- .safe_inverse(ZTZ_j, ridge = inv_ridge)

    Gamma_hat[j, ] <- as.numeric(ZTZ_inv_j %*% XTZ_j)

    XTZ_list[[j]] <- XTZ_j
    ZTZ_list[[j]] <- ZTZ_j
    ZTZ_inv_list[[j]] <- ZTZ_inv_j
  }

  list(
    Gamma_hat = Gamma_hat,
    XTZ_list = XTZ_list,
    ZTZ_list = ZTZ_list,
    ZTZ_inv_list = ZTZ_inv_list
  )
}


.estimate_lambda_hat <- function(Gamma_hat, S_idx, center_tol) {
  p <- nrow(Gamma_hat)

  Lambda_hat <- matrix(0, nrow = p, ncol = p)
  G_hat <- matrix(0, nrow = p, ncol = p)
  rownames(Lambda_hat) <- rownames(Gamma_hat)
  colnames(Lambda_hat) <- rownames(Gamma_hat)
  rownames(G_hat) <- rownames(Gamma_hat)
  colnames(G_hat) <- rownames(Gamma_hat)

  for (j in 1:p) {
    Sj_idx <- S_idx[[j]]

    g_j <- rowMeans(Gamma_hat[, Sj_idx, drop = FALSE])
    G_hat[, j] <- g_j

    if (!is.finite(g_j[j]) || abs(g_j[j]) < center_tol) {
      stop(sprintf(
        "Average reduced-form effect for node %d has near-zero j-th entry; cannot normalize.",
        j
      ))
    }

    Lambda_hat[, j] <- g_j / g_j[j]
  }

  diag(Lambda_hat) <- 1.0

  list(
    Lambda_hat = Lambda_hat,
    G_hat = G_hat
  )
}
