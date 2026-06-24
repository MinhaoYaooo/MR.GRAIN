.notears_projection <- function(B_init,
                                lam = 0.0,
                                max_iter = 100,
                                h_tol = 1e-8,
                                rho_max = 1e16,
                                w_threshold = 0.1) {

  B_init <- as.matrix(B_init)

  if (nrow(B_init) != ncol(B_init)) {
    stop("B_init must be a square matrix.")
  }

  if (any(!is.finite(B_init))) {
    stop("B_init must contain only finite values.")
  }

  p <- nrow(B_init)

  get_adj <- function(w) {
    w_pos <- w[1:(p * p)]
    w_neg <- w[(p * p + 1):(2 * p * p)]
    matrix(w_pos - w_neg, nrow = p, ncol = p)
  }

  calc_loss <- function(W) {
    D <- W - B_init
    loss <- 0.5 * sum(D^2)
    G_loss <- D
    list(loss = loss, G_loss = G_loss)
  }

  calc_h <- function(W) {
    E <- expm::expm(W * W)
    h_val <- sum(diag(E)) - p
    G_h <- t(E) * W * 2
    list(h = h_val, G_h = G_h)
  }

  rho <- 1.0
  alpha <- 0.0

  fn_obj <- function(w) {
    W <- get_adj(w)
    l_res <- calc_loss(W)
    h_res <- calc_h(W)
    l_res$loss + 0.5 * rho * h_res$h^2 + alpha * h_res$h + lam * sum(w)
  }

  gr_obj <- function(w) {
    W <- get_adj(w)
    l_res <- calc_loss(W)
    h_res <- calc_h(W)
    G_smooth <- l_res$G_loss + (rho * h_res$h + alpha) * h_res$G_h
    c(as.vector(G_smooth) + lam, as.vector(-G_smooth) + lam)
  }

  w_est <- numeric(2 * p * p)
  h_val <- Inf

  lower_b <- numeric(2 * p * p)
  upper_b <- rep(Inf, 2 * p * p)

  diag_idx <- (0:(p - 1)) * p + (1:p)
  lower_b[diag_idx] <- 0
  upper_b[diag_idx] <- 0
  lower_b[diag_idx + p * p] <- 0
  upper_b[diag_idx + p * p] <- 0

  for (iter in 1:max_iter) {
    while (rho < rho_max) {
      res <- stats::optim(
        par = w_est,
        fn = fn_obj,
        gr = gr_obj,
        method = "L-BFGS-B",
        lower = lower_b,
        upper = upper_b
      )

      w_new <- res$par
      h_new <- calc_h(get_adj(w_new))$h

      if (h_new > 0.25 * h_val) {
        rho <- rho * 10
      } else {
        break
      }
    }

    w_est <- w_new
    h_val <- h_new
    alpha <- alpha + rho * h_val

    if (h_val <= h_tol || rho >= rho_max) {
      break
    }
  }

  B_dag <- get_adj(w_est)

  B_dag[abs(B_dag) < w_threshold] <- 0.0
  diag(B_dag) <- 0.0

  return(B_dag)
}
