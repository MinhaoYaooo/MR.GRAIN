#' Project a Matrix onto the DAG Space (NOTEARS)
#'
#' This function solves a continuous optimization problem to find the weighted adjacency matrix
#' \eqn{W} representing a Directed Acyclic Graph (DAG) that is closest (in the Frobenius norm sense)
#' to an input initialization matrix \eqn{B_{init}}. It uses the augmented Lagrangian method to
#' enforce the acyclicity constraint.
#'
#' @param B_init A numeric matrix of shape \eqn{p \times p}. This is the initial estimate of the
#'   causal graph (e.g., obtained via IV regression) that likely contains cycles.
#' @param lam Numeric. The \eqn{\ell_1} penalty parameter for sparsity. Default is \code{0.0}.
#'   Higher values encourage a sparser graph.
#' @param max_iter Integer. The maximum number of iterations for the augmented Lagrangian loop.
#'   Default is \code{100}.
#' @param h_tol Numeric. The tolerance for the acyclicity constraint \eqn{h(W) = 0}.
#'   Default is \code{1e-8}.
#' @param rho_max Numeric. The maximum value for the penalty parameter \eqn{\rho}.
#'   Default is \code{1e16}.
#' @param threshold_mode Character. Thresholding strategy for sparsifying the final DAG weights.
#'   Either \code{"value"} (absolute cutoff) or \code{"quantile"} (data-driven cutoff based on
#'   the empirical distribution of edge magnitudes). Default is \code{"quantile"}.
#' @param w_threshold Numeric. Absolute threshold below which edge weights are set to zero in the
#'   final output when \code{threshold_mode = "value"}. Default is \code{0.1}.
#' @param q_threshold Numeric. Quantile threshold in \eqn{(0,1)} used when
#'   \code{threshold_mode = "quantile"}. Edges with magnitude below the
#'   \code{q_threshold}-quantile of the absolute off-diagonal weights are set to zero.
#'   Default is \code{0.95}.
#' @param quantile_ignore_zeros Logical. If \code{TRUE}, zero-valued weights are excluded when
#'   computing the quantile threshold. Default is \code{TRUE}.
#'
#' @return A numeric matrix \eqn{B^{\dagger}} (p x p) representing the projected DAG.
#'   The matrix is weighted and strictly acyclic (within numerical tolerance).
#'
#' @details
#' The optimization problem solves:
#' \deqn{\min_{W} \frac{1}{2} ||W - B_{init}||_F^2 + \lambda ||W||_1}
#' subject to the acyclicity constraint:
#' \deqn{h(W) = \text{tr}(e^{W \circ W}) - p = 0}
#'
#' After optimization, small edge weights are removed using either an absolute threshold
#' (\code{threshold_mode = "value"}) or a data-driven quantile threshold
#' (\code{threshold_mode = "quantile"}), applied to the absolute off-diagonal entries of
#' the estimated adjacency matrix.
#'
#' @importFrom expm expm
#' @importFrom stats optim quantile
#' @export
notears_projection <- function(B_init,
                               lam = 0.0,
                               max_iter = 100,
                               h_tol = 1e-8,
                               rho_max = 1e16,
                               threshold_mode = c("quantile", "value"),
                               w_threshold = 0.1,
                               q_threshold = 0.95,
                               quantile_ignore_zeros = TRUE) {

  p <- nrow(B_init)

  get_adj <- function(w) {
    w_pos <- w[1:(p*p)]
    w_neg <- w[(p*p + 1):(2*p*p)]
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

  diag_idx <- (0:(p-1)) * p + (1:p)
  lower_b[diag_idx] <- 0; upper_b[diag_idx] <- 0
  lower_b[diag_idx + p*p] <- 0; upper_b[diag_idx + p*p] <- 0

  for (iter in 1:max_iter) {
    while (rho < rho_max) {
      res <- stats::optim(par = w_est, fn = fn_obj, gr = gr_obj,
                          method = "L-BFGS-B", lower = lower_b, upper = upper_b)
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
    if (h_val <= h_tol || rho >= rho_max) break
  }

  B_dag <- get_adj(w_est)

  threshold_mode <- match.arg(threshold_mode)

  if (threshold_mode == "value") {
    thr <- w_threshold
  } else {
    v <- abs(as.vector(B_dag))
    v[diag_idx] <- NA_real_
    v <- v[is.finite(v)]
    if (quantile_ignore_zeros) v <- v[v > 0]
    if (length(v) == 0) {
      thr <- Inf
    } else {
      thr <- as.numeric(stats::quantile(v, probs = q_threshold, names = FALSE))
    }
  }

  B_dag[abs(B_dag) < thr] <- 0.0
  return(B_dag)
}
