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
#' @param w_threshold Numeric. A threshold below which edge weights are set to zero in the final output.
#'   Default is \code{0.1}.
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
#' @importFrom expm expm
#' @importFrom stats optim
#' @export
notears_projection <- function(B_init,
                               lam = 0.0,
                               max_iter = 100,
                               h_tol = 1e-8,
                               rho_max = 1e16,
                               w_threshold = 0.1) {

  p <- nrow(B_init)

  # Helper: Convert flattened weights (pos/neg split) back to adjacency matrix
  get_adj <- function(w) {
    w_pos <- w[1:(p*p)]
    w_neg <- w[(p*p + 1):(2*p*p)]
    matrix(w_pos - w_neg, nrow = p, ncol = p)
  }

  # Helper: Calculate Squared Loss and Gradient
  calc_loss <- function(W) {
    D <- W - B_init
    loss <- 0.5 * sum(D^2)
    G_loss <- D
    list(loss = loss, G_loss = G_loss)
  }

  # Helper: Calculate Acyclicity Constraint h(W) and Gradient
  calc_h <- function(W) {
    # h(W) = tr(expm(W*W)) - p
    E <- expm::expm(W * W)
    h_val <- sum(diag(E)) - p
    G_h <- t(E) * W * 2
    list(h = h_val, G_h = G_h)
  }

  rho <- 1.0
  alpha <- 0.0

  # Augmented Lagrangian Objective Function
  fn_obj <- function(w) {
    W <- get_adj(w)
    l_res <- calc_loss(W)
    h_res <- calc_h(W)
    obj <- l_res$loss + 0.5 * rho * h_res$h^2 + alpha * h_res$h + lam * sum(w)
    return(obj)
  }

  # Gradient of the Objective Function
  gr_obj <- function(w) {
    W <- get_adj(w)
    l_res <- calc_loss(W)
    h_res <- calc_h(W)
    G_smooth <- l_res$G_loss + (rho * h_res$h + alpha) * h_res$G_h
    c(as.vector(G_smooth) + lam, as.vector(-G_smooth) + lam)
  }

  # Initialization
  w_est <- numeric(2 * p * p)
  h_val <- Inf

  # Bounds for L-BFGS-B (Non-negative for split variables)
  lower_b <- numeric(2 * p * p)
  upper_b <- rep(Inf, 2 * p * p)

  # Enforce zero diagonal (no self-loops allowed directly)
  diag_idx <- (0:(p-1)) * p + (1:p)
  lower_b[diag_idx] <- 0; upper_b[diag_idx] <- 0
  lower_b[diag_idx + p*p] <- 0; upper_b[diag_idx + p*p] <- 0

  # Optimization Loop
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
  B_dag[abs(B_dag) < w_threshold] <- 0.0
  return(B_dag)
}
