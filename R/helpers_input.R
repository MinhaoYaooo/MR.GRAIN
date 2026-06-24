.safe_inverse <- function(M, ridge = 1e-6) {
  out <- tryCatch(solve(M), error = function(e) NULL)
  if (!is.null(out)) return(out)

  out <- tryCatch(solve(M + ridge * diag(nrow(M))), error = function(e) NULL)
  if (!is.null(out)) return(out)

  return(MASS::ginv(M))
}


.resolve_S_indices <- function(S, Z) {
  q <- ncol(Z)
  z_names <- colnames(Z)
  S_idx <- vector("list", length(S))

  for (j in seq_along(S)) {
    Sj <- S[[j]]

    if (length(Sj) == 0) {
      S_idx[[j]] <- integer(0)
      next
    }

    if (is.factor(Sj)) {
      Sj <- as.character(Sj)
    }

    if (is.numeric(Sj)) {
      if (any(!is.finite(Sj))) {
        stop(sprintf("S[[%d]] contains non-finite numeric values.", j))
      }
      if (any(abs(Sj - round(Sj)) > 0)) {
        stop(sprintf("S[[%d]] must contain integer column indices if numeric.", j))
      }

      Sj_idx <- as.integer(round(Sj))

      if (any(Sj_idx < 1 | Sj_idx > q)) {
        stop(sprintf("S[[%d]] contains invalid SNP indices.", j))
      }

      S_idx[[j]] <- unique(Sj_idx)
      next
    }

    if (is.character(Sj)) {
      if (is.null(z_names)) {
        stop(sprintf(
          "S[[%d]] contains character SNP IDs, but colnames(Z) is NULL.",
          j
        ))
      }

      match_idx <- match(Sj, z_names)

      if (any(is.na(match_idx))) {
        missing_ids <- unique(Sj[is.na(match_idx)])
        stop(sprintf(
          "S[[%d]] contains SNP IDs not found in colnames(Z): %s",
          j, paste(missing_ids, collapse = ", ")
        ))
      }

      S_idx[[j]] <- unique(as.integer(match_idx))
      next
    }

    stop(sprintf(
      "S[[%d]] must contain either numeric column indices or character SNP IDs.",
      j
    ))
  }

  return(S_idx)
}
