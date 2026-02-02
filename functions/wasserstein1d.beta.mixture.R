wasserstein1d.beta.mixture <- function(F_params, G_params,
                                       rel.tol = 1e-6, abs.tol = 1e-6) {
  is_single <- function(params) length(params) == 2
  
  f <- function(x) {
    pF <- if (is_single(F_params)) {
      pbeta(x, F_params[1], F_params[2])
    } else {
      pmix.beta(x,
                F_params[1:m0],
                F_params[(m0 + 1):(2 * m0)],
                F_params[(2 * m0 + 1):(3 * m0)])
    }
    
    pG <- if (is_single(G_params)) {
      pbeta(x, G_params[1], G_params[2])
    } else {
      pmix.beta(x,
                G_params[1:m0],
                G_params[(m0 + 1):(2 * m0)],
                G_params[(2 * m0 + 1):(3 * m0)])
    }
    
    abs(pF - pG)
  }
  
  tryCatch(
    integrate(f, lower = 0, upper = 1,
              rel.tol = rel.tol, abs.tol = abs.tol,
              stop.on.error = FALSE)$value,
    error = function(e) NA_real_
  )
}