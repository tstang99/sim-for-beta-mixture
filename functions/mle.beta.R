mle.beta.sub <- function(x, m0, para0, epsilon) {
  mix_prop <- para0[1:m0]
  alpha <- para0[(m0 + 1):(2 * m0)]
  beta <- para0[(2 * m0 + 1):(3 * m0)]
  theta <- para0[(m0 + 1):(3 * m0)]  
  
  n <- length(x)
  
  # E-step
  pdf.sub <- t(mapply(function(pp, aa, bb) pp * dbeta(x, aa, bb),
                      mix_prop, alpha, beta))
  pdf.mixture <- colSums(pdf.sub) + 1e-100
  ww <- sweep(pdf.sub, 2, pdf.mixture, FUN = "/")
  
  mix_prop <- (rowSums(ww) + epsilon) / (n + m0 * epsilon)
  
  # M-step
  theta <- M.Step(x, t(ww), mix_prop, theta)
  alpha <- theta[1:m0]
  beta <- theta[(m0 + 1):(2 * m0)]
  
  dens <- dmix.beta(x, mix_prop, alpha, beta)
  loglike <- sum(log(dens + 1e-100))
  
  ind <- order(alpha)
  
  c(mix_prop[ind], alpha[ind], beta[ind], loglike)
}

M.Step <- function(x, ww, mix_prop, theta) {
  k <- length(mix_prop)
  
  loglikelihood <- function(theta, x) {
    if (any(theta <= 0)) return(NA)
    alpha <- theta[1:k]
    beta <- theta[(k + 1):(2 * k)]
    dens_mat <- sapply(1:k, function(j) dbeta(x, alpha[j], beta[j]))
    ll <- sum(ww * log(dens_mat))
    ll
  }
  
  gradlik <- function(theta, x) {
    if (any(theta <= 0)) return(NA)
    alpha <- theta[1:k]
    beta <- theta[(k + 1):(2 * k)]
    
    g1j <- sapply(1:k, function(j) {
      sum(ww[, j] * (log(x) + digamma(alpha[j] + beta[j]) 
                     - digamma(alpha[j])))
    })
    
    g2j <- sapply(1:k, function(j) {
      sum(ww[, j] * (log(1 - x) + digamma(alpha[j] + beta[j]) 
                     - digamma(beta[j])))
    })
    c(g1j, g2j)
  }
  
  hesslik <- function(theta, x) {
    if (any(theta <= 0)) return(NA)
    alpha <- theta[1:k]
    beta <- theta[(k + 1):(2 * k)]
    Hess <- matrix(0, nrow = 2 * k, ncol = 2 * k)
    for (j in 1:k) {
      Hess[j, j] <- sum(ww[, j]
                        * (trigamma(alpha[j] + beta[j]) 
                           - trigamma(alpha[j]))) 
      Hess[j + k, j + k] <- sum(ww[, j]
                                * (trigamma(alpha[j] + beta[j]) 
                                   - trigamma(beta[j]))) 
      Hess[j, j + k] <- sum(ww[, j] 
                            * trigamma(alpha[j] 
                                       + beta[j]))
      Hess[j + k, j] <- Hess[j, j + k]
    }
    Hess
  }
  
  new_est_param <- maxLik::maxNR(fn = loglikelihood,
                                 grad = gradlik,
                                 hess = hesslik,
                                 start = theta,
                                 x = x)
  new_est_param$estimate
}


mle.beta <- function(x, m0, n.init=5, n.iter = 10, 
                     max.iter = 5000, tol = 1e-8,
                     epsilon = 1, maxit = 5000) {
  an <- length(x)^(0.5)
  sar_init <- sar.beta.mix(x,m0,epsilon=epsilon,an=an)
  sam_init <- sam.beta.mix(x,m0,epsilon=epsilon,an=an)
  all_init_params <- rbind(c(sar_init$mix_prop,
                             sar_init$alpha,sar_init$beta),
                           c(sam_init$mix_prop,
                             sam_init$alpha,sam_init$beta))
  for (i in seq_len(n.init)) {
    group_assign <- sample(1:m0, size = length(x), replace = TRUE)
    mix_prop <- as.numeric(table(factor(
      group_assign, levels = 1:m0))) / length(x)
    shapes_sam <- sam.calculation(x, group_assign, m0)
    all_init_params <- rbind(all_init_params, c(mix_prop, shapes_sam))
  }
  for (i in seq_len(n.init)) {
    group_assign <- sample(1:m0, size = length(x), replace = TRUE)
    mix_prop <- as.numeric(table(factor(
      group_assign, levels = 1:m0))) / length(x)
    shapes_sar <- sar.calculation(x, group_assign, m0)
    all_init_params <- rbind(all_init_params, c(mix_prop, shapes_sar))
  }
  results <- vector("list", nrow(all_init_params))
  for (i in seq_len(nrow(all_init_params))) {
    para0 <- all_init_params[i,]
    out <- NULL
    for (j in seq_len(n.iter)) {
      output <- mle.beta.sub(x, m0, para0, epsilon)
      para0 <- output[1:(3 * m0)]
      out <- output
    }
    results[[i]] <- out
  }
  candidates <- do.call(rbind, results)
  index <- which.max(candidates[, 3 * m0 + 1])
  para0 <- candidates[index, 1:(3 * m0)]
  loglike0 <- candidates[index, 3 * m0 + 1]
  increment <- Inf
  tt <- 0
  while (increment > tol && tt < max.iter) {
    step <- mle.beta.sub(x, m0, para0, epsilon)
    para0 <- step[1:(3 * m0)]
    loglike1 <- step[3 * m0 + 1]
    increment <- loglike1 - loglike0
    loglike0 <- loglike1
    tt <- tt + 1
  }
  mix_prop <- rousignif(para0[1:m0])
  mix_prop <- rousignif(mix_prop/sum(mix_prop))
  alpha <- rousignif(para0[(m0 + 1):(2 * m0)])
  beta <- rousignif(para0[(2 * m0 + 1):(3 * m0)])
  pdf.sub <- matrix(NA_real_, nrow = m0, ncol = length(x))
  for (k in seq_len(m0)) 
    pdf.sub[k, ] <- mix_prop[k] * dbeta(x, alpha[k], beta[k])
  pdf.mixture <- colSums(pdf.sub) + 1e-100
  ww <- sweep(pdf.sub, 2, pdf.mixture, "/")
  list(
    mix_prop = unname(mix_prop),
    alpha = unname(alpha),
    beta = unname(beta),
    loglik = unname(loglike0),
    iter.n = tt,
    classification = apply(t(ww), 1, which.max)
  )
}