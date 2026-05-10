pmle.beta.rand <- function(x, m0, n.init = 5, n.iter = 10, max.iter = 5000, tol = 1e-6,
                           epsilon = 1, an = NULL) {
  if (is.null(an)) an <- length(x)^(1/2)
  
  all_init_params <- NULL
  for (i in seq_len(2)) {
    mix_prop <- runif(m0)
    mix_prop <- mix_prop / sum(mix_prop)
    alpha <- exp(rnorm(m0))
    beta <- exp(rnorm(m0))
    all_init_params <- rbind(all_init_params, c(mix_prop, alpha, beta))
  }
  
  for (i in seq_len(n.init)) {
    group_assign <- sample(1:m0, size = length(x), replace = TRUE)
    mix_prop <- as.numeric(table(factor(group_assign, levels = 1:m0))) / length(x)
    shapes_sam <- sam.calculation(x, group_assign, m0)
    all_init_params <- rbind(all_init_params, c(mix_prop, shapes_sam))
  }
  
  for (i in seq_len(n.init)) {
    group_assign <- sample(1:m0, size = length(x), replace = TRUE)
    mix_prop <- as.numeric(table(factor(group_assign, levels = 1:m0))) / length(x)
    shapes_sar <- sar.calculation(x, group_assign, m0)
    all_init_params <- rbind(all_init_params, c(mix_prop, shapes_sar))
  }
  
  results <- vector("list", nrow(all_init_params))
  
  for (i in seq_len(nrow(all_init_params))) {
    para0 <- all_init_params[i,]
    out <- NULL
    for (j in seq_len(n.iter)) {
      output <- pmle.beta.sub(x, m0, para0, an, epsilon)
      para0 <- output[1:(3 * m0)]
      out <- output
    }
    results[[i]] <- out
  }
  
  candidates <- do.call(rbind, results)
  index <- which.max(candidates[, 3 * m0 + 2])
  para0 <- candidates[index, 1:(3 * m0)]
  ploglik <- candidates[index, 3 * m0 + 2]
  
  increment <- Inf
  tt <- 0
  
  while (increment > tol && tt < max.iter) {
    step <- pmle.beta.sub(x, m0, para0, an, epsilon)
    para0 <- step[1:(3 * m0)]
    increment <- step[3 * m0 + 2] - ploglik
    ploglik <- step[3 * m0 + 2]
    tt <- tt + 1
  }
  
  mix_prop <- para0[1:m0]
  alpha <- para0[(m0 + 1):(2 * m0)]
  beta <- para0[(2 * m0 + 1):(3 * m0)]
  
  pdf.sub <- matrix(NA_real_, nrow = m0, ncol = length(x))
  for (k in seq_len(m0)) pdf.sub[k, ] <- mix_prop[k] * dbeta(x, alpha[k], beta[k])
  
  pdf.mixture <- colSums(pdf.sub) + 1e-100
  ww <- sweep(pdf.sub, 2, pdf.mixture, "/")
  classification <- apply(t(ww), 1, which.max)
  
  loglik <- sum(log(dmix.beta(x, mix_prop, alpha, beta) + 1e-100))
  
  list(
    mix_prop = rousignif(mix_prop),
    alpha = rousignif(alpha),
    beta = rousignif(beta),
    loglik = rousignif(loglik),
    ploglik = rousignif(ploglik),
    iter.n = tt,
    classification = classification
  )
}