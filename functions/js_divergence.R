js_divergence <- function(p, q) {
  m  <- (p + q) / 2
  kl <- function(x, y) ifelse(x == 0, 0, x * log(x / y, 2))
  0.5 * sum(kl(p, m)) + 0.5 * sum(kl(q, m))
}