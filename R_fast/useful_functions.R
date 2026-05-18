#' Classic expit function.
expit <- function(x) {
  return(exp(x) / (1 + exp(x)))
}

expit2 <- function(x) {
  return(1 / (1 + exp( - x)))
}

logit <- function(x) {
  return(log(x / (1 - x)))
}


log_dtrunc01 <- function(x, mean, sd, logZ) {
  z <- (x - mean) / sd
  dnorm(z, log = TRUE) - log(sd) - logZ
}

