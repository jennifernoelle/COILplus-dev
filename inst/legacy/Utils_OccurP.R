
library(truncnorm)
library(mvtnorm)

# Simulate from plant occurrences from prior
simulate_plants <- function(nS, nM, nP, mh_pprior_sd, pi_L = 0.5){
  
  
  # True interaction matrix
  L <- matrix(rbinom(nM*nP, 1, pi_L), nM, nP)
  
  # Detection probabilities
  p <- runif(nM, 0, 1)
  q <- runif(nP, 0, 1)
  pq <- outer(p, q)
  
  # Study focus: all studies are animal-focused
  Focus <- array(0, dim = c(nM, nP, nS))
  p_occur_others <- matrix(0, nM, nS) # probability of occurrence for other species type (mammals, so in this case we can let it be 0/1)
  for(s in 1:nS){
    # First focus animal
    these.mammals <- sample(1:nM, 2)
    Focus[these.mammals[1],,s] <- 1
    p_occur_others[these.mammals[1], s] <- 1
    # Second focus animal (sometimes)
    u <- runif(1)
    if(u > 0.85){
      Focus[these.mammals[2],,s] <- 1
      p_occur_others[these.mammals[2],s] <- 1
    }
  }
  
  occur_others <- p_occur_others # since this is binary, they're equivalent
  
  # Plant occurrences and occurrence probabilities
  probs = seq(0,1, 0.1) # just some random probs to sample from
  prior_probs <- matrix(sample(probs, size = nP*nS, replace = TRUE), nP, nS) # Prior occurrence probabilities to use as input and to simulate true probs
  Pi_p <- matrix(rtruncnorm(n = nP*nS,  mean = prior_probs, sd = mh_pprior_sd, a = 0, b = 1), nP, nS) # True occurrence probabilities for this species type
  O_p <- matrix(rbinom(nP*nS, 1, Pi_p), nP, nS) # True occurrence indicators
  
  # Observed interactions
  A <- array(data = NA, dim = c(nM, nP, nS))
  for(s in 1:nS){
    this.prob <- L * Focus[,,s] * pq *  matrix(1, nM, 1) %*% O_p[,s]
    A[,,s] <- matrix(rbinom(nM*nP, 1, this.prob), nM, nP)
  }
  
  # Replace prior occurrence probs and true occurrence prob with 1 for detected species
  detected.p <- t(apply(A,2, function(x) colSums(x) > 0)) # plants are rows, studies are columns
  detected.m <- t(apply(A,1, function(x) colSums(x) > 0)) # mammals are rows, studies are columns
  
  Pi_p[which(detected.p == 1)] <- 1
  prior_probs[which(detected.p == 1)] <- 1  
  
  # Randomly initialize current states to be stored and provided as input
  theta_curr <- array(rnorm(nP*nS*2,0,1), dim = c(nP,nS,2)) # Store in mcmc function: must retrieve from previous iters, Tricky because of probit link
  p_curr <- pnorm(theta_curr[,,1]) # the current state of occurence prob to start from, these should converge to the true values in Pi_p
  occ_curr <- ifelse(theta_curr[,,2] >0 , 1, 0) # the current state of occurrence ind to start from, these should converge to the true values in O_p
  
  ret <- list(q = q,p = p, A= A, L = L, Focus = Focus, Pi_p = Pi_p, O_p = O_p, p_occur_others = p_occur_others,
              prior_probs = prior_probs, theta_curr = theta_curr, p_curr = p_curr, 
              occ_curr = occ_curr, occur_others = occur_others)
  return(ret)
}


# Log likelihood to be used for updating occurrence indicators 
# pi_prop_st: all new proposal occurrence probs for this study
# occur_prop_st: all new proposal occurrence indicators for this study
# prod1: (1-p_i p_j)^LFO
loglik <- function(pi_prop_st, prod1, occur_prop_st){
  prob1 <- pi_prop_st * prod1 # Multiply by occurrence probabilities to get prob occurrence == 1
  prob0 <- 1 - pi_prop_st # Prob occurrence == 0
  prob1 <- prob1 / (prob1 + prob0) # Standardize
  this_prob <- prob1*occur_prop_st + (1-prob1)*(1-occur_prop_st) # p = p1 if occ = 1, (1-p1) if occ = 0
  log(this_prob) # this gives a vector length num_obs
}

# Faster, stable replacement for loglik()
fast_loglik <- function(pi_prop_st, prod1, occur_prop_st) {
  eps <- 1e-12
  pi  <- pmin(pmax(pi_prop_st, eps), 1 - eps)
  pr1 <- pmin(pmax(prod1,       eps), 1 - eps)
  
  a <- log(pi) + log(pr1)        # log numerator
  b <- log1p(-pi)                # log(1 - pi)
  
  d <- b - a
  log_p1 <- -log1p(exp(d))       # log(prob1)
  log_p0 <- -log1p(exp(-d))      # log(1 - prob1)
  
  occur_prop_st * log_p1 + (1 - occur_prop_st) * log_p0
}






