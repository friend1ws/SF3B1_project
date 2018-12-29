library(VGAM)

bb_mloglLikelihood <- function(params, values) {
  
  p_alpha <- params[1]
  p_beta <- params[2]
  
  v_n <- values[[1]]
  v_k <- values[[2]]
  

  ML <- 0
  ML <- ML + sum(lgamma(v_n + 1) - lgamma(v_k + 1) - lgamma(v_n - v_k + 1))
  ML <- ML + sum(lgamma(p_alpha + v_k) + lgamma(p_beta + v_n - v_k) - lgamma(p_alpha + p_beta + v_n))
  ML <- ML + length(v_n) * (lgamma(p_alpha + p_beta) - lgamma(p_alpha) - lgamma(p_beta))
    
  return(-ML)
    
}
  

zibb_mloglLikelihood <- function(params, values) {
  
  p_alpha <- params[1]
  p_beta <- params[2]
  p_pi <- params[3]
  
  v_n <- values[[1]]
  v_k <- values[[2]]
  
  zero_ind <- which(v_k == 0)
  non_zero_ind <- setdiff(1:length(v_k), zero_ind)
  
  ML <- 0
  ML <- ML + length(non_zero_ind) * log(1 - p_pi) + sum( bb_loglikelihood(p_alpha, p_beta, v_n[non_zero_ind], v_k[non_zero_ind]) )
  ML <- ML + sum(log(p_pi + (1 - p_pi) * exp( bb_loglikelihood(p_alpha, p_beta, v_n[zero_ind], v_k[zero_ind])) ))
  
  return(-ML) 
  
}

zip_mlogLikelihood <- function(params, values) {
  
  p_lambda <- params[1]
  p_pi <- params[2]

  v_k <- values 
  
  zero_ind <- which(v_k == 0) 
  non_zero_ind <- setdiff(1:length(v_k), zero_ind)
  
  ML <- 0
  ML <- ML + length(non_zero_ind) * log(1 - p_pi) - length(non_zero_ind) * p_lambda + sum(v_k) * log(p_lambda)
  ML <- ML + length(zero_ind) * log(p_pi + (1 - p_pi) * exp(-p_lambda))
  
}

zib_mlogLikelihood <- function(params, values) {
  
  p_prob <- params[1]
  p_pi <- params[2]
  
  v_n <- values[[1]]
  v_k <- values[[2]]

  zero_ind <- which(v_k == 0) 
  non_zero_ind <- setdiff(1:length(v_k), zero_ind)
  
  ML <- 0
  ML <- ML + length(non_zero_ind) * log(1 - p_pi) + sum(v_k[non_zero_ind]) * log(p_prob) + sum(v_n[non_zero_ind] - v_k[non_zero_ind]) * log(1 - p_prob)
  ML <- ML + sum(log(p_pi + (1 - p_pi) * (1 - p_prob)^v_n[zero_ind]))
  
  return(-ML)

}

bb_loglikelihood <- function(p_alpha, p_beta, v_n, v_k) {
  
  ML <- 0
  ML <- ML + lgamma(v_n + 1) - lgamma(v_k + 1) - lgamma(v_n - v_k + 1)
  ML <- ML + lgamma(p_alpha + v_k) + lgamma(p_beta + v_n - v_k) - lgamma(p_alpha + p_beta + v_n)
  ML <- ML + lgamma(p_alpha + p_beta) - lgamma(p_alpha) - lgamma(p_beta)
  
  return(ML)
  
}

zibb_optim <- function(v_n, v_k) {
  
  cret <- constrOptim(c(5, 5, 0.2), zibb_mloglLikelihood, grad=NULL, 
                      ui = rbind(diag(3), -diag(3)), ci=c(0.1, 0.1, 0.01, -1000, -1000, -0.99), 
                      values = list(v_n, v_k), outer.iterations = 1000, outer.eps = 1e-8)
  return(cret)
}

bb_optim <- function(v_n, v_k) {
  
  cret <- constrOptim(c(5, 5), bb_mloglLikelihood, grad=NULL, 
                      ui = rbind(diag(2), -diag(2)), ci=c(0.1, 0.1, -1000, -1000), 
                      values = list(v_n, v_k), outer.iterations = 1000, outer.eps = 1e-8)
  return(cret)
}

zib_optim <- function(v_n, v_k) {
  
  cret <- constrOptim(c(0.5, 0.2), zib_mlogLikelihood, grad=NULL, 
                      ui = rbind(diag(2), -diag(2)), ci=c(0, 0.01, -1, -0.99), 
                      values = list(v_n, v_k), outer.iterations = 1000, outer.eps = 1e-8)
  return(cret)
  
}

get_prob_for_zibb <- function(params, s_n, s_k) {
  
  p_alpha <- params[1]
  p_beta <- params[2]
  p_pi <- params[3]
  
  if (s_k == 0) {
    return(p_pi + (1 - p_pi) * dbetabinom.ab(0, s_n, p_alpha, p_beta))
  } else {
    return((1 - p_pi) * dbetabinom.ab(s_k, s_n, p_alpha, p_beta))
  }
  
}

get_prob_for_bb <- function(params, s_n, s_k) {
  
  p_alpha <- params[1]
  p_beta <- params[2]
  
  dbetabinom.ab(s_k, s_n, p_alpha, p_beta)

}

get_prob_for_zib <- function(params, s_n, s_k) {
  
  p_prob <- params[1]
  p_pi <- params[2]
  
  if (s_k == 0) {
    return(p_pi + (1 - p_pi) * dbinom(0, s_n, p_prob))
  } else {
    return((1 - p_pi) * dbinom(s_k, s_n, p_prob))
  }

}
