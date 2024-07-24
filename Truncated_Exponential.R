# 05/08/2019
# Truncated Exponential Distribution

# 1. a function to generate samples from Exp(rate) truncated at (a,b)
sam_texp <- function(n, rate, a=0, b=Inf){
  quans = runif(n)
  sams = -(log(exp(-rate*a)-quans*(exp(-rate*a)-exp(-rate*b))))/rate
  return(sams)
}


# 2. a function to evaluate density of truncated exponential
dtexp <- function(x, rate, a=0, b=Inf, log=F){
  if(any(x < a | x > b)){
    stop("Truncated Exponential: input values not within the bounds!\n")
  }
  dx = rate * exp(-rate * x)/(exp(-rate * a) - exp(-rate * b))
  if(length(x)==1){
    res = dx
  }else{
    res = prod(dx)
  }
  if(log){
    return(log(res))
  }else{
    return(res)
  }
}
sam_texp2 = function(n, rate, a = 0, b = 0, c = Inf){
  quans = runif(n)
  p1 = exp(-rate[, 1]*a) - exp(-rate[, 1]*b)
  p2 = exp(-rate[, 2]*b) - exp(-rate[, 2]*c)
  p_sum = p1 + p2
  p = cbind(p1, p2)/p_sum
  sams = rep(0, n)
  idx1 = which(quans <= p[, 1])
  idx2 = which(quans > p[, 1])
  sams[idx1] = -log(exp(-rate[idx1, 1]*a[idx1]) - quans[idx1]*p_sum[idx1])/rate[idx1, 1]
  sams[idx2] = -log(exp(-rate[idx2, 2]*b[idx2]) - (quans[idx2] - p[idx2, 1])*p_sum[idx2])/rate[idx2, 2]
  return(sams)
}
sam_texp3 = function(n, rate, a = 0, b = 0, c = 0, d = Inf){
  quans = runif(n)
  p1 = exp(-rate[, 1]*a) - exp(-rate[, 1]*b)
  p2 = exp(-rate[, 2]*b) - exp(-rate[, 2]*c)
  p3 = exp(-rate[, 3]*c) - exp(-rate[, 3]*d)
  p_sum = p1 + p2 + p3
  p = cbind(p1, p2, p3)/p_sum
  sams = rep(0, n)
  idx1 = which(quans <= p[, 1])
  idx2 = which(p[, 1] < quans & quans <= (p[, 1]+p[, 2]))
  idx3 = which((p[, 1]+p[, 2]) < quans & quans <= 1)
  sams[idx1] = -log(exp(-rate[idx1, 1]*a[idx1]) - quans[idx1]*p_sum[idx1])/rate[idx1, 1]
  sams[idx2] = -log(exp(-rate[idx2, 2]*b[idx2]) - (quans[idx2] - p[idx2, 1])*p_sum[idx2])/rate[idx2, 2]
  sams[idx3] = -log(exp(-rate[idx3, 3]*c[idx3]) - (quans[idx3] - (p[idx3, 1] + p[idx3, 2]))*p_sum[idx3])/rate[idx3, 3]
  return(sams)
}


# # checking
# x = sam_texp(20,0.5,1,7)
# log(prod(TruncatedDistributions::dtexp(x,0.5,1,7)))
# dtexp(x,0.5,1,7,log = T)
# # seems to be working
