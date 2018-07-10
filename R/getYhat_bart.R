#'Map fitted latents to categorical prediction 
#'
#'Map fitted latents to categorical prediction for the jth posterior draw,
#'@param j Index of the current round of posterior draws,
#'@param ndim Number of latents,
#'@param n Sample size,
#'@param releveled The index-category dictionary returned from multinomial probit BART fit,
#'@param maxy The index of the reference level from within the multinomial probit BART core,
#'@param mu Sum of trees fit for the unscaled latent variables in 2 dimensions - (n x ndraws) x ndim,
#'@param sigmas Posteriors of the covariance matrix in 2 dimensions - (ndim x ndim) x ndraws; the jth column is the covariance matrix estimated from jth round of posterior draw by column,
#'@param alpha The estimated working parameter for scaling the sum of trees fit, 
#'@return yhat Simulated outcome as a vector of length n,
#'@import MASS bayesm mlbench mlogit cvTools stats
#'@export

getYhat_bart = function(j, ndim, n, releveled, maxy, mu, sigmas, alpha){
  
  mu_j = mu[(n*(j-1)+1):(n*j),]
  mu_j = mu_j / alpha
  Sig_j = matrix(sigmas[,j],ncol = ndim)
  # generate mvnorm by mean equals to vector of zeros, then translate to mean at mu_j
  tmp = mvrnorm(n, rep(0, ndim), Sig_j)
  tmp = tmp + mu_j
  
  pclass = max.col(tmp)
  maxw = apply(tmp,1,max)
  pclass[which(maxw<0)] = maxy
  
  yhat =  releveled[pclass] #vector of length n
  return(yhat)
}