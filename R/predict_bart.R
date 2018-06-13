#'Predict Method for BART Fits
#'
#'Predicted values based on Bayesian Additive Regression Trees model object,
#'@param obj Fitted model object from BART_mod,
#'@param newdata The data matrix to look for variables with which to predict, typically without the response; if Gcomp equals TRUE, number of rows should equal n x ndraws, where n is the number of subject. For the ith subject, (i, i+n, ..., i+(ndraws-1)n) rows are its simulated posterior outcomes from the previous simulation,
#'@param Gcomp Make predictions in the format of dynamic G-computation if TRUE; the default is FALSE,
#'@return treefit ndraws x n posterior matrix of the sum of trees fit,
#'@return samp_y ndraws x n posterior matrix of the simulated outcome,
#'@examples
#'##simulate data (example from Friedman MARS paper)
#'f = function(x){
#'  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
#'}
#'sigma = 1.0 #y = f(x) + sigma*z , z~N(0,1)
#'n = 100 #number of observations
#'set.seed(99)
#'x=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
#'newx=matrix(runif(n*10),n,10) 
#'Ey = f(x)
#'y=Ey+sigma*rnorm(n)
#'dat = data.frame(x,y)
#'fml = as.formula("y ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10")
#'bfit = BART_call(fml, data = dat, test.data = NULL,
#'                 Prior = list(nu = 3, sigq = 0.9,
#'                              ntrees = 100,
#'                              kfac = 2,
#'                              pswap = 0.1, pbd = 0.5, pb = 0.25,
#'                              alpha = 0.95, beta = 2.0,
#'                              nc = 100, minobsnode = 10),
#'                 Mcmc = list(burn=100, ndraws = 1000))
#'
#'newy = predict_bart(bfit, newx)
#'@import bayesm mlbench mlogit cvTools stats
#'@export
#'@useDynLib GcompBART
predict_bart  <- function(obj, newdata = NULL, Gcomp = FALSE)
{
  if(is.data.frame(newdata)) newdata = as.matrix(newdata)
  
  xcolnames = obj$xcolnames
  
  if (!is.null(newdata)){
    if(length(xcolnames) == 1 ){
      Xtest <- data.frame(newdata[,xcolnames])
      names(Xtest) <- xcolnames[1]
    } else {
      Xtest <- newdata[,xcolnames]
    }
  } else {
    stop("newdata can not be NULL.")
  }
  
  testn = nrow(Xtest)
  
  ntrees = obj$ntrees
  burn = obj$burn
  ndraws = obj$ndraws
      
  sigmasample = obj$sigmasample #on the scale of original data
  
  if(Gcomp){
    testn = testn / ndraws
    if(testn %% 1 != 0)
      stop(paste("Error: testn does not equal to n x ndraws."))
  }
  
  
  cat("Number of trees: ", ntrees, ".\n\n", sep="")
  cat("Number of draws: ", ndraws, ".\n\n", sep="")
  cat("burn-in: ", burn, ".\n\n", sep="")
  
  L1 = obj$TreeMod[[1]]
  L2 = obj$TreeMod[[2]]
  L3 = obj$TreeMod[[3]]
  L4 = obj$TreeMod[[4]]
  L5 = obj$TreeMod[[5]]
  L6 = obj$TreeMod[[6]]
  L7 = obj$TreeMod[[7]]
  L8 = obj$TreeMod[[8]]
  
  xi =  matrix(unlist(obj$xi), ncol = length(obj$xi[[1]]), byrow = TRUE)
  
  if(obj$type == "continuous"){
    rgy = obj$rgy
    
    res =   mybartpred(Gcomp,
                       L1,
                       L2,
                       L3,
                       L4,
                       L5,
                       L6,
                       L7,
                       L8,
                       Xtest,
                       as.integer(ncol(Xtest)),
                       as.integer(testn),
                       as.integer(ndraws),
                       as.integer(burn),
                       as.integer(ntrees),
                       xi)
    
    #vec_test is the sum of tree fits
    
    treefit = (rgy[2]-rgy[1])*(res$vec_test+.5) + rgy[1]
    stmp = rep(sigmasample,each = testn)
    samp_y = rnorm(testn*ndraws, treefit, stmp)
    samp_y = matrix(samp_y, nrow = testn)
    treefit = matrix(treefit, nrow = testn)

  } else if(obj$type == "binary"){
    
    res =   mybartpred(Gcomp,
                       L1,
                       L2,
                       L3,
                       L4,
                       L5,
                       L6,
                       L7,
                       L8,
                       Xtest,
                       as.integer(ncol(Xtest)),
                       as.integer(testn),
                       as.integer(ndraws),
                       as.integer(burn),
                       as.integer(ntrees),
                       xi)
    treefit = res$vec_test
    samp_y = rnorm(testn*ndraws, treefit, 1)
    samp_y = 1*(samp_y > 0)
    samp_y = matrix(samp_y, nrow = testn)
       
    treefit = matrix(treefit, nrow = testn)
    
  } else {
    
    treefit = c()
    
    for(i in 1:(obj$ndim)){
      Lk1 = obj$TreeMod[[1]][[i]]
      Lk2 = obj$TreeMod[[2]][[i]]
      Lk3 = obj$TreeMod[[3]][[i]]
      Lk4 = obj$TreeMod[[4]][[i]]
      Lk5 = obj$TreeMod[[5]][[i]]
      Lk6 = obj$TreeMod[[6]][[i]]
      Lk7 = obj$TreeMod[[7]][[i]]
      Lk8 = obj$TreeMod[[8]][[i]]
      
      res =   mybartpred(Gcomp,
                         Lk1,
                         Lk2,
                         Lk3,
                         Lk4,
                         Lk5,
                         Lk6,
                         Lk7,
                         Lk8,
                         Xtest,
                         as.integer(ncol(Xtest)),
                         as.integer(testn),
                         as.integer(ndraws),
                         as.integer(burn),
                         as.integer(ntrees),
                         xi)
      treefit = cbind(treefit, res$vec_test)
    }
    
    sigmas = matrix(obj$sigmasample, nrow = (obj$ndim)^2, ncol = testn)
    
    samp_y = lapply(1:ndraws, function(j) getYhat_bart(j, obj$ndim, testn, obj$releveled, obj$maxy, treefit, sigmas))

    samp_y = simplify2array(samp_y)
    #tmp = c(t(treefit)) # z_draw_id_dim: e.g. dim = 2, (z111,z112,z121,z122,..., z nd n 1, z nd n 2)
    
    #res =   mympbartpred(as.integer(obj$ndim),
    #                   as.integer(testn),
    #                   as.integer(ndraws),
    #                   tmp, 
    #                   obj$sigmasample, 
    #                   obj$working_param,
    #                   as.integer(obj$maxy))
    #samp_y = matrix(res$samp_y, nrow = testn)
    
  }
    
  
  ret = list(treefit = treefit,
             samp_y = samp_y);
  
  return(ret)
  
}
