# GcompBART
Accelerated Multinomial Probit BART and an Integrated Tool  for Bayesian Prediction and Causal Comparison of Dynamic Policies

This package combines two work:

1. Bayesian Framework for Predictive and Causal Modeling with Application to HIV care Cascade 

This work studies the Bayesian formulation of G computation algorithm, and incorporates Bayesian additive regression trees (BART) into the generative components of the causal framework for posterior sampling of counterfactual longitudinal paths over time under certain dynamic policies of interest.
This integrated tool can be used when the time-varying confounders and the outcomes are continuous, binary, or multinomial. 

2. Accelerated Multinomial Probit Bayesian Additive Regression Trees (MPBART)

The model for multinomial response uses the multinomial probit (MNP) regression framework (Imai and van Dyk 2005). In this package I improved the original MPBART (Kindo et al 2016) so the MCMC convergence is significantly faster.


** For  OS X 10.11 and higher, R version 4.3.0 and higher, install gfortran-12.2-universal.pkg  from https://cran.r-project.org/bin/macosx/tools/ **

```
install.packages("devtools") # if you have not installed "devtools" package
devtools::install_github("yizhenxu/GcompBART")

library(GcompBART)

nb = 100 # number of burn-in samples
nd = 1000 # number of post-burn-in posterior draws 
nt = 100 # nubmer of trees
```

# simulate binomial/continuous data

** simulate data for continuous outcome, binary outcome, and covariates (example from Friedman MARS paper)**

```  
f = function(x){
  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma = 1.0 #y = f(x) + sigma*z , z~N(0,1)
n = 100 #number of observations

set.seed(99)

x = matrix(runif(n*10),n,10) #10 variables, only first 5 matter
Ey = f(x)
y = Ey+sigma*rnorm(n)
z = 1*(y > 10)
dat = data.frame(x,y,z)
```

# continuous outcome model

```
Prior_cont = list(nu = 3, sigq = 0.9,
                  ntrees = nt,
                  kfac = 2,
                  pswap = 0.1, pbd = 0.5, pb = 0.25,
                  alpha = 0.95, beta = 2.0,
                  nc = 100, minobsnode = 10)
Mcmc_cont = list(burn=nb, ndraws = nd)

fml = "y ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10"

bfit_con = model_bart(as.formula(fml), data = dat, type = "continuous",
                      Prior = Prior_cont,
                      Mcmc = Mcmc_cont)

pred_con = predict_bart(obj = bfit_con, newdata = dat) # prediction using new data
mean((pred_con$treefit - dat$y)^2) # mean squared error of using the sum of trees as predictions
```

# binary outcome model

```
Prior_binary = list(ntrees = nt,
                    kfac = 2,
                    pswap = 0.1, pbd = 0.5, pb = 0.25,
                    alpha = 0.95, beta = 2.0,
                    nc = 100, minobsnode = 10)
Mcmc_binary = list(burn=nb, ndraws = nd)

fml = "z ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10"

bfit_bin = model_bart(as.formula(fml), data = dat, type = "binary",
                      Prior = Prior_binary,
                      Mcmc = Mcmc_binary)

pred_bin = predict_bart(obj = bfit_bin, newdata = dat) # prediction using new data
mean(pred_bin$samp_y != dat$z) # misclassification rate of the binary predictions (randomly generated based on the estimated probabilities)
mean(1*(pred_bin$treefit>0) != dat$z) #misclassification rate of the 1(sum of trees > 0) as binary predictions
```

# simulate multinomial data

** simulate data for covariates and multinomial outcome (simulation setting 2 of the Accelerated MPBART manuscript) **

```
n = 1000 # sample size of the training + test sets
ntr = 500 # sample size of the training set
f1 = function(x){
  15*sin(pi*x[,1]*x[,2]) + (x[,3]-.5)^2 - 10*x[,4] - 5*x[,5]
}

f2 = function(x,v){
  (x[,3]-.5)^2 - x[,4]*x[,5] + 4*v
}

set.seed(2024)
u=matrix(runif(n*5),n,5) #10 variables, only first 5 matter
v = runif(n)*2

Sig = matrix(c(1, 0.5, 0.5, 1), 2, 2)

z1 = f1(u); z2 = f2(u,v)
library(MASS)
z = mvrnorm(n, c(0,0), Sig)
z = z + cbind(z1, z2)
y = 1*(z[,2] >= z[,1]) +1
y[z[,1] < 0 & z[,2] < 0] = 3
table(y)

newdat = data.frame(u,v,y)
colnames(newdat)[6] = "X6"
trd = newdat[1:ntr,]
ted = newdat[(ntr+1):n,]
table(trd$y)
table(ted$y)
```

# multinomial outcome model

```
Prior_mult = function(p, ntree){
  return(list(nu = p-1+1, V = diag(p-1),
              ntrees = ntree,
              kfac = 2,
              pswap = 0.1, pbd = 0.5, pb = 0.25,
              alpha = 0.95, beta = 2.0,
              nc = 100, minobsnode = 10))
}

Mcmc_mult = function(p,w0 = NULL,sig0 = NULL,nb, nd){
  res = list(sigma0 = diag(p-1), burn = nb, ndraws = nd, nSigDr = 50)  
  if(!is.null(w0)){
    res = append(res, list("w0" = w0), length(res))
  } 
  if(!is.null(sig0)){
    res$sigma0 = sig0
  }
  return(res)
}

fml = "y ~ X1 + X2 + X3 + X4 + X5 + X6"

baseyl = 2 # reference level

p=length(unique(trd$y)) # number of outcome levels

bfit_mult = model_bart(as.formula(fml), data = trd, type = "multinomial",
                   base = baseyl,
                   Prior = Prior_mult(p = p, ntree = nt),
                   Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                   correction = FALSE, Kindo = FALSE, do_alpha2_prior = FALSE)

```

** calculate the training and test set predictive accuracy **

```
ymode = function(vec){
  names(which.max(table(vec))) 
}

accuracyfun = function(bmpfit,testd,trdy,tedy){
  # posterior prediction by the entire distribution
  # mode accuracy, mean accuracy
  res = rep(NA, 4); names(res) = c("Train Mean", "Train Mode", "Test Mean", "Test Mode")
  nd = bmpfit$ndraws
  # train prediction
  BTr = bmpfit$samp_y
  yTr = matrix(rep(trdy, nd),ncol = nd)
  PmodeTr = apply(BTr, 1, ymode)
  res[1] = mean(BTr == yTr)
  res[2] = mean(PmodeTr == trdy)
  # test prediction
  BTe = predict_bart(obj = bmpfit, newdata = testd)$samp_y
  yTe = matrix(rep(tedy, nd),ncol = nd)
  PmodeTe = apply(BTe, 1, ymode)
  res[3] = mean(BTe == yTe)
  res[4] = mean(PmodeTe == tedy)
  return(res)
}

accuracyfun(bfit_mult, ted, trd$y, ted$y)
#Train Mean Train Mode  Test Mean  Test Mode 
#0.873856   0.950000   0.828820   0.890000 
```
