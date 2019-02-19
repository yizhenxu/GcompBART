#'Diagnostic Plots
#'
#'Plot average M-H acceptance percentage, tree nodes, tree leaves, tree depth, and predictor inclusion proportions.
#'@param bfit Fitted object from model_bart call,
#'@param plot_type 0 for predictor inclusion proportions plot; 1 for diagnostic plots of average M-H acceptance percentage, tree nodes, tree leaves, and tree depth,
#'@param byrow Default is TRUE, each row of the diagnostic plot represents one latent; column if FALSE. This is only for multinomial response,
#'@return Percent_Acceptance Percent acceptance of Metropolis-Hastings proposals across the ntrees number of trees for each posterior draw after burn-in,
#'@return Tree_Num_Nodes Average number of tree nodes across the ntrees number of trees for each posterior draw after burn-in,
#'@return Tree_Num_Leaves Average number of leaves across the ntrees number of trees for each posterior draw after burn-in,
#'@return Tree_Depth Average tree depth across the ntrees number of trees for each posterior draw after burn-in,
#'@return Inclusion_Proportions Predictor inclusion frequencies. Smaller value of ntrees (such as 10, 20, 50, 100) is recommended for the purposes of variable selection.
#'@examples
#'##simulate data (example from Friedman MARS paper)
#'f = function(x){
#'  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
#'}
#'sigma = 1.0 #y = f(x) + sigma*z , z~N(0,1)
#'n = 100 #number of observations
#'set.seed(99)
#'x=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
#'Ey = f(x)
#'y=Ey+sigma*rnorm(n)
#'fml = as.formula("y ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10")
#'bfit = BART_call(fml, data = dd, test.data = NULL,
#'                 Prior = list(nu = 3, sigq = 0.9,
#'                              ntrees = 200,
#'                              kfac = 2,
#'                              pswap = 0.1, pbd = 0.5, pb = 0.25,
#'                              alpha = 0.95, beta = 2.0,
#'                              nc = 100, minobsnode = 10),
#'                 Mcmc = list(burn=100, ndraws = 88))
#'DiagPlot(bfit, 0)
#'@import stats
#'@export
DiagPlot <- function(bfit, plot_type, byrow = TRUE){

  response_type = bfit$type
  if(response_type %in% c("continuous","binary","multinomial")){
    nd = bfit$ndraws + bfit$burn  
  } else {
    nd = bfit$ndraws + bfit$burn + bfit$fitMNP
  }
  
  if(plot_type){ #diagnostic plots
    if(response_type %in% c("multinomial","multinomial2","multinomial3")){
      pm1 = bfit$ndim
      par(mfrow=c(pm1,4))
      if(byrow == FALSE){
        par(mfcol=c(4,pm1))
      }
      for(i in 1:pm1){
        scatter.smooth(1:nd, bfit$Percent_Acceptance[i,], lpars =
                         list(col = "red", lwd = 3, lty = 3)
                       , xlab = "MCMC Iteration", ylab = "% of Tree Acceptance",main = paste0("latent ",i))
        scatter.smooth(1:nd, bfit$Tree_Num_Nodes[i,], lpars =
                         list(col = "red", lwd = 3, lty = 3)
                       , xlab = "MCMC Iteration", ylab = "Average Num of Tree Nodes")
        scatter.smooth(1:nd, bfit$Tree_Num_Leaves[i,], lpars =
                         list(col = "red", lwd = 3, lty = 3)
                       , xlab = "MCMC Iteration", ylab = "Average Num of Tree Leaves")
        scatter.smooth(1:nd, bfit$Tree_Depth[i,], lpars =
                         list(col = "red", lwd = 3, lty = 3)
                       , xlab = "MCMC Iteration", ylab = "Average Tree Depth")
      }
      
      par(mfrow=c(1,1))
    } else {
      par(mfrow=c(2,2))
      scatter.smooth(1:nd, bfit$Percent_Acceptance, lpars =
                       list(col = "red", lwd = 3, lty = 3)
                     , xlab = "MCMC Iteration", ylab = "% of Tree Acceptance")
      scatter.smooth(1:nd, bfit$Tree_Num_Nodes, lpars =
                       list(col = "red", lwd = 3, lty = 3)
                     , xlab = "MCMC Iteration", ylab = "Average Num of Tree Nodes")
      scatter.smooth(1:nd, bfit$Tree_Num_Leaves, lpars =
                       list(col = "red", lwd = 3, lty = 3)
                     , xlab = "MCMC Iteration", ylab = "Average Num of Tree Leaves")
      scatter.smooth(1:nd, bfit$Tree_Depth, lpars =
                       list(col = "red", lwd = 3, lty = 3)
                     , xlab = "MCMC Iteration", ylab = "Average Tree Depth")
      par(mfrow=c(1,1))
    }
    
  } else { #inclusion proportions
    if(response_type %in% c("multinomial","multinomial2","multinomial3")){
      pm1 = bfit$ndim
      par(mfrow=c(1,pm1))
      
      for(i in 1:pm1){
        tmp = bfit$Inclusion_Proportions[i,]
        barplot(sort(tmp,decreasing = T), ylab = "Inclusion Proportion", las=2, main=paste0("latent ",i))
      }
      
      par(mfrow=c(1,1))
    } else {
      tmp = bfit$Inclusion_Proportions
      barplot(sort(tmp,decreasing = T), ylab = "Inclusion Proportion", las=2)
    }
    
  }#plot_type
  
}
