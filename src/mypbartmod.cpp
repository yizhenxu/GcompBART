#include <Rcpp.h>
#include <iostream>
#include <string.h>
#include <stddef.h>
#include <stdio.h>
#include <vector>
#include <ctime>
#include <algorithm>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include "funs.h"
#include "tree.h"
#include "info.h"
using namespace Rcpp;

// [[Rcpp::export]]
List mypbartmod(NumericMatrix XMat,
               int nn, NumericVector y,
               int n_cov, 
               int nd,
               int burn,
               int numtrees,
               double kfac,
               double pswap,
               double pbd,
               double pb,
               double alpha,
               double beta,
               int nc,
               int minobs,
               NumericVector binaryX,
               int dgn) {
  
  size_t m = (size_t) numtrees;
  
  dinfo di;
  di.n_samp = (size_t)nn; di.n_cov = (size_t)n_cov; di.n_dim = 0;
  
  xinfo xi; /* The cutpoint matrix (ncov x nc) */
  getcutpoints(nc, (int)di.n_cov, (int)di.n_samp, XMat, xi);
  
  size_t minobsnode = (size_t) minobs;
  
  std::vector<double> allfit; /* The sum of fit of all trees for train.data (nsub x 1)*/
  allfit.resize(di.n_samp);

  std::vector<double> w; /* The latent variable for train.data (nsub x 1)*/
  w.resize(di.n_samp);
  
  for(size_t i=0;i<di.n_samp;i++) {
    allfit[i] = 0.0;
    w[i] = 2*(y[i]) - 1.0; //y is 0 and 1, w is -1 and 1
  }
  
  double u,Z; //for updating latent w
  
  /* fit of current tree */
  //std::vector<double> ftemp; /* for allfit  (nsub x 1)*/
  //ftemp.resize(di.n_samp);
  
  std::vector<std::vector<double> > ftemp; /*m x nsub, (i,j) is the current fit of tree i person j*/
  ftemp.resize(m);  
  for(size_t i=0; i<m; i++){
    ftemp[i].resize(di.n_samp);
    for(size_t j=0; j<di.n_samp; j++){
      ftemp[i][j] = 0;
    }
  }  
  
  
  //partial residuals from backfitting
  std::vector<double> rtemp;
  rtemp.resize(di.n_samp);
  
  pinfo pi;
  pi.pswap = pswap; //prob of swap move, default 0.1
  pi.pbd = pbd; //prob of birth/death move, default 0.5
  pi.pb = pb; //prob of birth given  birth/death, default 0.25
  
  pi.alpha = alpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
  pi.beta = beta; //
  pi.tau=3.0/(kfac*sqrt((double)m)); // categorical outcome
  pi.sigma= 1; //error standard deviation of the latent W
  
  //initialize tree
  
  std::vector<tree>  t; /* ntree vector of trees */
  t.resize(m);
  for(size_t i=0; i<m; i++){
    t[i].setm(0.00);
  }

  //percA[i] = average over all percAtmp for draw i
  std::vector<double> percAtmp;
  percAtmp.resize(m);
  //numNodes[i] = average over all numNodestmp for draw i
  std::vector<double> numNodestmp;
  numNodestmp.resize(m);
  //numLeaves[i] = average over all numLeavestmp for draw i
  std::vector<double> numLeavestmp;
  numLeavestmp.resize(m);
  //treeDepth[i] = average over all treeDepthtmp for draw i
  std::vector<double> treeDepthtmp;
  treeDepthtmp.resize(m);
  //nLtD[i] = number of leaves at tree depth of the ith tree for the current round of draw
  std::vector<double> nLtDtmp;
  nLtDtmp.resize(m);
  
  for(size_t i=0; i<m; i++){
    percAtmp[i] = 0.0;
    numNodestmp[i] = 1.0;
    numLeavestmp[i] = 1.0;
    treeDepthtmp[i] = 0.0;
    nLtDtmp[i] = 1.0;
  }
  
  //incProptmp[j] = count occurence of xj in splitting rules for the current round of draw
  //incProp[j] = average of incProptmp[j] over all draws
  std::vector<double> incProptmp;
  incProptmp.resize(di.n_cov);
  
  for(size_t i=0; i<di.n_cov; i++){
    incProptmp[i] = 0.0;
  }
  int tnd = burn + nd;
  NumericVector vec_train(nn*nd), percA(tnd), numNodes(tnd), numLeaves(tnd), treeDepth(tnd), incProp(n_cov);
  
  //TreeMod
  std::vector<double> L1; //action type
  L1.reserve(m*tnd);
  std::vector<double> L2; //node id when action type > 0
  L2.reserve(floor(0.5*m*tnd));
  std::vector<double> L3; //GROW info (v, c, mul, mur) when action type == 1
  L3.reserve(floor(4*0.5*pb*m*tnd));
  std::vector<double> L4; //PRUNE info (mu) when action type = =2
  L4.reserve(floor(0.5*(pbd - pb)*m*tnd));
  std::vector<double> L5; //SWAP info (lind, rind, vd, cd, vk, ck, nbnc.size()) when action type == 3
  L5.reserve(floor(7*0.5*pswap*m*tnd));
  std::vector<double> L6; //SWAP nx bottom nodes when action type == 3
  L6.reserve(floor(5*0.5*pswap*m*tnd));
  std::vector<double> L7; //CHANGE info (v, c, nbnv.size()) when action type == 4
  L7.reserve(floor(3*0.5*(1 - pbd - pswap)*m*tnd));
  std::vector<double> L8; //CHANGE nx bottom nodes when action type == 4
  L8.reserve(floor(5*0.5*(1 - pbd - pswap)*m*tnd));
  
  //MCMC
  
  time_t tp;
  int time1 = time(&tp);
  
  /* Initialize counters for outputs sigmasample, vec_test and vec_train */
  int countvectrain = 0;
  
  
  for(int loop=0;loop < tnd;loop++) { /* Start posterior draws */
  
  if(loop%100==0) Rprintf("\n iteration: %d of %d \n",loop, nd+burn);
  
  /* Step 1 */
  /* See tree sampling theory.doc for explanation */
  for(size_t ntree = 0 ; ntree <m; ntree++){
    //fit(t[ntree], XMat, di, xi, ftemp);
    
    for(size_t i=0;i<di.n_samp;i++) {
      allfit[i] -= ftemp[ntree][i];
      rtemp[i] = w[i] - allfit[i];
    }
    
    di.y = &rtemp[0];
    
    if(dgn){
      bd1(XMat, t[ntree], xi, di, pi, minobsnode, binaryX, &nLtDtmp[ntree], &percAtmp[ntree], &numNodestmp[ntree], &numLeavestmp[ntree], &treeDepthtmp[ntree], incProptmp, L1, L2, L3, L4, L5, L6, L7, L8);
    } else {
      bd(XMat, t[ntree], xi, di, pi, minobsnode, binaryX, L1, L2, L3, L4, L5, L6, L7, L8);
    }
    
    fit(t[ntree], XMat, di, xi, ftemp[ntree]);
    for(size_t i=0;i<di.n_samp;i++) {
      allfit[i] += ftemp[ntree][i];
    }
    
  }//ntree
  
  
  //done sampling (T,M)
  
  /* Step 2 update latent variable w*/
  for(size_t i=0;i<di.n_samp;i++) {
    u = unif_rand();
    if(y[i] > 0) {
      Z = R::qnorm((1.0-u)*R::pnorm(-allfit[i],0.0,1.0,1,0) + u,0.0,1.0,1,0);
    } else {
      Z = -R::qnorm((1.0-u)*R::pnorm(allfit[i],0.0,1.0,1,0) + u,0.0,1.0,1,0);
    }
    w[i] = allfit[i] + Z;
  }
  
  
  if(dgn){
    for(size_t i=0;i<m;i++) {
      percA[loop] += percAtmp[i];
      numNodes[loop] += numNodestmp[i];
      numLeaves[loop] += numLeavestmp[i];
      treeDepth[loop] += treeDepthtmp[i];
    }
    percA[loop] /= m;
    numNodes[loop] /= m;
    numLeaves[loop] /= m;
    treeDepth[loop] /= m;
  }
    
  if(loop>=burn){
    
    if(dgn){
      double numsplits = 0;
      for(size_t i=0; i<di.n_cov; i++){
        numsplits += incProptmp[i];
      }
      for(size_t i=0; i<di.n_cov; i++){
        incProp[i] += incProptmp[i]/numsplits;//proportion of all splitting rules that uses xi
      }
    }

    for(size_t k = 0; k <di.n_samp; k++){
      vec_train[countvectrain] = allfit[k];
      countvectrain++;
    }//end prediction for train
    
  }//end prediction for current loop
  
  } //end of loop
  
  if(dgn){
    for(size_t i=0; i<di.n_cov; i++){
      incProp[i] /= nd;
    }
  }
  
  int time2 = time(&tp);
  Rprintf("time for mcmc loop %d secs", time2-time1);
  
  for(size_t i=0; i<m; i++){
    t[i].tonull();
  }//delete trees
  
  List TreeMod = List::create(L1,L2,L3,L4,L5,L6,L7,L8);
  
  List z;
  
  if(dgn){
    z = List::create( Rcpp::Named("treefit") = vec_train, 
                      Rcpp::Named("Percent_Acceptance") = percA, 
                      Rcpp::Named("Tree_Num_Nodes") = numNodes, 
                      Rcpp::Named("Tree_Num_Leaves") = numLeaves, 
                      Rcpp::Named("Tree_Depth") = treeDepth, 
                      Rcpp::Named("Inclusion_Proportions") = incProp,
                      Rcpp::Named("TreeMod") = TreeMod,
                      Rcpp::Named("xi") = Rcpp::wrap(xi)) ;
  } else {
    z = List::create( Rcpp::Named("treefit") = vec_train,
                      Rcpp::Named("TreeMod") = TreeMod,
                      Rcpp::Named("xi") = Rcpp::wrap(xi)) ;
  }
  
  
  return z ;
  
}
