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
#include "bd.h"
using namespace Rcpp;

// [[Rcpp::export]]
List mympbartmod_Kindos(NumericMatrix& XMat1,
                  NumericMatrix& sigmai1,
                  NumericMatrix& V1,
                  int nu,
                  int nn, int nndim,
                  NumericVector& y1,
                  IntegerVector& ncov, 
                  int nd,
                  int burn,
                  IntegerVector& m,
                  int nSigDr,
                  double kfac,
                  double pswap,
                  double pbd,
                  double pb,
                  double alpha,
                  double beta,
                  int nc,
                  int minobs,
                  NumericVector& binaryX1,
                  int dgn,
                  int Jiao,
                  NumericMatrix& w1,
                  int fitMNP,
                  int bylatent) {
  
  int mm; /* The maximum number of trees */
  mm = R_NegInf;
  for(int i=0; i < nndim;i++){
    if(mm < m[i]) mm = m[i];
  }
  
  dinfo di;
  di.n_samp = (size_t)nn; di.n_dim = (size_t)nndim;
  di.n_cov = 0;
  
  std::vector<double> V, sigmai;
  
  V.resize(nndim*nndim);
  sigmai.resize(nndim*nndim);
  
  int itemp = 0;
  for(int i = 0; i < nndim; i++){
    for(int j = 0; j < nndim; j++){
      V[itemp] = V1(i,j);
      sigmai[itemp] = sigmai1(i,j);
      
      itemp++;
    }
  }
  
  std::vector<int> y;
  y.resize(nn);
  for(int i = 0; i < nn; i++){
    y[i] = y1[i];
  }
  
  
  /* allow different X for each dim */
  std::vector<std::vector<std::vector<double> > > XMat;
  XMat.resize(di.n_dim);
  for(int j=0; j < nndim; j++){
    XMat[j].resize(di.n_samp);
  }
  
  for(int j=0; j < nndim; j++){
    for(int i=0; i < nn; i++){
      XMat[j][i].resize(ncov[j]);
    }
  }
  
  itemp = 0;
  for(int j=0; j < nndim; j++){
    for(int k=0; k< ncov[j]; k++){
      for(int i=0; i < nn; i++){
        XMat[j][i][k] = XMat1(i,itemp);
      }
      itemp++ ;
    }
  }
  
  Rprintf("\n Finish preparing design matrix XMat. \n");
  
  std::vector<xinfo> xi;
  xi.resize(nndim);
  for(int j=0; j< nndim; j++){
    getcutpoints(nc, ncov[j], nn, XMat[j], xi[j]);
  }
  
  /*define binaryX as a list of numericvectors*/
  std::vector<std::vector<double> > binaryX;
  binaryX.resize(nndim);
  itemp = 0;
  for(int j=0; j < nndim; j++){
    binaryX[j].resize(ncov[j]);
    for(int k=0; k < ncov[j]; k++){
      binaryX[j][k] = binaryX1[itemp];
      itemp++;
    }
  }
  
  size_t minobsnode = (size_t) minobs;
  
  double maxy; /* The reference level */
  maxy = R_NegInf;
  for(size_t i=0; i<di.n_samp;i++){
    if(maxy<y[i]) maxy = y[i];
  }
  
  std::vector<std::vector<double> > allfit; /* The sum of fit of all trees for train.data (nlatent x nsub)*/
  std::vector<std::vector<double> > wtilde; /* (nlatent x nsub) */
  std::vector<std::vector<double> > stemp0; /* The second term in pseudo outcome for bytree*/
  allfit.resize(di.n_dim);
  wtilde.resize(di.n_dim);
  stemp0.resize(di.n_dim);
  for(size_t k=0; k<di.n_dim; k++){
    allfit[k].resize(di.n_samp);
    wtilde[k].resize(di.n_samp);
    stemp0[k].resize(di.n_samp);
  }
  
  for(size_t k=0; k<di.n_dim; k++){
    for(size_t i=0;i<di.n_samp;i++) {
      allfit[k][i] = 0.0;
      wtilde[k][i] = 0.0;
      stemp0[k][i] = 0.0;
    }
  }
  
  std::vector<double> stemp; /* The second term in pseudo outcome for bylatent*/
  stemp.resize(nn);
  
  /* fit of current tree */
  std::vector<std::vector<std::vector<double> > > ftemp; /*m x nsub, (i,j) is the current fit of tree i person j*/
  ftemp.resize(nndim);
  
  for(int i=0; i < nndim; i++){
    ftemp[i].resize(m[i]);
  }
  
  for(int k=0; k < nndim; k++){
    for(int i=0; i < m[k]; i++){
      ftemp[k][i].resize(nn);
      for(int j=0; j < nn; j++){
        ftemp[k][i][j] = 0;
      }
    }  
  }
  
  
  std::vector<double> mvnsample,mvnmean; /* temp storage for mvn samples (nlatent vector) */
  mvnsample.resize(di.n_dim);
  mvnmean.resize(di.n_dim);
  
  //initialize PerSigInv for bylatent
  std::vector<std::vector<double> > PerSigInv; /*Inverse of permuted Sigma*/
  PerSigInv.resize(nndim);
  
  for(int i=0; i < nndim; i++){
    PerSigInv[i].resize(nndim);
  }
  
  //initialize PerSigInv List for bytree
  std::vector<std::vector<std::vector<double> > > PerSigInvList; /*Inverse of permuted Sigma*/
  PerSigInvList.resize(nndim);
  
  for(int i=0; i < nndim; i++){
    PerSigInvList[i].resize(nndim);
  }
  
  for(int i=0; i < nndim; i++){
    for(int j=0; j < nndim; j++){
      PerSigInvList[i][j].resize(nndim);
    }
  }
  
  std::vector<std::vector<double> > r; /* pseudoresponse (nlatent x nsub) */
  r.resize(di.n_dim);
  for(size_t k=0; k<di.n_dim; k++){
    r[k].resize(di.n_samp);
  }
  
  //initialize vtemp = (w - allfit)[-k] 
  std::vector<double> vtemp;
  vtemp.resize(nndim - 1);
  
  std::vector<double> condsig;
  condsig.resize(di.n_dim);
  
  double cvar;//squared condsig at each latent
  
  pinfo pi;
  pi.pswap = pswap; //prob of swap move, default 0.1
  pi.pbd = pbd; //prob of birth/death move, default 0.5
  pi.pb = pb; //prob of birth given  birth/death, default 0.25
  
  pi.alpha = alpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
  pi.beta = beta; //
  //pi.tau=3.0/(kfac*sqrt((double)m)); // categorical outcome
  pi.sigma= 1; //error standard deviation of the latent W
  
  std::vector<double> tauList;
  tauList.resize(nndim);
  for(int i=0; i < nndim; i++){
    tauList[i] = 3.0/(kfac*sqrt((double)m[i]));
  }
  
  //storage for ouput
  //in sample fit
  
  double max_temp = 0.0;// used to find class membership
  int pclass = 0;// used to find class membership
  
  //initialize tree
  
  std::vector<std::vector<tree> > t; /* ntree x nlatent matrix of trees */
  t.resize(nndim);
  for(int i=0; i < nndim; i++){
    t[i].resize(m[i]);
  }
  for(int i=0; i < nndim; i++){
    for(int k=0;k < m[i]; k++) t[i][k].setm(0.00);
  }
  
  ////////
  std::vector<std::vector<double> > mtemp1;
  std::vector<std::vector<double> > WishSample,WishSampleTilde,WishSampleInv, SigmaTmp, SigmaTmpInv;
  WishSampleInv.resize(di.n_dim);
  WishSample.resize(di.n_dim);
  SigmaTmp.resize(di.n_dim);
  SigmaTmpInv.resize(di.n_dim);
  
  WishSampleTilde.resize(di.n_dim);
  
  WishSampleInv.resize(di.n_dim);
  mtemp1.resize(di.n_dim);
  for(size_t j=0;j<di.n_dim;j++){
    WishSample[j].resize(di.n_dim);
    WishSampleTilde[j].resize(di.n_dim);
    mtemp1[j].resize(di.n_dim);
    WishSampleInv[j].resize(di.n_dim);
    SigmaTmp[j].resize(di.n_dim);
    SigmaTmpInv[j].resize(di.n_dim);
    
  }
  
  /* Initialize Zi(r) in Algorithm 3.2, Page 20 of Jiao & van Dyk 2015 */
  std::vector<std::vector<double> > zir;
  zir.resize(di.n_dim);
  for(size_t j=0;j<di.n_dim;j++){
    zir[j].resize(di.n_samp);
  }
  
  
  std::vector<std::vector<double> > WishMat1, WishMat1Inv, WishSampleTildeInv;
  
  WishMat1.resize(di.n_dim);
  WishSampleTildeInv.resize(di.n_dim);
  for(size_t j=0;j<di.n_dim;j++){
    WishMat1[j].resize(di.n_dim);
    WishSampleTildeInv[j].resize(di.n_dim);
  }
  
  //parameters for correction from Jiao & van Dyk 2015
  size_t check_temp;
  size_t max_class_zir = 0;
  double max_zir;
  int itercnt;
  
  
  double alpha2, alpha2old, ss;
  int sigdrawcounter = 0;
  ////////
  
  //percA[k][i] = percA[k*m + i] = average over all percAtmp for draw i latent k
  std::vector<std::vector<double> > percAtmp;
  percAtmp.resize(di.n_dim);
  //numNodes[k][i] =  numNodes[k*m + i] = average over all numNodestmp for draw i latent k
  std::vector<std::vector<double> > numNodestmp;
  numNodestmp.resize(di.n_dim);
  //numLeaves[k][i] = numLeaves[k*m + i] = average over all numLeavestmp for draw i latent k
  std::vector<std::vector<double> > numLeavestmp;
  numLeavestmp.resize(di.n_dim);
  //treeDepth[k][i] = treeDepth[k*m + i] = average over all treeDepthtmp for draw i latent k
  std::vector<std::vector<double> > treeDepthtmp;
  treeDepthtmp.resize(di.n_dim);
  //nLtD[k][i] = nLtD[k*m + i] = number of leaves at tree depth of the ith tree for the current round of draw of latent k
  std::vector<std::vector<double> > nLtDtmp;
  nLtDtmp.resize(di.n_dim);
  
  
  for(int i=0;i < nndim;i++) {
    percAtmp[i].resize(m[i]);
    numNodestmp[i].resize(m[i]);
    numLeavestmp[i].resize(m[i]);
    treeDepthtmp[i].resize(m[i]);
    nLtDtmp[i].resize(m[i]);
  }
  
  for(int k=0; k < nndim; k++){
    for(int i=0;i < m[k];i++) {
      percAtmp[k][i] = 0.0;
      numNodestmp[k][i] = 1.0;
      numLeavestmp[k][i] = 1.0;
      treeDepthtmp[k][i] = 0.0;
      nLtDtmp[k][i] = 1.0;
    }
  }
  
  //incProptmp[j] = count occurence of xj in splitting rules for the current round of draw
  //incProp[j] = average of incProptmp[j] over all draws
  std::vector<std::vector<double> > incProptmp;
  std::vector<std::vector<double> > incProp;
  
  incProptmp.resize(di.n_dim);
  incProp.resize(di.n_dim);
  
  for(int i=0; i < nndim; i++){
    incProptmp[i].resize(ncov[i]);
    incProp[i].resize(ncov[i]);
  }
  
  for(int k=0; k < nndim; k++){
    for(int i=0;i < ncov[k]; i++) {
      incProptmp[k][i] = 0.0;
      incProp[k][i] = 0.0;
    }
  }
  
  ////////
  int tnd = nd + fitMNP + burn;
  NumericVector sigmasample(nndim*nndim*nd), vec_train(nn*nd), wp(nd);
  NumericMatrix percA(nndim, tnd), numNodes(nndim, tnd), numLeaves(nndim, tnd), treeDepth(nndim, tnd); //incProp(nndim, n_cov)
  
  //TreeMod
  std::vector<std::vector<double> > L1; //action type
  L1.resize(nndim);
  std::vector<std::vector<double> > L2; //node id when action type > 0
  L2.resize(nndim);
  std::vector<std::vector<double> > L3; //GROW info (v, c, mul, mur) when action type == 1
  L3.resize(nndim);
  std::vector<std::vector<double> > L4; //PRUNE info (mu) when action type = =2
  L4.resize(nndim);
  std::vector<std::vector<double> > L5; //SWAP info (lind, rind, vd, cd, vk, ck, nbnc.size()) when action type == 3
  L5.resize(nndim);
  std::vector<std::vector<double> > L6; //SWAP nx bottom nodes when action type == 3
  L6.resize(nndim);
  std::vector<std::vector<double> > L7; //CHANGE info (v, c, nbnv.size()) when action type == 4
  L7.resize(nndim);
  std::vector<std::vector<double> > L8; //CHANGE nx bottom nodes when action type == 4
  L8.resize(nndim);
  
  for(int i = 0; i<nndim; i++){
    L1[i].reserve(m[i]*tnd);
    L2[i].reserve(floor(0.5*m[i]*tnd));
    L3[i].reserve(floor(4*0.5*pb*m[i]*tnd));
    L4[i].reserve(floor(0.5*(pbd - pb)*m[i]*tnd));
    L5[i].reserve(floor(7*0.5*pswap*m[i]*tnd));
    L6[i].reserve(floor(5*0.5*pswap*m[i]*tnd));
    L7[i].reserve(floor(3*0.5*(1 - pbd - pswap)*m[i]*tnd));
    L8[i].reserve(floor(5*0.5*(1 - pbd - pswap)*m[i]*tnd));
  }
  
  std::vector<double> w, mu;
  w.resize(nndim*nn);
  mu.resize(nndim*nn);
  for(int i = 0; i < (nndim*nn); i++){
    //w[i] = 0;
    mu[i] = 0;
  }
  
  //MCMC
  
  time_t tp;
  int time1 = time(&tp);
  
  /* Initialize counters for outputs sigmasample, vec_test and vec_train */
  int countvectrain = 0;
  
  /* If fitMNP==0, initialize w to Mcmc$w0 (matrix of 0 if no input) */
  /* If fitMNP > 0, Mcmc$w0 is result from MNP, fitted in preheat */
  for(size_t k=0; k<di.n_dim; k++){
    for(size_t i=0;i<di.n_samp;i++) {
      w[i*di.n_dim + k] = w1(i,k);
    }
  }
  
  if(fitMNP > 0){
    
    //skip Step 1 (a) of sampling w, use w from MNP
    
    /* Step 1 (b) */
    /* mtemp1 = V x inverse(Sigma) */
    ss=0;
    for(size_t j=0;j<di.n_dim;j++){
      for(size_t k=0;k<di.n_dim;k++) {
        mtemp1[j][k]=0;
      }
    }
    
    for(size_t i=0;i<di.n_dim;i++){
      for(size_t j=0;j<di.n_dim;j++){
        for(size_t k=0;k<di.n_dim;k++){
          mtemp1[j][k]+=V[j*di.n_dim + i]*sigmai[i*di.n_dim + k];
        }
      }
    }
    /* ss = trace(V x inverse(Sigma)) */
    for(size_t j=0;j<di.n_dim;j++) ss+=mtemp1[j][j];
    /* alpha^2 = trace(V x inverse(Sigma)) / rchisq */
    alpha2=ss/(double)R::rchisq((double)nu*di.n_dim);
    
    /* Step 1 (c) */
    for(size_t k=0; k<di.n_dim; k++){
      for(size_t i=0;i<di.n_samp;i++) {
        wtilde[k][i] = sqrt(alpha2) * (w[i*di.n_dim + k]);
      }
    }
    
    /* Step 2 */
    for(int ploop=0;ploop<fitMNP;ploop++){
      
      if(ploop%100==0) Rprintf("\n MNP tree fit iteration: %d of %d \n",ploop, fitMNP);
      
      for(int k=0; k < nndim; k++){
        for(int i=0;i < m[k];i++) {
          percAtmp[k][i] = 0.0;
        }
      }
      
      if(bylatent > 0){
        
        
        for(int k=0; k < nndim; k++){
          
          dinvperm( &sigmai[0], nndim, k, PerSigInv);
          cvar=1/PerSigInv[nndim-1][nndim-1];
          condsig[k]=sqrt(cvar);
          
          for(int i=0; i < nn; i++){
            itemp = 0;
            for(int j=0;j < nndim; j++){
              if(j != k)
                vtemp[itemp++] = w[i*nndim + j]-allfit[j][i];
            }
            stemp[i] = 0;
            for(int j=0; j < (nndim-1); j++)
              stemp[i] += PerSigInv[nndim-1][j]*vtemp[j]*cvar;
          }
          
          for(int ntree = 0 ; ntree < m[k]; ntree++){
            
            for(int i=0;i < nn ;i++) {
              allfit[k][i] -= ftemp[k][ntree][i];
              //get pseudo response
              r[k][i] = w[i*nndim + k] - allfit[k][i] + stemp[i];
              
            }
            
            di.y = &r[k][0];
            //pi.sigma = sqrt(alpha2)*condsig[k]; //sqrt psi_k tilde
            pi.sigma = condsig[k];
            pi.tau=tauList[k];
            
            if(dgn){
              bd1F(XMat[k], (size_t) ncov[k], t[k][ntree], xi[k], di, pi, minobsnode, binaryX[k], &nLtDtmp[k][ntree], &percAtmp[k][ntree], &numNodestmp[k][ntree], &numLeavestmp[k][ntree], &treeDepthtmp[k][ntree], incProptmp[k], L1[k], L2[k], L3[k], L4[k], L5[k], L6[k], L7[k], L8[k]);
            } else {
              bdF(XMat[k], (size_t) ncov[k], t[k][ntree], xi[k], di, pi, minobsnode, binaryX[k], L1[k], L2[k], L3[k], L4[k], L5[k], L6[k], L7[k], L8[k]);
            }
            
            fit(t[k][ntree], XMat[k], nn, ncov[k], xi[k], ftemp[k][ntree]);
            for(size_t i=0;i<di.n_samp;i++) {
              allfit[k][i] += ftemp[k][ntree][i]	;
            }
          }//ntree
          
        }//ndim
        
        
        
      } else { 
        //bytrees
        //for drawing sth tree for kth level in tth round, 
        //T(k,s)^t|w^t, T(j,l)^t, [j=1...p-1, l<s], T(i,h)^(t-1), [i=1...p-1, h>=s] 
        
        for(int k=0; k < nndim; k++){
          dinvperm( &sigmai[0], nndim, k, PerSigInvList[k]);
          cvar=1/PerSigInvList[k][nndim-1][nndim-1];
          condsig[k]=sqrt(cvar);//condsig[k] is sqrt psi_k
        }
        
        for(int ntree = 0 ; ntree < mm; ntree++){
          
          for(int k=0; k < nndim; k++){
            
            if(ntree < m[k]){
              
              for(int i=0;i < nn; i++) {
                itemp = 0;
                for(int j=0;j < nndim; j++){
                  if(j != k)
                    vtemp[itemp++] = w[i*nndim + j]-allfit[j][i];
                }
                stemp0[k][i] = 0;
                for(int j=0; j < (nndim-1); j++){
                  stemp0[k][i] += PerSigInvList[k][nndim-1][j]*vtemp[j]*condsig[k]*condsig[k];//pow() requires cmath.h
                }
              }//i
              
            }//if ntree
          }//k
          
          for(int k=0; k < nndim; k++){//start the tree update for ntree at dim k
            
            if(ntree < m[k]){
              
              for(int i=0; i < nn; i++){          
                allfit[k][i] -= ftemp[k][ntree][i];
                r[k][i] = w[i*di.n_dim + k] - allfit[k][i] + stemp0[k][i];
              }
              
              //update tree
              di.y = &r[k][0];
              //pi.sigma = sqrt(alpha2)*condsig[k]; //sqrt psi_k tilde
              pi.sigma = condsig[k];
              pi.tau=tauList[k];
              
              if(dgn){
                bd1F(XMat[k], (size_t) ncov[k], t[k][ntree], xi[k], di, pi, minobsnode, binaryX[k], &nLtDtmp[k][ntree], &percAtmp[k][ntree], &numNodestmp[k][ntree], &numLeavestmp[k][ntree], &treeDepthtmp[k][ntree], incProptmp[k], L1[k], L2[k], L3[k], L4[k], L5[k], L6[k], L7[k], L8[k]);
              } else {
                bdF(XMat[k], (size_t) ncov[k], t[k][ntree], xi[k], di, pi, minobsnode, binaryX[k], L1[k], L2[k], L3[k], L4[k], L5[k], L6[k], L7[k], L8[k]);
              }
              
              fit(t[k][ntree], XMat[k], nn, ncov[k], xi[k], ftemp[k][ntree]);
              
              for(int i=0;i < nn; i++) {
                allfit[k][i] += ftemp[k][ntree][i]	;
              }
              
            }//if ntree
            
          }//for nndim
          
        }//ntree
        
      }//end bytree
      
      
      if(dgn){
        for(int k=0;k < nndim; k++){
          for(int i=0; i < m[k]; i++) {
            percA(k, ploop) += percAtmp[k][i];
            numNodes(k, ploop) += numNodestmp[k][i];
            numLeaves(k, ploop) += numLeavestmp[k][i];
            treeDepth(k, ploop) += treeDepthtmp[k][i];
          }
          percA(k, ploop) /= m[k];
          numNodes(k, ploop) /= m[k];
          numLeaves(k, ploop) /= m[k];
          treeDepth(k, ploop) /= m[k];
        }
      }
    }//finish fitMNP rounds of tree fittings
    
    /* Step 3 */
    
    /* WishMat1 = V+ sum_i Zi*Zi^T*/
    for(size_t j=0;j<di.n_dim;j++){
      for(size_t k=0;k<di.n_dim;k++){
        WishMat1[j][k]=0;
      }
    }
    
    for(size_t i=0;i<di.n_samp;i++){
      for(size_t j=0;j<di.n_dim;j++){
        for(size_t k=0;k<di.n_dim;k++){
          WishMat1[j][k] += (wtilde[j][i]-sqrt(alpha2)*allfit[j][i])* (wtilde[k][i] - sqrt(alpha2)*allfit[k][i]);
        }
      }
    }
    
    for(size_t j=0;j<di.n_dim;j++){
      for(size_t k=0;k<di.n_dim;k++){
        WishMat1[j][k] += V[k*di.n_dim + j] ;
      }
    }
    
    
    dinv(WishMat1 ,di.n_dim,WishMat1Inv);
    
    /* keep the old alpha2 from step 1 sampling (w,alpha)*/
    alpha2old = alpha2;
    
    if(Jiao){
      //correction from Jiao & van Dyk 2015
      check_temp = 0;
      itercnt = 1;
      
      while(check_temp != di.n_samp && itercnt<= nSigDr){
        
        /* Step 3 (a) */
        // generate inverse Sigma
        rWish(WishSampleTildeInv, WishMat1Inv, (int)(nu+di.n_samp),(int)di.n_dim);
        // invert to get Sigma
        dinv(WishSampleTildeInv ,di.n_dim, WishSampleTilde);
        
        // now check condition 10
        
        /* Step 3 (b) */
        // definition of zir needs new alpha2 based on Sigma, alpha2 = trace(Sigma)/p, y has p+1 levels
        alpha2 = 0;
        for(size_t j=0; j< di.n_dim; j++) alpha2 += (WishSampleTilde[j][j])/double(di.n_dim);
        
        /* Step 3 (c) */
        // difine zir
        for(size_t i=0;i<di.n_samp;i++){
          for(size_t j=0;j<di.n_dim;j++){
            zir[j][i] = wtilde[j][i]- (1 - sqrt(alpha2/alpha2old))*sqrt(alpha2old)*allfit[j][i]; //see pseudocode for explanation
          }
        }
        // condition 10 should exist for every training sample
        check_temp = 0;// count samples that satisfies
        for(size_t i=0;i<di.n_samp;i++){
          
          max_zir = R_NegInf;
          for(size_t j=0;j<di.n_dim;j++){
            if(zir[j][i] > max_zir){
              max_zir = zir[j][i];
              max_class_zir = j+1;
            }
          }
          if(max_zir <= 0){
            max_class_zir = (size_t)maxy;
          }
          if((size_t)y[i] == max_class_zir){
            check_temp++;
          }
        }
        
        itercnt++;
      }//while
      
      if(itercnt>= nSigDr){
        std::cout << "\n iteration on Sigma reached upper limit\n";
      }
      //correction end
    } else {
      /* Step 3 (a) */
      // generate inverse Sigma Tilde
      rWish(WishSampleTildeInv, WishMat1Inv, (int)(nu+di.n_samp),(int)di.n_dim);
      // invert to get Sigma Tilde
      dinv(WishSampleTildeInv ,di.n_dim, WishSampleTilde);
      
      /* Step 3 (b) */
      // definition of zir needs new alpha2 based on Sigma, alpha2 = trace(Sigma)/p, y has p+1 levels
      alpha2 = 0;
      for(size_t j=0; j< di.n_dim; j++) alpha2 += (WishSampleTilde[j][j]);
      alpha2 = alpha2/double(di.n_dim);
      
    }
    
    
    /* Step 3 (e) and (f) */
    for(size_t i=0; i<di.n_samp; i++){
      for(size_t k=0; k < di.n_dim; k++){
        mu[i*di.n_dim + k] = allfit[k][i]*sqrt(alpha2old)/sqrt(alpha2); //divide allfit this to transform
        w[i*di.n_dim +k] = allfit[k][i] + (wtilde[k][i]-sqrt(alpha2old)*allfit[k][i]) /sqrt(alpha2) ;
      }
    }
    
    /* Step 3 (d) */
    for(size_t j=0;j<di.n_dim;j++){
      for(size_t k=0;k<di.n_dim;k++){
        sigmai[j*di.n_dim + k] = WishSampleTildeInv[j][k]*alpha2;
      }
    }
    
    
  }//finish prefit MNP
  
  ////////////////
  
  for(int loop=0;loop<(nd+burn);loop++) { /* Start posterior draws */
    
    for(int k=0; k < nndim; k++){
      for(int i=0;i < m[k];i++) {
        percAtmp[k][i] = 0.0;
      }
    }
    
    if(loop%100==0) Rprintf("\n iteration: %d of %d \n",loop, nd+burn);
    
    /* Step 1 (a) */
    draww(&w[0], &mu[0], &sigmai[0], &nn,&nndim, &y[0]);
    //if(loop==0) draww(&w[0], &mu[0], &sigmai[0], &nn,&nndim, &y[0]); /* regenerate w for the initial draw */
    
    /* Step 1 (b) */
    /* mtemp1 = V x inverse(Sigma) */
    ss=0;
    for(size_t j=0;j<di.n_dim;j++){
      for(size_t k=0;k<di.n_dim;k++) {
        mtemp1[j][k]=0;
      }
    }
    
    for(size_t i=0;i<di.n_dim;i++){
      for(size_t j=0;j<di.n_dim;j++){
        for(size_t k=0;k<di.n_dim;k++){
          mtemp1[j][k]+=V[j*di.n_dim + i]*sigmai[i*di.n_dim + k];
        }
      }
    }
    /* ss = trace(V x inverse(Sigma)) */
    for(size_t j=0;j<di.n_dim;j++) ss+=mtemp1[j][j];
    /* alpha^2 = trace(V x inverse(Sigma)) / rchisq */
    alpha2=ss/(double)R::rchisq((double)nu*di.n_dim);
    
    /* Step 1 (c) */
    for(size_t k=0; k<di.n_dim; k++){
      for(size_t i=0;i<di.n_samp;i++) {
        wtilde[k][i] = sqrt(alpha2) * (w[i*di.n_dim + k]);
      }
    }
    
    //done sampling alpha2, w
    
    /* Step 2 */
    /* See tree sampling theory.doc for explanation */
    
    if(bylatent > 0){
      
      for(int k=0; k < nndim; k++){
        
        dinvperm( &sigmai[0], nndim, k, PerSigInv);
        cvar=1/PerSigInv[nndim-1][nndim-1];
        condsig[k]=sqrt(cvar);
        
        for(int i=0; i < nn; i++){
          itemp = 0;
          for(int j=0;j < nndim; j++){
            if(j != k)
              vtemp[itemp++] = w[i*nndim + j]-allfit[j][i];
          }
          stemp[i] = 0;
          for(int j=0; j < (nndim-1); j++)
            stemp[i] += PerSigInv[nndim-1][j]*vtemp[j]*cvar;
        }
        
        for(int ntree = 0 ; ntree < m[k]; ntree++){
          
          for(int i=0;i < nn ;i++) {
            allfit[k][i] -= ftemp[k][ntree][i];
            //get pseudo response
            r[k][i] = w[i*nndim + k] - allfit[k][i] + stemp[i];
            
          }
          
          di.y = &r[k][0];
          //pi.sigma = sqrt(alpha2)*condsig[k]; //sqrt psi_k tilde
          pi.sigma = condsig[k];
          pi.tau=tauList[k];
          
          if(dgn){
            bd1F(XMat[k], (size_t) ncov[k], t[k][ntree], xi[k], di, pi, minobsnode, binaryX[k], &nLtDtmp[k][ntree], &percAtmp[k][ntree], &numNodestmp[k][ntree], &numLeavestmp[k][ntree], &treeDepthtmp[k][ntree], incProptmp[k], L1[k], L2[k], L3[k], L4[k], L5[k], L6[k], L7[k], L8[k]);
          } else {
            bdF(XMat[k], (size_t) ncov[k], t[k][ntree], xi[k], di, pi, minobsnode, binaryX[k], L1[k], L2[k], L3[k], L4[k], L5[k], L6[k], L7[k], L8[k]);
          }
          
          fit(t[k][ntree], XMat[k], nn, ncov[k], xi[k], ftemp[k][ntree]);
          for(size_t i=0;i<di.n_samp;i++) {
            allfit[k][i] += ftemp[k][ntree][i]	;
          }
        }//ntree
        
      }//ndim
      
    } else { 
      //bytrees
      //for drawing sth tree for kth level in tth round, 
      //T(k,s)^t|w^t, T(j,l)^t, [j=1...p-1, l<s], T(i,h)^(t-1), [i=1...p-1, h>=s] 
      
      for(int k=0; k < nndim; k++){
        dinvperm( &sigmai[0], nndim, k, PerSigInvList[k]);
        cvar=1/PerSigInvList[k][nndim-1][nndim-1];
        condsig[k]=sqrt(cvar);//condsig[k] is sqrt psi_k
      }
      
      for(int ntree = 0 ; ntree < mm; ntree++){
        
        for(int k=0; k < nndim; k++){
          
          if(ntree < m[k]){
            
            for(int i=0;i < nn; i++) {
              itemp = 0;
              for(int j=0;j < nndim; j++){
                if(j != k)
                  vtemp[itemp++] = w[i*nndim + j]-allfit[j][i];
              }
              stemp0[k][i] = 0;
              for(int j=0; j < (nndim-1); j++){
                stemp0[k][i] += PerSigInvList[k][nndim-1][j]*vtemp[j]*condsig[k]*condsig[k];//pow() requires cmath.h
              }
            }//i
            
          }//if ntree
        }//k
        
        for(int k=0; k < nndim; k++){//start the tree update for ntree at dim k
          
          if(ntree < m[k]){
            
            for(int i=0; i < nn; i++){          
              allfit[k][i] -= ftemp[k][ntree][i];
              r[k][i] = w[i*di.n_dim + k] - allfit[k][i] + stemp0[k][i];
            }
            
            //update tree
            di.y = &r[k][0];
            //pi.sigma = sqrt(alpha2)*condsig[k]; //sqrt psi_k tilde
            pi.sigma = condsig[k];
            pi.tau=tauList[k];
            
            if(dgn){
              bd1F(XMat[k], (size_t) ncov[k], t[k][ntree], xi[k], di, pi, minobsnode, binaryX[k], &nLtDtmp[k][ntree], &percAtmp[k][ntree], &numNodestmp[k][ntree], &numLeavestmp[k][ntree], &treeDepthtmp[k][ntree], incProptmp[k], L1[k], L2[k], L3[k], L4[k], L5[k], L6[k], L7[k], L8[k]);
            } else {
              bdF(XMat[k], (size_t) ncov[k], t[k][ntree], xi[k], di, pi, minobsnode, binaryX[k], L1[k], L2[k], L3[k], L4[k], L5[k], L6[k], L7[k], L8[k]);
            }
            
            fit(t[k][ntree], XMat[k], nn, ncov[k], xi[k], ftemp[k][ntree]);
            
            for(int i=0;i < nn; i++) {
              allfit[k][i] += ftemp[k][ntree][i]	;
            }
            
          }//if ntree
          
        }//for nndim
        
      }//ntree
      
    }//end bytree
    
    
    //done sampling (T,M)
    
    /* Step 3 */
    
    /* WishMat1 = V+ sum_i Zi*Zi^T*/
    for(size_t j=0;j<di.n_dim;j++){
      for(size_t k=0;k<di.n_dim;k++){
        WishMat1[j][k]=0;
      }
    }
    
    for(size_t i=0;i<di.n_samp;i++){
      for(size_t j=0;j<di.n_dim;j++){
        for(size_t k=0;k<di.n_dim;k++){
          WishMat1[j][k] += (wtilde[j][i]-sqrt(alpha2)*allfit[j][i])* (wtilde[k][i] - sqrt(alpha2)*allfit[k][i]);
        }
      }
    }
    
    for(size_t j=0;j<di.n_dim;j++){
      for(size_t k=0;k<di.n_dim;k++){
        WishMat1[j][k] += V[k*di.n_dim + j] ;
      }
    }
    
    
    dinv(WishMat1 ,di.n_dim,WishMat1Inv);
    
    /* keep the old alpha2 from step 1 sampling (w,alpha)*/
    alpha2old = alpha2;
    
    if(Jiao){
      //correction from Jiao & van Dyk 2015
      check_temp = 0;
      itercnt = 1;
      
      while(check_temp != di.n_samp && itercnt<= nSigDr){
        
        /* Step 3 (a) */
        // generate inverse Sigma
        rWish(WishSampleTildeInv, WishMat1Inv, (int)(nu+di.n_samp),(int)di.n_dim);
        // invert to get Sigma
        dinv(WishSampleTildeInv ,di.n_dim, WishSampleTilde);
        
        // now check condition 10
        
        /* Step 3 (b) */
        // definition of zir needs new alpha2 based on Sigma, alpha2 = trace(Sigma)/p, y has p+1 levels
        alpha2 = 0;
        for(size_t j=0; j< di.n_dim; j++) alpha2 += (WishSampleTilde[j][j])/double(di.n_dim);
        
        /* Step 3 (c) */
        // difine zir
        for(size_t i=0;i<di.n_samp;i++){
          for(size_t j=0;j<di.n_dim;j++){
            zir[j][i] = wtilde[j][i]- (1 - sqrt(alpha2/alpha2old))* sqrt(alpha2old)*allfit[j][i]; //see pseudocode for explanation
          }
        }
        // condition 10 should exist for every training sample
        check_temp = 0;// count samples that satisfies
        for(size_t i=0;i<di.n_samp;i++){
          
          max_zir = R_NegInf;
          for(size_t j=0;j<di.n_dim;j++){
            if(zir[j][i] > max_zir){
              max_zir = zir[j][i];
              max_class_zir = j+1;
            }
          }
          if(max_zir <= 0){
            max_class_zir = (size_t)maxy;
          }
          if((size_t)y[i] == max_class_zir){
            check_temp++;
          }
        }
        
        itercnt++;
      }//while
      
      if(itercnt>= nSigDr){
        std::cout << "\n iteration on Sigma reached upper limit\n";
      }
      //correction end
    } else {
      /* Step 3 (a) */
      // generate inverse Sigma Tilde
      rWish(WishSampleTildeInv, WishMat1Inv, (int)(nu+di.n_samp),(int)di.n_dim);
      // invert to get Sigma Tilde
      dinv(WishSampleTildeInv ,di.n_dim, WishSampleTilde);
      
      /* Step 3 (b) */
      // definition of zir needs new alpha2 based on Sigma, alpha2 = trace(Sigma)/p, y has p+1 levels
      alpha2 = 0;
      for(size_t j=0; j< di.n_dim; j++) alpha2 += (WishSampleTilde[j][j]);
      alpha2 = alpha2/double(di.n_dim);
      
    }
    
    
    /* Step 3 (e) and (f) */
    for(size_t i=0; i<di.n_samp; i++){
      for(size_t k=0; k < di.n_dim; k++){
        mu[i*di.n_dim + k] = allfit[k][i]*sqrt(alpha2old)/sqrt(alpha2); //divide allfit this to transform
        w[i*di.n_dim +k] = allfit[k][i] + (wtilde[k][i]-sqrt(alpha2old)*allfit[k][i]) /sqrt(alpha2) ;
      }
    }
    
    /* Step 3 (d) */
    for(size_t j=0;j<di.n_dim;j++){
      for(size_t k=0;k<di.n_dim;k++){
        sigmai[j*di.n_dim + k] = WishSampleTildeInv[j][k]*alpha2;
        SigmaTmpInv[j][k] = WishSampleTildeInv[j][k]*alpha2;
        //SigmaTmp[j][k] = WishSampleTilde[j][k]/alpha2;
      }
    }
    
    if(dgn){
      for(int k=0;k < nndim; k++){
        for(int i=0; i < m[k]; i++) {
          percA(k, fitMNP+loop) += percAtmp[k][i];
          numNodes(k, fitMNP+loop) += numNodestmp[k][i];
          numLeaves(k, fitMNP+loop) += numLeavestmp[k][i];
          treeDepth(k, fitMNP+loop) += treeDepthtmp[k][i];
        }
        percA(k, fitMNP+loop) /= m[k];
        numNodes(k, fitMNP+loop) /= m[k];
        numLeaves(k, fitMNP+loop) /= m[k];
        treeDepth(k, fitMNP+loop) /= m[k];
      }
    }
    
    if(loop>=burn){
      
      wp[loop-burn] = sqrt(alpha2);
      
      if(dgn){
        for(size_t k=0;k<di.n_dim;k++){
          
          double numsplits = 0;
          for(size_t i=0; i< (size_t)ncov[k]; i++){
            //std::cout<< "loop"<<loop<<" "<<incProptmp[i] <<" \n";
            numsplits += incProptmp[k][i];
          }
          for(size_t i=0; i< (size_t)ncov[k]; i++){
            incProp[k][i] += incProptmp[k][i]/numsplits;//proportion of all splitting rules that uses xi
          }
        }//k
      }//dgn
      
      dinv(SigmaTmpInv ,di.n_dim,SigmaTmp);
      for(size_t j=0;j<di.n_dim;j++){
        for(size_t k=0;k<di.n_dim;k++){
          sigmasample[sigdrawcounter] = SigmaTmp[j][k];
          sigdrawcounter++;
        }
      }
      
      
      for(size_t k = 0; k <di.n_samp; k++){
        max_temp = R_NegInf;
        for(size_t l=0; l<di.n_dim; l++){
          mvnmean[l] = mu[k*di.n_dim + l];
        }
        
        rMVN(mvnsample, mvnmean, SigmaTmp,di.n_dim);
        
        for(int l = 0 ; l < nndim; l++){
          if(mvnsample[l] > max_temp){
            max_temp = mvnsample[l];
            pclass = l+1;
          }
        }
        if(max_temp <=0) {
          pclass = (int)maxy;
        }
        vec_train[countvectrain] = pclass;
        countvectrain++;
        //cout << "pclass: " << pclass << endl;
      }//end prediction for train
      
    }//end prediction for current loop
    
  } //end of loop
  
  if(dgn){
    for(size_t k=0; k<di.n_dim; k++){
      for(size_t i=0; i< (size_t)ncov[k]; i++){
        incProp[k][i] /= nd;
      }
    }
  }
  
  int time2 = time(&tp);
  Rprintf("time for mcmc loop %d secs", time2-time1);
  
  for(int i=0; i < nndim; i++){
    for(int k=0;k < m[i];k++) t[i][k].tonull();
  }//delete trees
  
  List TreeMod = List::create(L1,L2,L3,L4,L5,L6,L7,L8);
  
  List z;
  
  if(dgn){
    z = List::create( Rcpp::Named("samp_y") = vec_train, 
                      Rcpp::Named("sigmasample") = sigmasample,
                      Rcpp::Named("Percent_Acceptance") = percA, 
                      Rcpp::Named("Tree_Num_Nodes") = numNodes, 
                      Rcpp::Named("Tree_Num_Leaves") = numLeaves, 
                      Rcpp::Named("Tree_Depth") = treeDepth, 
                      Rcpp::Named("Inclusion_Proportions") = Rcpp::wrap(incProp),
                      Rcpp::Named("TreeMod") = TreeMod,
                      Rcpp::Named("xi") = Rcpp::wrap(xi),
                      Rcpp::Named("working_param") = wp,
                      Rcpp::Named("maxy") = maxy) ;
  } else {
    z = List::create( Rcpp::Named("samp_y") = vec_train,
                      Rcpp::Named("sigmasample") = sigmasample,
                      Rcpp::Named("TreeMod") = TreeMod,
                      Rcpp::Named("xi") = Rcpp::wrap(xi),
                      Rcpp::Named("working_param") = wp,
                      Rcpp::Named("maxy") = maxy) ;
  }
  
  
  return z ;
  
}
