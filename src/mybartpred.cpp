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
List mybartpred(int Gcomp,
                NumericVector L1,
                NumericVector L2,
                NumericVector L3,
                NumericVector L4,
                NumericVector L5,
                NumericVector L6,
                NumericVector L7,
                NumericVector L8,
                NumericMatrix testXMat,
                int n_cov, 
                int testn, int nthin,
                int nd, int npost,
                int burn,
                int numtrees,
                NumericMatrix xi1) {
  
  dinfo dip;
  dip.n_samp = (size_t)testn; dip.n_cov = (size_t)n_cov; dip.n_dim = 0; dip.y=0;
  
  xinfo xi;
  xi.resize(n_cov);
  for(int i=0;i<n_cov;i++) {
    xi[i].resize(xi1.ncol());
    for(int j=0;j<xi1.ncol();j++) xi[i][j] = xi1(i,j);
  }

  std::vector<double> ppredmeanvec; /* The sum of fit of all trees for test.data (nsub_test x 1)*/
  ppredmeanvec.resize(dip.n_samp);
  
  std::vector<double> fpredtemp; /* for ppredmeanvec (nsub_test x 1) */
  //temporary fit vector to compute prediction
  fpredtemp.resize(dip.n_samp);
  
  //initialize tree
  size_t m = (size_t) numtrees;
  std::vector<tree>  t; /* ntree vector of trees */
  t.resize(m);
  for(size_t i=0; i<m; i++){
    t[i].setm(0.00);
  }
  
  int i1 = 0;
  int i2 = 0;
  int i3 = 0;
  int i4 = 0;
  int i5 = 0;
  int i6 = 0;
  int i7 = 0;
  int i8 = 0;
  tree::npv nbnv;
  
  std::vector<double> vec_test; /* The sum of fit of all trees for test.data (nsub_test x 1)*/
  vec_test.resize(nd*testn);
  
  Rcpp::NumericMatrix testXMat1;
  
  //MCMC
  
  time_t tp;
  int time1 = time(&tp);
  
  /* Initialize counters for outputs psigmasample, vec_test */
  int countvectest = 0;
  int thincount = nthin;
  int countdraw  = 0
  
  for(int loop=0;loop<(npost+burn);loop++) { /* Start posterior draws */
    
    if(loop%100==0) Rprintf("\n iteration: %d of %d \n",loop, npost+burn);
    
    /* Step 1 */
    /* See tree sampling theory.doc for explanation */
    for(size_t ntree = 0 ; ntree <m; ntree++){
      
      if(L1[i1] > 0 ){
        if(L1[i1] == 1 ){
          t[ntree].birth((size_t) L2[i2], (size_t) L3[i3], (size_t) L3[i3+1], L3[i3+2], L3[i3+3]);
          i3 += 4;
        }
        if(L1[i1] == 2 ){
          t[ntree].death((size_t) L2[i2], L4[i4]);
          i4++;
        }
        if(L1[i1] == 3 ){
          
          tree::tree_p nx =  t[ntree].getptr((size_t) L2[i2]);
          
          nx->v = (size_t)L5[i5+2]; nx->c = (size_t)L5[i5+3];
          if(L5[i5] == 1){
            (nx->l)->v = (size_t)L5[i5+4];(nx->l)->c = (size_t)L5[i5+5];
          }
          if(L5[i5+1] == 1){
            (nx->r)->v = (size_t)L5[i5+4];(nx->r)->c = (size_t)L5[i5+5];
          }
          
          nbnv.clear();
          nx->getbots(nbnv);//all the bottom nodes of subtree nx
          for(size_t k = 0; k != (size_t)L5[i5+6]; k++){
            nbnv[k]->mu = L6[i6];
            i6++;
          }
          
          i5 += 7;
        }
        if(L1[i1] == 4 ){
          
          tree::tree_p nx =  t[ntree].getptr((size_t) L2[i2]);
          
          nx->v = (size_t)L7[i7]; nx->c = (size_t)L7[i7+1];
          
          nbnv.clear();
          nx->getbots(nbnv);//all the bottom nodes of subtree nx
          for(size_t k = 0; k != (size_t)L7[i7+2]; k++){
            nbnv[k]->mu = L8[i8];
            i8++;
          }
          
          i7 += 3;
        }
      
      i2++;  
        
      }
      
      i1++;
      
    }//ntree
    
    if(loop>=burn){
      
      if(thincount==nthin){
        for(size_t k=0; k<dip.n_samp; k++){
          ppredmeanvec[k] = 0.0;
        }
        
        if(Gcomp == 1){
          testXMat1 = testXMat(Rcpp::Range(countdraw*testn, (countdraw+1)*testn - 1), Rcpp::Range(0, n_cov - 1));
          countdraw++;
          for(size_t j=0;j<m;j++) {
            fit(t[j], testXMat1, dip, xi, fpredtemp);
            for(size_t k=0;k<dip.n_samp;k++) ppredmeanvec[k] += fpredtemp[k];
          }
        } else {
          for(size_t j=0;j<m;j++) {
            fit(t[j], testXMat, dip, xi, fpredtemp);
            for(size_t k=0;k<dip.n_samp;k++) ppredmeanvec[k] += fpredtemp[k];
          }
        }
        
        for(size_t k = 0; k <dip.n_samp; k++){
          vec_test[countvectest] = ppredmeanvec[k];
          countvectest++;
        }//end prediction for test
        
        thincount = 0;
      } else {
        thincount++;
      }//end thincount==nthin
      
      
    }//end loop>=burn
    
    
  } //end of loop
  

  int time2 = time(&tp);
  Rprintf("time for mcmc loop %d secs", time2-time1);
  
  for(size_t i=0; i<m; i++){
    t[i].tonull();
  }//delete trees
  
  List z;
  z = List::create( Rcpp::Named("vec_test") = vec_test) ;
  
  
  return z ;
  
}
