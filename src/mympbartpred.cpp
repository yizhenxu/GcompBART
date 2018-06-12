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
List mympbartpred(int ndim,
                  int testn,
                  int nd,
                  NumericVector means,
                  NumericVector sigmasample,
                  NumericVector wp,
                  int maxy) {
  
  std::vector<double> mvnsample,mvnmean; /* temp storage for mvn samples (nlatent vector) */
  mvnsample.resize(ndim);
  mvnmean.resize(ndim);
  
  double max_temp = 0.0;// used to find class membership
  int pclass = 0;// used to find class membership
  int mtemp = 0;
  int stemp = 0;
  int countvectest= 0;
  
  std::vector<std::vector<double> > SigmaTmp;
  SigmaTmp.resize(ndim);
  for(int j=0;j<ndim;j++){
    SigmaTmp[j].resize(ndim);
  }
  
  std::vector<double> samp_y; /* The sum of fit of all trees for test.data (nsub_test x 1)*/
  samp_y.resize(nd*testn);
  
  for(int loop = 0; loop < nd; loop++){
    
    for(int j=0;j<ndim;j++){
      for(int k=0;k<ndim;k++){
        SigmaTmp[j][k] = sigmasample[stemp];
        stemp++;
      }
    }
    
    for(int k = 0; k < testn; k++){
      max_temp = R_NegInf;
      
      for(int l=0; l<ndim; l++){
        mvnmean[l] = means[mtemp]/wp[loop] ;
        mtemp++;
      }
      
      rMVN(mvnsample, mvnmean, SigmaTmp, ndim);
      
      for(int l = 0 ; l < ndim; l++){
        if(mvnsample[l] > max_temp){
          max_temp = mvnsample[l];
          pclass = l+1;
        }
      }
      if(max_temp <=0) {
        pclass = maxy;
      }
      samp_y[countvectest] = pclass;
      countvectest++;
      
    }//k
    
  }//loop
  
  
  
  List z;
  z = List::create( Rcpp::Named("samp_y") = samp_y) ;
  
  
  return z ;
  
}
