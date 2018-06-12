#ifndef GUARD_bd_h
#define GUARD_bd_h

#include "info.h"
#include "tree.h"

#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>


bool bd(Rcpp::NumericMatrix& X, tree& x, xinfo& xi, dinfo& di, pinfo& pi, size_t minobsnode, Rcpp::NumericVector& binaryX,
        std::vector<double>& L1,std::vector<double>& L2,std::vector<double>& L3,std::vector<double>& L4,
        std::vector<double>& L5,std::vector<double>& L6,std::vector<double>& L7,std::vector<double>& L8);

bool bd1(Rcpp::NumericMatrix& X, tree& x, xinfo& xi, dinfo& di, pinfo& pi, size_t minobsnode, Rcpp::NumericVector& binaryX, double* nLtD, double* pA, double* nN, double* nL, double* tD, std::vector<double>& iP,
         std::vector<double>& L1,std::vector<double>& L2,std::vector<double>& L3,std::vector<double>& L4,
         std::vector<double>& L5,std::vector<double>& L6,std::vector<double>& L7,std::vector<double>& L8);

#endif
