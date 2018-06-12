#include <iostream>
#include "tree.h"
#include "info.h"

#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

using std::cout;
using std::endl;


#define LTPI 1.83787706640934536

void condmom_mnp(double *x, double *mu, double *sigi, int p, int j, double *m, double *csig);
double rtrun(double mu, double sigma,double trunpt, int above);
void drawwi(double *w, double *mu, double *sigmai,int *p, int *y);
void draww(double *w, double *mu, double *sigmai, int *n, int *p, int *y);
// get cut-points
void getcutpoints(int nc, int n_cov, int n_samp,
                  Rcpp::NumericMatrix& X, xinfo& xi);


bool cansplit(tree::tree_p n, xinfo& xi);
//
void fit(tree& t, Rcpp::NumericMatrix& X, dinfo di, xinfo& xi, std::vector<double>& fv);

// get pseudo response
void getpseudoresponse(dinfo& di, std::vector<std::vector<double> >& ftemp,
                       std::vector<std::vector<double> >& rtemp, double *sigmai,
                       std::vector<std::vector<double> >& r, std::vector<double>& condsig);
//--------------------------------------------------
//log of the integrated likelihood
double lil(size_t n, double sy, double sy2, double sigma, double tau);
//--------------------------------------------------
//get sufficient stats for pair of bottom children nl(left) and nr(right) in tree x
void getsuff(Rcpp::NumericMatrix& X,
             tree& x, tree::tree_cp nl, tree::tree_cp nr,
             xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr);

//--------------------------------------------------
//get sufficient stats for children (v,c) of node nx in tree x
void getsuff(Rcpp::NumericMatrix& X,
             tree& x, tree::tree_cp nx, size_t v, size_t c,
             xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr);


//--------------------------------------------------
//get sufficients stats for all bottom nodes
void allsuff(std::vector<std::vector<double> >& X,
             tree& x, xinfo& xi, dinfo& di, tree::npv& bnv,
             std::vector<sinfo>& sv);
//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else alpha/(1+d)^beta
double pgrow(tree::tree_p n, xinfo& xi, pinfo& pi);
//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars, Rcpp::NumericVector& binaryX);
//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree& t, xinfo& xi, pinfo& pi, tree::npv& goodbots);

void drmu(std::vector<std::vector<double> >& X, tree& t, xinfo& xi, dinfo& di, pinfo& pi);


void dinv(std::vector<std::vector<double> >& X,
          int	size,
          std::vector<std::vector<double> >& X_inv);
void dcholdc(std::vector<std::vector<double> >& X, int size, std::vector<std::vector<double> >& L);

void rWish(std::vector<std::vector<double> >& Sample,
           std::vector<std::vector<double> >& S,
           int df,
           int size);


void DrawSigma(dinfo& di, double *V, std::vector<std::vector<double> >& allfit,
               double *w, std::vector<std::vector<double> >& WishSample, int nu);

void readx(std::vector<std::vector<std::vector<double> > >& XMat,dinfo& di, double *pX);

void SWP(std::vector<std::vector<double> >& X,
         size_t k,
         size_t size);

void rMVN(
    std::vector<double>& Sample,
    std::vector<double>& mean,
    std::vector<std::vector<double> >& Var,
    size_t size);

void swapkidInd(tree::tree_cp n1, tree::tree_cp n2, int* lind, int* rind);

bool validswap(tree& xstar,  tree::tree_p nstar, xinfo& xi, int lind, int rind, Rcpp::NumericVector& binaryX);

bool pathcheck(tree::tree_p n, size_t v, int L, int U);

double lpiT(tree::tree_p n, xinfo& xi, pinfo& pi, Rcpp::NumericVector& binaryX);

double lilT(std::vector<std::vector<double> >& X, tree::tree_p x,
            tree::tree_p n, xinfo& xi, pinfo& pi, dinfo& di, int upd);

double lilT1(Rcpp::NumericMatrix& X, tree& x,
             tree::tree_p n, xinfo& xi, pinfo& pi, dinfo& di, tree::npv& nbnv, int upd);
