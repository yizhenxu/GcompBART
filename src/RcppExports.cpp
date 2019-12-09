// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// mybartmod
List mybartmod(NumericMatrix XMat, double sigest, int nn, NumericVector y, int n_cov, int nu, double lambda, int nd, int burn, int numtrees, double kfac, double pswap, double pbd, double pb, double alpha, double beta, int nc, int minobs, NumericVector binaryX, int dgn);
RcppExport SEXP _GcompBART_mybartmod(SEXP XMatSEXP, SEXP sigestSEXP, SEXP nnSEXP, SEXP ySEXP, SEXP n_covSEXP, SEXP nuSEXP, SEXP lambdaSEXP, SEXP ndSEXP, SEXP burnSEXP, SEXP numtreesSEXP, SEXP kfacSEXP, SEXP pswapSEXP, SEXP pbdSEXP, SEXP pbSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP ncSEXP, SEXP minobsSEXP, SEXP binaryXSEXP, SEXP dgnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type XMat(XMatSEXP);
    Rcpp::traits::input_parameter< double >::type sigest(sigestSEXP);
    Rcpp::traits::input_parameter< int >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n_cov(n_covSEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type nd(ndSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type numtrees(numtreesSEXP);
    Rcpp::traits::input_parameter< double >::type kfac(kfacSEXP);
    Rcpp::traits::input_parameter< double >::type pswap(pswapSEXP);
    Rcpp::traits::input_parameter< double >::type pbd(pbdSEXP);
    Rcpp::traits::input_parameter< double >::type pb(pbSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type minobs(minobsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type binaryX(binaryXSEXP);
    Rcpp::traits::input_parameter< int >::type dgn(dgnSEXP);
    rcpp_result_gen = Rcpp::wrap(mybartmod(XMat, sigest, nn, y, n_cov, nu, lambda, nd, burn, numtrees, kfac, pswap, pbd, pb, alpha, beta, nc, minobs, binaryX, dgn));
    return rcpp_result_gen;
END_RCPP
}
// mybartpred
List mybartpred(int Gcomp, NumericVector L1, NumericVector L2, NumericVector L3, NumericVector L4, NumericVector L5, NumericVector L6, NumericVector L7, NumericVector L8, NumericMatrix testXMat, int n_cov, int testn, int nthin, int nd, int npost, int burn, int numtrees, NumericMatrix xi1);
RcppExport SEXP _GcompBART_mybartpred(SEXP GcompSEXP, SEXP L1SEXP, SEXP L2SEXP, SEXP L3SEXP, SEXP L4SEXP, SEXP L5SEXP, SEXP L6SEXP, SEXP L7SEXP, SEXP L8SEXP, SEXP testXMatSEXP, SEXP n_covSEXP, SEXP testnSEXP, SEXP nthinSEXP, SEXP ndSEXP, SEXP npostSEXP, SEXP burnSEXP, SEXP numtreesSEXP, SEXP xi1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Gcomp(GcompSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type L1(L1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type L2(L2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type L3(L3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type L4(L4SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type L5(L5SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type L6(L6SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type L7(L7SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type L8(L8SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type testXMat(testXMatSEXP);
    Rcpp::traits::input_parameter< int >::type n_cov(n_covSEXP);
    Rcpp::traits::input_parameter< int >::type testn(testnSEXP);
    Rcpp::traits::input_parameter< int >::type nthin(nthinSEXP);
    Rcpp::traits::input_parameter< int >::type nd(ndSEXP);
    Rcpp::traits::input_parameter< int >::type npost(npostSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type numtrees(numtreesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xi1(xi1SEXP);
    rcpp_result_gen = Rcpp::wrap(mybartpred(Gcomp, L1, L2, L3, L4, L5, L6, L7, L8, testXMat, n_cov, testn, nthin, nd, npost, burn, numtrees, xi1));
    return rcpp_result_gen;
END_RCPP
}
// mympbartmod
List mympbartmod(NumericMatrix XMat, NumericMatrix sigmai1, NumericMatrix V1, int nu, int nn, int nndim, NumericVector y1, int n_cov, int nd, int burn, int numtrees, int nSigDr, double kfac, double pswap, double pbd, double pb, double alpha, double beta, int nc, int minobs, NumericVector binaryX, int dgn, int Jiao, int pMDA);
RcppExport SEXP _GcompBART_mympbartmod(SEXP XMatSEXP, SEXP sigmai1SEXP, SEXP V1SEXP, SEXP nuSEXP, SEXP nnSEXP, SEXP nndimSEXP, SEXP y1SEXP, SEXP n_covSEXP, SEXP ndSEXP, SEXP burnSEXP, SEXP numtreesSEXP, SEXP nSigDrSEXP, SEXP kfacSEXP, SEXP pswapSEXP, SEXP pbdSEXP, SEXP pbSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP ncSEXP, SEXP minobsSEXP, SEXP binaryXSEXP, SEXP dgnSEXP, SEXP JiaoSEXP, SEXP pMDASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type XMat(XMatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sigmai1(sigmai1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type V1(V1SEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< int >::type nndim(nndimSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< int >::type n_cov(n_covSEXP);
    Rcpp::traits::input_parameter< int >::type nd(ndSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type numtrees(numtreesSEXP);
    Rcpp::traits::input_parameter< int >::type nSigDr(nSigDrSEXP);
    Rcpp::traits::input_parameter< double >::type kfac(kfacSEXP);
    Rcpp::traits::input_parameter< double >::type pswap(pswapSEXP);
    Rcpp::traits::input_parameter< double >::type pbd(pbdSEXP);
    Rcpp::traits::input_parameter< double >::type pb(pbSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type minobs(minobsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type binaryX(binaryXSEXP);
    Rcpp::traits::input_parameter< int >::type dgn(dgnSEXP);
    Rcpp::traits::input_parameter< int >::type Jiao(JiaoSEXP);
    Rcpp::traits::input_parameter< int >::type pMDA(pMDASEXP);
    rcpp_result_gen = Rcpp::wrap(mympbartmod(XMat, sigmai1, V1, nu, nn, nndim, y1, n_cov, nd, burn, numtrees, nSigDr, kfac, pswap, pbd, pb, alpha, beta, nc, minobs, binaryX, dgn, Jiao, pMDA));
    return rcpp_result_gen;
END_RCPP
}
// mympbartmod1
List mympbartmod1(NumericMatrix XMat, NumericMatrix sigmai1, NumericMatrix V1, int nu, int nn, int nndim, NumericVector y1, int n_cov, int nd, int burn, int numtrees, int nSigDr, double kfac, double pswap, double pbd, double pb, double alpha, double beta, int nc, int minobs, NumericVector binaryX, int dgn, int Jiao, int pMDA);
RcppExport SEXP _GcompBART_mympbartmod1(SEXP XMatSEXP, SEXP sigmai1SEXP, SEXP V1SEXP, SEXP nuSEXP, SEXP nnSEXP, SEXP nndimSEXP, SEXP y1SEXP, SEXP n_covSEXP, SEXP ndSEXP, SEXP burnSEXP, SEXP numtreesSEXP, SEXP nSigDrSEXP, SEXP kfacSEXP, SEXP pswapSEXP, SEXP pbdSEXP, SEXP pbSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP ncSEXP, SEXP minobsSEXP, SEXP binaryXSEXP, SEXP dgnSEXP, SEXP JiaoSEXP, SEXP pMDASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type XMat(XMatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sigmai1(sigmai1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type V1(V1SEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< int >::type nndim(nndimSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< int >::type n_cov(n_covSEXP);
    Rcpp::traits::input_parameter< int >::type nd(ndSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type numtrees(numtreesSEXP);
    Rcpp::traits::input_parameter< int >::type nSigDr(nSigDrSEXP);
    Rcpp::traits::input_parameter< double >::type kfac(kfacSEXP);
    Rcpp::traits::input_parameter< double >::type pswap(pswapSEXP);
    Rcpp::traits::input_parameter< double >::type pbd(pbdSEXP);
    Rcpp::traits::input_parameter< double >::type pb(pbSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type minobs(minobsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type binaryX(binaryXSEXP);
    Rcpp::traits::input_parameter< int >::type dgn(dgnSEXP);
    Rcpp::traits::input_parameter< int >::type Jiao(JiaoSEXP);
    Rcpp::traits::input_parameter< int >::type pMDA(pMDASEXP);
    rcpp_result_gen = Rcpp::wrap(mympbartmod1(XMat, sigmai1, V1, nu, nn, nndim, y1, n_cov, nd, burn, numtrees, nSigDr, kfac, pswap, pbd, pb, alpha, beta, nc, minobs, binaryX, dgn, Jiao, pMDA));
    return rcpp_result_gen;
END_RCPP
}
// mympbartmod2
List mympbartmod2(NumericMatrix XMat, NumericMatrix sigmai1, NumericMatrix V1, int nu, int nn, int nndim, NumericVector y1, int n_cov, int nd, int burn, int numtrees, int nSigDr, double kfac, double pswap, double pbd, double pb, double alpha, double beta, int nc, int minobs, NumericVector binaryX, int dgn, int Jiao, NumericMatrix w1, int fitMNP);
RcppExport SEXP _GcompBART_mympbartmod2(SEXP XMatSEXP, SEXP sigmai1SEXP, SEXP V1SEXP, SEXP nuSEXP, SEXP nnSEXP, SEXP nndimSEXP, SEXP y1SEXP, SEXP n_covSEXP, SEXP ndSEXP, SEXP burnSEXP, SEXP numtreesSEXP, SEXP nSigDrSEXP, SEXP kfacSEXP, SEXP pswapSEXP, SEXP pbdSEXP, SEXP pbSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP ncSEXP, SEXP minobsSEXP, SEXP binaryXSEXP, SEXP dgnSEXP, SEXP JiaoSEXP, SEXP w1SEXP, SEXP fitMNPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type XMat(XMatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sigmai1(sigmai1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type V1(V1SEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< int >::type nndim(nndimSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< int >::type n_cov(n_covSEXP);
    Rcpp::traits::input_parameter< int >::type nd(ndSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type numtrees(numtreesSEXP);
    Rcpp::traits::input_parameter< int >::type nSigDr(nSigDrSEXP);
    Rcpp::traits::input_parameter< double >::type kfac(kfacSEXP);
    Rcpp::traits::input_parameter< double >::type pswap(pswapSEXP);
    Rcpp::traits::input_parameter< double >::type pbd(pbdSEXP);
    Rcpp::traits::input_parameter< double >::type pb(pbSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type minobs(minobsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type binaryX(binaryXSEXP);
    Rcpp::traits::input_parameter< int >::type dgn(dgnSEXP);
    Rcpp::traits::input_parameter< int >::type Jiao(JiaoSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w1(w1SEXP);
    Rcpp::traits::input_parameter< int >::type fitMNP(fitMNPSEXP);
    rcpp_result_gen = Rcpp::wrap(mympbartmod2(XMat, sigmai1, V1, nu, nn, nndim, y1, n_cov, nd, burn, numtrees, nSigDr, kfac, pswap, pbd, pb, alpha, beta, nc, minobs, binaryX, dgn, Jiao, w1, fitMNP));
    return rcpp_result_gen;
END_RCPP
}
// mympbartmod3
List mympbartmod3(NumericMatrix& XMat1, NumericMatrix& sigmai1, NumericMatrix& V1, int nu, int nn, int nndim, NumericVector& y1, IntegerVector& ncov, int nd, int burn, IntegerVector& m, int nSigDr, double kfac, double pswap, double pbd, double pb, double alpha, double beta, int nc, int minobs, NumericVector& binaryX1, int dgn, int Jiao, NumericMatrix& w1, int fitMNP, int bylatent);
RcppExport SEXP _GcompBART_mympbartmod3(SEXP XMat1SEXP, SEXP sigmai1SEXP, SEXP V1SEXP, SEXP nuSEXP, SEXP nnSEXP, SEXP nndimSEXP, SEXP y1SEXP, SEXP ncovSEXP, SEXP ndSEXP, SEXP burnSEXP, SEXP mSEXP, SEXP nSigDrSEXP, SEXP kfacSEXP, SEXP pswapSEXP, SEXP pbdSEXP, SEXP pbSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP ncSEXP, SEXP minobsSEXP, SEXP binaryX1SEXP, SEXP dgnSEXP, SEXP JiaoSEXP, SEXP w1SEXP, SEXP fitMNPSEXP, SEXP bylatentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type XMat1(XMat1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type sigmai1(sigmai1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type V1(V1SEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< int >::type nndim(nndimSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type ncov(ncovSEXP);
    Rcpp::traits::input_parameter< int >::type nd(ndSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type nSigDr(nSigDrSEXP);
    Rcpp::traits::input_parameter< double >::type kfac(kfacSEXP);
    Rcpp::traits::input_parameter< double >::type pswap(pswapSEXP);
    Rcpp::traits::input_parameter< double >::type pbd(pbdSEXP);
    Rcpp::traits::input_parameter< double >::type pb(pbSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type minobs(minobsSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type binaryX1(binaryX1SEXP);
    Rcpp::traits::input_parameter< int >::type dgn(dgnSEXP);
    Rcpp::traits::input_parameter< int >::type Jiao(JiaoSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type w1(w1SEXP);
    Rcpp::traits::input_parameter< int >::type fitMNP(fitMNPSEXP);
    Rcpp::traits::input_parameter< int >::type bylatent(bylatentSEXP);
    rcpp_result_gen = Rcpp::wrap(mympbartmod3(XMat1, sigmai1, V1, nu, nn, nndim, y1, ncov, nd, burn, m, nSigDr, kfac, pswap, pbd, pb, alpha, beta, nc, minobs, binaryX1, dgn, Jiao, w1, fitMNP, bylatent));
    return rcpp_result_gen;
END_RCPP
}
// mympbartpred
List mympbartpred(int ndim, int testn, int nd, NumericVector means, NumericVector sigmasample, NumericVector wp, int maxy);
RcppExport SEXP _GcompBART_mympbartpred(SEXP ndimSEXP, SEXP testnSEXP, SEXP ndSEXP, SEXP meansSEXP, SEXP sigmasampleSEXP, SEXP wpSEXP, SEXP maxySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ndim(ndimSEXP);
    Rcpp::traits::input_parameter< int >::type testn(testnSEXP);
    Rcpp::traits::input_parameter< int >::type nd(ndSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type means(meansSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmasample(sigmasampleSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wp(wpSEXP);
    Rcpp::traits::input_parameter< int >::type maxy(maxySEXP);
    rcpp_result_gen = Rcpp::wrap(mympbartpred(ndim, testn, nd, means, sigmasample, wp, maxy));
    return rcpp_result_gen;
END_RCPP
}
// mypbartmod
List mypbartmod(NumericMatrix XMat, int nn, NumericVector y, int n_cov, int nd, int burn, int numtrees, double kfac, double pswap, double pbd, double pb, double alpha, double beta, int nc, int minobs, NumericVector binaryX, int dgn);
RcppExport SEXP _GcompBART_mypbartmod(SEXP XMatSEXP, SEXP nnSEXP, SEXP ySEXP, SEXP n_covSEXP, SEXP ndSEXP, SEXP burnSEXP, SEXP numtreesSEXP, SEXP kfacSEXP, SEXP pswapSEXP, SEXP pbdSEXP, SEXP pbSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP ncSEXP, SEXP minobsSEXP, SEXP binaryXSEXP, SEXP dgnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type XMat(XMatSEXP);
    Rcpp::traits::input_parameter< int >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n_cov(n_covSEXP);
    Rcpp::traits::input_parameter< int >::type nd(ndSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type numtrees(numtreesSEXP);
    Rcpp::traits::input_parameter< double >::type kfac(kfacSEXP);
    Rcpp::traits::input_parameter< double >::type pswap(pswapSEXP);
    Rcpp::traits::input_parameter< double >::type pbd(pbdSEXP);
    Rcpp::traits::input_parameter< double >::type pb(pbSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type minobs(minobsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type binaryX(binaryXSEXP);
    Rcpp::traits::input_parameter< int >::type dgn(dgnSEXP);
    rcpp_result_gen = Rcpp::wrap(mypbartmod(XMat, nn, y, n_cov, nd, burn, numtrees, kfac, pswap, pbd, pb, alpha, beta, nc, minobs, binaryX, dgn));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _GcompBART_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GcompBART_mybartmod", (DL_FUNC) &_GcompBART_mybartmod, 20},
    {"_GcompBART_mybartpred", (DL_FUNC) &_GcompBART_mybartpred, 18},
    {"_GcompBART_mympbartmod", (DL_FUNC) &_GcompBART_mympbartmod, 24},
    {"_GcompBART_mympbartmod1", (DL_FUNC) &_GcompBART_mympbartmod1, 24},
    {"_GcompBART_mympbartmod2", (DL_FUNC) &_GcompBART_mympbartmod2, 25},
    {"_GcompBART_mympbartmod3", (DL_FUNC) &_GcompBART_mympbartmod3, 26},
    {"_GcompBART_mympbartpred", (DL_FUNC) &_GcompBART_mympbartpred, 7},
    {"_GcompBART_mypbartmod", (DL_FUNC) &_GcompBART_mypbartmod, 17},
    {"_GcompBART_rcpp_hello_world", (DL_FUNC) &_GcompBART_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_GcompBART(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
