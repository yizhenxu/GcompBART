#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include "tree.h"
#include "info.h"
#include "bd.h"
#include "funs.h"

/*function to compute the inverse of permuted Sigma */
void dinvperm(double *sigi,
          int	p, int j,
          std::vector<std::vector<double> >& PerSigInv){
  int k, l;
  std::vector<std::vector<double> > SigmaInv, Sigma, PerSig;
  SigmaInv.resize(p);
  Sigma.resize(p);
  PerSig.resize(p);
  
  for(k=0;k<p;k++){
    SigmaInv[k].resize(p);
    Sigma[k].resize(p);
    PerSig[k].resize(p);
  }
  
  int itemp = 0;
  for (k = 0; k < p; k++){
    for (l = 0; l < p; l++){
      SigmaInv[l][k] = sigi[itemp++];
    }
  }
  dinv(SigmaInv, p, Sigma);
  
  /** permutation of Sigma **/
  int kindx, lindx;          /* Indexes for the permutation of Sigma */
  kindx = 0;
  for(k=0;k<p;k++){
    lindx=0;
    for(l=0;l<p;l++){
      if(j!=k)
        if(j!=l)
          PerSig[k+kindx][l+lindx]=Sigma[k][l];
        else{
          lindx=-1;
          PerSig[k+kindx][p-1]=Sigma[k][l];
        }
        else{
          kindx=-1;
          if(j==l){
            lindx=-1;
            PerSig[p-1][p-1]=Sigma[j][j];
          }
          else
            PerSig[p-1][l+lindx]=Sigma[k][l];
        }
    }
  }
  dinv(PerSig, p, PerSigInv);
  
}
  
/* function to compute moments of x[j] | x[-j], j starts from 0, adapted from MNP */
void condmom_mnp(double *x, double *mu, double *sigi, int p, int j, double *m, double *csig)
{
  int k, l;
  std::vector<std::vector<double> > SigmaInv, Sigma, PerSig, PerSigInv;
  SigmaInv.resize(p);
  Sigma.resize(p);
  PerSigInv.resize(p);
  PerSig.resize(p);

  for(k=0;k<p;k++){
    SigmaInv[k].resize(p);
    Sigma[k].resize(p);
    PerSigInv[k].resize(p);
    PerSig[k].resize(p);
  }

  int itemp = 0;
  for (k = 0; k < p; k++){
    for (l = 0; l < p; l++){
      SigmaInv[l][k] = sigi[itemp++];
    }
  }
  dinv(SigmaInv, p, Sigma);

  /** permutation of Sigma **/
  int kindx, lindx;          /* Indexes for the permutation of Sigma */
  kindx = 0;
  for(k=0;k<p;k++){
    lindx=0;
    for(l=0;l<p;l++){
      if(j!=k)
        if(j!=l)
          PerSig[k+kindx][l+lindx]=Sigma[k][l];
        else{
          lindx=-1;
          PerSig[k+kindx][p-1]=Sigma[k][l];
        }
        else{
          kindx=-1;
          if(j==l){
            lindx=-1;
            PerSig[p-1][p-1]=Sigma[j][j];
          }
          else
            PerSig[p-1][l+lindx]=Sigma[k][l];
        }
    }
  }
  dinv(PerSig, p, PerSigInv);

  double cvar;
  cvar=1/PerSigInv[p-1][p-1];
  *csig=sqrt(cvar);

  *m  = mu[j];
  double* vtemp = new double[p-1];
  itemp = 0;
  for(k=0;k<p;k++){
    if(k != j)
      vtemp[itemp++] = x[k]-mu[k];
  }
  for(k=0;k<(p-1);k++)
    *m -= PerSigInv[p-1][k]*vtemp[k]*cvar;

  delete[] vtemp;
  vtemp = 0;
}

double rtrun(double mu, double sigma,double trunpt, int above)
{
  double FA,FB,rnd,result,arg ;
  if (above) {
    FA=0.0; FB=R::pnorm(((trunpt-mu)/(sigma)),0.0,1.0,1,0);
  }
  else {
    FB=1.0; FA=R::pnorm(((trunpt-mu)/(sigma)),0.0,1.0,1,0);
  }

  
  rnd=unif_rand();
  arg=rnd*(FB-FA)+FA;
  if(arg > .999999999) arg=.999999999;
  if(arg < .0000000001) arg=.0000000001;
  result = mu + sigma * R::qnorm(arg,0.0,1.0,1,0);
  
  return result;
}
void drawwi(double *w, double *mu, double *sigmai,int *p, int *y)
{
  /*	function to draw w_i by Gibbing's thru p vector   */

  int i,j,above;
  double bound;
  double mean, csig;

  for (i=0; i < *p; ++i)
  {
    bound=0.0;
    for (j=0; j < *p ; ++j)
      { if (j != i) {bound=R::fmax2(bound,w[j]); }}
    if (*y == i+1)
      above = 0;
    else {
      above = 1;
      if(*y == (*p+1)) bound = 0; //when reference level
    }

    condmom_mnp(w,mu,sigmai,*p, i,&mean,&csig);
    w[i]=rtrun(mean,csig,bound,above);

  }
}

void draww(double *w, double *mu, double *sigmai, int *n, int *p, int *y)
{
  /*	function to gibbs down entire w vector for all n obs  */
  int i, ind;
  for (i=0; i < *n; ++i)
  {
    ind= *p * i;
    drawwi(w+ind,mu+ind,sigmai,p,y+i);
  }
}


//---------------------------------------------------------------
// get cut-points
void getcutpoints(int nc, int n_cov, int n_samp,
                  std::vector<std::vector<double> >& X, xinfo& xi){
  double xinc; //increments
  double xx;
  
  
  std::vector<double> minx(n_cov,R_PosInf); // to store the minimum of each of the individual specific pred
  std::vector<double> maxx(n_cov,R_NegInf);// to store the max of each of the individual specific pred
  
  for(int j=0;j<n_cov;j++) {
    for(int i=0;i<n_samp;i++) {
      xx = X[i][j];
      if(xx < minx[j]) minx[j]=xx;
      if(xx > maxx[j]) maxx[j]=xx;
    }
  }
  
  
  
  //make grid of nc cutpoints between min and max for each x.
  xi.resize(n_cov);
  for(int i=0;i<n_cov;i++) {
    xinc = (maxx[i]-minx[i])/(nc+1.0);
    xi[i].resize(nc);
    for(int j=0;j<nc;j++) xi[i][j] = minx[i] + (j+1)*xinc;
  }
  
  
}

void getcutpoints(int nc, int n_cov, int n_samp,
                  Rcpp::NumericMatrix& X, xinfo& xi){
  double xinc; //increments
  double xx;


  std::vector<double> minx(n_cov,R_PosInf); // to store the minimum of each of the individual specific pred
  std::vector<double> maxx(n_cov,R_NegInf);// to store the max of each of the individual specific pred

  for(int j=0;j<n_cov;j++) {
    for(int i=0;i<n_samp;i++) {
      xx = X(i,j);
      if(xx < minx[j]) minx[j]=xx;
      if(xx > maxx[j]) maxx[j]=xx;
    }
  }



  //make grid of nc cutpoints between min and max for each x.
  xi.resize(n_cov);
  for(int i=0;i<n_cov;i++) {
    xinc = (maxx[i]-minx[i])/(nc+1.0);
    xi[i].resize(nc);
    for(int j=0;j<nc;j++) xi[i][j] = minx[i] + (j+1)*xinc;
  }


}

//--------------------------------------------------
//does this bottom node n have any variables it can split on.
bool cansplit(tree::tree_p n, xinfo& xi)
{
  int L,U;
  bool v_found = false; //have you found a variable you can split on
  size_t v=0;
  while(!v_found && (v < xi.size())) { //invar: splitvar not found, vars left
    L=0; U = xi[v].size()-1;
    n->rg(v,&L,&U);
    if(U>=L) v_found=true;
    v++;
  }
  return v_found;
}
//std X, n_samp, n_cov
void fit(tree& t, std::vector<std::vector<double> >& X, int n_samp, int n_cov, xinfo& xi, std::vector<double>& fv)
{
  double* xx = new double[n_cov];
  tree::tree_cp bn;
  
  for(int i=0;i < n_samp;i++) {
    for(int j=0;j < n_cov; j++){
      xx[j] = X[i][j];
    }
    
    bn = t.bn(xx,xi);
    fv[i] = bn->getm();
  }
  delete[] xx;
  xx = 0;
}
//Numeric X, di
void fit(tree& t, Rcpp::NumericMatrix& X, dinfo di, xinfo& xi, std::vector<double>& fv)
{
  double* xx = new double[di.n_cov];
  tree::tree_cp bn;

  for(size_t i=0;i<di.n_samp;i++) {
    for(size_t j=0;j<di.n_cov; j++){
      xx[j] = X(i,j);
    }

    bn = t.bn(xx,xi);
    fv[i] = bn->getm();
  }
  delete[] xx;
  xx = 0;
}


// get pseudo response


void getpseudoresponse(dinfo& di, std::vector<std::vector<double> >& ftemp,
                       std::vector<std::vector<double> >& rtemp, double *sigmai,
                       std::vector<std::vector<double> >& r, std::vector<double>& condsig){
  double mean, csig;
  double* tempres = new double[di.n_dim];
  double* tempmean = new double[di.n_dim];
  int itemp = 0;
  for(size_t i=0; i<di.n_samp; i++){
    itemp = 0;
    //prediction from current tree for current i
    for(size_t k=0; k<di.n_dim; k++){
      tempmean[itemp++] = ftemp[k][i];
    }
    itemp = 0;
    //simulated latent for current i
    for(size_t k=0; k<di.n_dim; k++){
      tempres[itemp++] = rtemp[k][i];
    }

    for(size_t k=0; k<di.n_dim; k++){
      condmom_mnp(tempres,tempmean,sigmai,(int)di.n_dim,(int)k,&mean,&csig);
      r[k][i] = tempres[k]- mean + tempmean[k];
      if(i==0) condsig[k] = csig;
    }

  }

  delete[] tempres;
  tempres = 0;
  delete[] tempmean;
  tempmean = 0;

}



//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree& t, xinfo& xi, pinfo& pi, tree::npv& goodbots)
{
  double pb;  //prob of birth to be returned
  tree::npv bnv; //all the bottom nodes
  t.getbots(bnv);
  for(size_t i=0;i!=bnv.size();i++)
    if(cansplit(bnv[i],xi)) goodbots.push_back(bnv[i]);
    if(goodbots.size()==0) { //are there any bottom nodes you can split on?
      pb=0.0;
    } else {
      if(t.treesize()==1) pb=1.0; //is there just one node?
      else pb=pi.pb;
    }
    return pb;
}
//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars, Rcpp::NumericVector& binaryX)
{
  int L,U;
  for(size_t v=0;v!=xi.size();v++) {//try each variable
    L=0; U = xi[v].size()-1;
    n->rg(v,&L,&U);
    if(binaryX[v]==0 && U>=L) goodvars.push_back(v);
    if(binaryX[v]==1 && L==0 && U==(int)(xi[v].size()-1) ) goodvars.push_back(v);//binary variable not used in the current path to the root
  }
}
//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else alpha/(1+d)^beta
double pgrow(tree::tree_p n, xinfo& xi, pinfo& pi)
{
  if(cansplit(n,xi)) {
    return pi.alpha/pow(1.0+n->depth(),pi.beta);
  } else {
    return 0.0;
  }
}
//--------------------------------------------------
//get sufficients stats for all bottom nodes
void allsuff(std::vector<std::vector<double> >& X,
             tree& x, xinfo& xi, dinfo& di, tree::npv& bnv,
             std::vector<sinfo>& sv)
{
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations
  size_t ni;         //the  index into vector of the current bottom node
  double* xx = new double[di.n_cov];
  double y;          //current y

  bnv.clear();
  x.getbots(bnv);

  typedef tree::npv::size_type bvsz;
  bvsz nb = bnv.size();
  sv.resize(nb);

  std::map<tree::tree_cp,size_t> bnmap;
  for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;

  for(size_t i=0;i<di.n_samp;i++) {
    for(size_t j=0;j<di.n_cov; j++){
      xx[j] = X[i][j];
    }
    y=di.y[i];

    tbn = x.bn(xx,xi);
    ni = bnmap[tbn];

    ++(sv[ni].n);
    sv[ni].sy += y;
    sv[ni].sy2 += y*y;
  }
  delete[] xx;
  xx = 0;
}
//--------------------------------------------------
//get sufficient stats for children (v,c) of node nx in tree x
//std X, di, ncov
void getsuff(std::vector<std::vector<double> >& X, size_t ncov,
             tree& x, tree::tree_cp nx, size_t v, size_t c,
             xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
  double* xx = new double[ncov];
  double y;  //current y
  sl.n=0;sl.sy=0.0;sl.sy2=0.0;
  sr.n=0;sr.sy=0.0;sr.sy2=0.0;
  
  for(size_t i=0;i<di.n_samp;i++) {
    for(size_t j=0;j<ncov; j++){
      xx[j] = X[i][j];
    }
    
    if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node
      y = di.y[i];
      if(xx[v] < xi[v][c]) {
        sl.n++;
        sl.sy += y;
        sl.sy2 += y*y;
      } else {
        sr.n++;
        sr.sy += y;
        sr.sy2 += y*y;
      }
    }
  }
  // free memory
  delete []xx;
  xx = 0;
  
}
//Numeric X, di
void getsuff(Rcpp::NumericMatrix& X,
             tree& x, tree::tree_cp nx, size_t v, size_t c,
             xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
  double* xx = new double[di.n_cov];
  double y;  //current y
  sl.n=0;sl.sy=0.0;sl.sy2=0.0;
  sr.n=0;sr.sy=0.0;sr.sy2=0.0;

  for(size_t i=0;i<di.n_samp;i++) {
    for(size_t j=0;j<di.n_cov; j++){
      xx[j] = X(i,j);
    }

    if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node
      y = di.y[i];
      if(xx[v] < xi[v][c]) {
        sl.n++;
        sl.sy += y;
        sl.sy2 += y*y;
      } else {
        sr.n++;
        sr.sy += y;
        sr.sy2 += y*y;
      }
    }
  }
  // free memory
  delete []xx;
  xx = 0;

}
//--------------------------------------------------
//get sufficient stats for pair of bottom children nl(left) and nr(right) in tree x
//std X, di, ncov
void getsuff(std::vector<std::vector<double> >& X, size_t ncov,
             tree& x, tree::tree_cp nl, tree::tree_cp nr,
             xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
  
  
  double* xx = new double[ncov];
  double y;  //current y
  sl.n=0;sl.sy=0.0;sl.sy2=0.0;
  sr.n=0;sr.sy=0.0;sr.sy2=0.0;
  
  for(size_t i=0;i<di.n_samp;i++) {
    for(size_t j=0;j<ncov; j++){
      xx[j] = X[i][j];
    }
    tree::tree_cp bn = x.bn(xx,xi);
    if(bn==nl) {
      y = di.y[i];
      sl.n++;
      sl.sy += y;
      sl.sy2 += y*y;
    }
    if(bn==nr) {
      y = di.y[i];
      sr.n++;
      sr.sy += y;
      sr.sy2 += y*y;
    }
  }
  // free memory
  delete []xx;
  xx = 0;
  
}

//Numeric X, di
void getsuff(Rcpp::NumericMatrix& X,
             tree& x, tree::tree_cp nl, tree::tree_cp nr,
             xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{


  double* xx = new double[di.n_cov];
  double y;  //current y
  sl.n=0;sl.sy=0.0;sl.sy2=0.0;
  sr.n=0;sr.sy=0.0;sr.sy2=0.0;

  for(size_t i=0;i<di.n_samp;i++) {
    for(size_t j=0;j<di.n_cov; j++){
      xx[j] = X(i,j);
    }
    tree::tree_cp bn = x.bn(xx,xi);
    if(bn==nl) {
      y = di.y[i];
      sl.n++;
      sl.sy += y;
      sl.sy2 += y*y;
    }
    if(bn==nr) {
      y = di.y[i];
      sr.n++;
      sr.sy += y;
      sr.sy2 += y*y;
    }
  }
  // free memory
  delete []xx;
  xx = 0;

}
//--------------------------------------------------
//log of the integrated likelihood
double lil(size_t n, double sy, double sy2, double sigma, double tau)
{
  double yb,yb2,S,sig2,d;
  double sum, rv;

  yb = sy/n;
  yb2 = yb*yb;
  S = sy2 - (n*yb2);
  sig2 = sigma*sigma;
  d = n*tau*tau + sig2;
  sum = S/sig2 + (n*yb2)/d;
  rv = -(n*LTPI/2.0) - (n-1)*log(sigma) -log(d)/2.0;
  rv = rv -sum/2.0;
  return rv;
}


void drmu(std::vector<std::vector<double> >& X, tree& t, xinfo& xi, dinfo& di, pinfo& pi)
{
  

  tree::npv bnv;
  std::vector<sinfo> sv;
  allsuff(X, t,xi,di,bnv,sv);

  double a = 1.0/(pi.tau * pi.tau);
  double sig2 = pi.sigma * pi.sigma;
  double b,ybar;

  for(tree::npv::size_type i=0;i!=bnv.size();i++) {
    b = sv[i].n/sig2;
    ybar = sv[i].sy/sv[i].n;
    bnv[i]->setm(b*ybar/(a+b) + norm_rand()/sqrt(a+b));
  }


  

}

/* Inverting a matrix (MNP package) */
void dinv(std::vector<std::vector<double> >& X,
          int	size,
          std::vector<std::vector<double> >& X_inv)
{
  int i,j, k, errorM;
  double* pdInv = new double[(int)(size * size)];
  X_inv.resize(size);
  for (j = 0; j < size; j++) X_inv[j].resize(size);


  for (i = 0, j = 0; j < size; j++)
    for (k = 0; k <= j; k++)
      pdInv[i++] = X[k][j];
  F77_CALL(dpptrf)("U", &size, pdInv, &errorM);
  if (!errorM) {
    F77_CALL(dpptri)("U", &size, pdInv, &errorM);
    if (errorM) {
      Rprintf("LAPACK dpptri failed, %d\n", errorM);
      Rcpp::stop("Exiting from dinv().\n");

    }
  }
  else {
    Rprintf("LAPACK dpptrf failed, %d\n", errorM);
    Rcpp::stop("Exiting from dinv().\n");
  }
  for (i = 0, j = 0; j < size; j++) {
    for (k = 0; k <= j; k++) {
      X_inv[j][k] = pdInv[i];
      X_inv[k][j] = pdInv[i++];
    }
  }

  // free memory
  delete [] pdInv;
  pdInv = 0;

}

/* Cholesky decomposition, returns lower triangular matrix (MNP package) */
void dcholdc(std::vector<std::vector<double> >& X, int size, std::vector<std::vector<double> >& L)
{
  int i, j, k, errorM;
  double* pdTemp = new double[(int)(size * size)];
  L.resize(size);

  for (j = 0; j < size; j++) L[j].resize(size);
  for (j = 0, i = 0; j < size; j++)
    for (k = 0; k <= j; k++)
      pdTemp[i++] = X[k][j];
  F77_CALL(dpptrf)("U", &size, pdTemp, &errorM);
  if (errorM) {
    Rprintf("LAPACK dpptrf failed, %d\n", errorM);
    Rcpp::stop("Exiting from dcholdc().\n");
  }
  for (j = 0, i = 0; j < size; j++) {
    for (k = 0; k < size; k++) {
      if(j<k)
        L[j][k] = 0.0;
      else
        L[j][k] = pdTemp[i++];
    }
  }

  // free memory
  delete [] pdTemp;
  pdTemp = 0;

}



void rWish(std::vector<std::vector<double> >& Sample,        /* The matrix with to hold the sample */
std::vector<std::vector<double> >& S,             /* The parameter */
int df,                 /* the degrees of freedom */
int size)               /* The dimension */
{
  

  int i,j,k;

  double* V = new double[(int)size];
  std::vector<std::vector<double> > B, C, N, mtemp;
  B.resize(size); C.resize(size); N.resize(size); mtemp.resize(size);
  for (j = 0; j < size; j++){
    B[j].resize(size); C[j].resize(size); N[j].resize(size); mtemp[j].resize(size);
  }

  for(i=0;i<size;i++) {
    V[i]=R::rchisq((double) df-i-1);
    B[i][i]=V[i];
    for(j=(i+1);j<size;j++)
      N[i][j]=norm_rand();
  }

  for(i=0;i<size;i++) {
    for(j=i;j<size;j++) {
      Sample[i][j]=0;
      Sample[j][i]=0;
      mtemp[i][j]=0;
      mtemp[j][i]=0;
      if(i==j) {
        if(i>0)
          for(k=0;k<j;k++)
            B[j][j]+=N[k][j]*N[k][j];
      }
      else {
        B[i][j]=N[i][j]*sqrt(V[i]);
        if(i>0)
          for(k=0;k<i;k++)
            B[i][j]+=N[k][i]*N[k][j];
      }
      B[j][i]=B[i][j];
    }
  }

  dcholdc(S, size, C);
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
      for(k=0;k<size;k++)
        mtemp[i][j]+=C[i][k]*B[k][j];
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
      for(k=0;k<size;k++)
        Sample[i][j]+=mtemp[i][k]*C[j][k];
  
  // free memory
  delete[] V;
  V = 0;
}

void DrawSigma(dinfo& di, double *V, std::vector<std::vector<double> >& allfit,
               double *w, std::vector<std::vector<double> >& WishSample, int nu)
{

  std::vector<std::vector<double> > WishMat1;
  std::vector<std::vector<double> > epsilon;
  std::vector<std::vector<double> > WishMat1Inv;

  epsilon.resize(di.n_dim);
  WishMat1.resize(di.n_dim);
  for(size_t j=0;j<di.n_dim;j++){
    WishMat1[j].resize(di.n_dim);
    epsilon[j].resize(di.n_samp);
  }

  for(size_t i=0; i<di.n_samp; i++){
    for(size_t k=0; k<di.n_dim; k++){
      epsilon[k][i] = w[i*di.n_dim + k] - allfit[k][i];
    }
  }

  for(size_t j=0;j<di.n_dim;j++){
    for(size_t k=0;k<di.n_dim;k++){
      WishMat1[j][k] = V[j*di.n_dim + k];
    }
  }


  for(size_t i=0; i<di.n_samp; i++){
    for(size_t j=0;j<di.n_dim;j++){
      for(size_t k=0;k<di.n_dim;k++){
        WishMat1[j][k] +=epsilon[j][i]*epsilon[k][i];
      }
    }
  }

  dinv(WishMat1 ,di.n_dim,WishMat1Inv);
  rWish(WishSample, WishMat1Inv, (int)(nu+di.n_samp),(int)di.n_dim);

}


//read X

void readx(std::vector<std::vector<std::vector<double> > >& XMat,dinfo& di, double *pX){
  XMat.resize(di.n_dim);
  for(size_t j=0; j < di.n_dim; j++){
    XMat[j].resize(di.n_samp);
  }

  for(size_t j=0; j < di.n_dim; j++){
    for(size_t i=0; i < di.n_samp; i++){
      XMat[j][i].resize(di.n_cov);
    }
  }


  int itemp = 0;
  for(size_t i=0; i < di.n_samp; i++){
    for(size_t j=0; j <  di.n_dim; j++){
      for(size_t k=0; k< di.n_cov; k++){
        XMat[j][i][k] = pX[itemp++];
      }
    }
  }
}



/*  The Sweep operator  (MNP package) */
void SWP(
    std::vector<std::vector<double> >& X,             // The Matrix to work on
    size_t k,                  //The row to sweep
    size_t size)               // The dim. of X
{

  if (X[k][k] < 10e-20)
    Rcpp::stop("SWP: singular matrix.\n");
  else
    X[k][k]=-1/X[k][k];
  for(size_t i=0;i<size;i++)
    if(i!=k){
      X[i][k]=-X[i][k]*X[k][k];
      X[k][i]=X[i][k];
    }
    for(size_t i=0;i<size;i++)
      for(size_t j=0;j<size;j++)
        if(i!=k && j!=k)
          X[i][j]=X[i][j]+X[i][k]*X[k][j]/X[k][k];

}

// draw from MVN -- adapted from R package MNP
void rMVN(
    std::vector<double>& Sample,
    std::vector<double>& mean,
    std::vector<std::vector<double> >& Var,
    size_t size)
{
  

  std::vector<std::vector<double> > Model;
  Model.resize(size +1);
  for(size_t j=0; j<= size; j++){
    Model[j].resize(size + 1);
  }

  double cond_mean;

  /* draw from mult. normal using SWP */
  for(size_t j=1;j<=size;j++){
    for(size_t k=1;k<=size;k++) {
      Model[j][k]=Var[j-1][k-1];
    }
    Model[0][j]=mean[j-1];
    Model[j][0]=mean[j-1];
  }
  Model[0][0]=-1;
  Sample[0]=(double)norm_rand()*sqrt(Model[1][1])+Model[0][1];
  for(size_t j=2;j<=size;j++){
    SWP(Model,j-1,size+1);
    cond_mean=Model[j][0];
    for(size_t k=1;k<j;k++) cond_mean+=Sample[k-1]*Model[j][k];
    Sample[j-1]=(double)norm_rand()*sqrt(Model[j][j])+cond_mean;
  }

  
}

/*Get indicators lind and rind for the to-swap kid(s)*/
void swapkidInd(tree::tree_cp n1, tree::tree_cp n2, int* lind, int* rind){
  
  double u;
  if( (n1->l) && (n2->l) ){ //both nodes have rules
    if( (n1->v == n2->v) && (n1->c == n2->c) ){ //both nodes have the same rule
      *lind = 1; *rind = 1;
    } else { //both nodes have rules but the rules are different
      u=unif_rand();
      if(u<.5) { //pick left
        *lind = 1; *rind = 0;
      } else { //pick right
        *lind = 0; *rind = 1;
      }
    }
  } else { //only one node has rule
      if(n1->l){ //left child has rule
        *lind = 1; *rind = 0;
      }
      if(n2->l){ //right child has rule
        *lind = 0; *rind = 1;
      }
  }
  
}


  /*Return boolean indicator of swap feasibility*/
  /*n points to the to-swap node (not going to swap at n here)*/
  /*xstar is a to-swap copy of the entire tree*/
  /*nstar points to the to-swap node of xstar*/
  /*xstar is the swapped copy from validswap regardless of the feasibility*/
  bool validswap(tree& xstar,  tree::tree_p nstar, xinfo& xi, int lind, int rind, Rcpp::NumericVector& binaryX){

    //swap rules for nstar
    size_t kidv, kidc; //kid’s rule
    if(lind + rind == 2){
      (nstar->r)->v = nstar->v;  (nstar->r)->c = nstar->c;
      nstar->v = (nstar->l)->v; nstar->c = (nstar->l)->c;
      (nstar->l)->v = (nstar->r)->v;  (nstar->l)->c = (nstar->r)->c;
      kidv = (nstar->r)->v;//splitting var of the child after swapping
    } else {
      if(lind == 1){
        kidv = (nstar->l)->v; kidc = (nstar->l)->c;
        (nstar->l)->v = nstar->v; (nstar->l)->c = nstar->c;
        nstar->v = kidv; nstar->c = kidc;
        kidv = (nstar->l)->v;
      }
      if(rind == 1){
        kidv = (nstar->r)->v; kidc = (nstar->r)->c;
        (nstar->r)->v = nstar->v; (nstar->r)->c = nstar->c;
        nstar->v = kidv; nstar->c = kidc;
        kidv = (nstar->r)->v;
      }
    } //end swap in nstar

    //check range of the splitting variables involved in swapping
    //along subtree paths
    bool valid = false;
    int L, U;
    size_t dadv = nstar->v;//splitting var of the parent after swapping

    if(binaryX[dadv] == 0){
      L = 0; U = xi[dadv].size() -1;
      nstar->rg(dadv, &L, &U); //[L,U] = range of cutoff for dadv at nstar from root
      valid = pathcheck(nstar, dadv, L, U);//is dadv valid in subtree after swap?
    }
    if(binaryX[dadv] == 1){
      //binary dadv (at nstar) should not be used in its descendents (will cause empty nodes)
      if(nstar->nuse(dadv) == 1) valid = true;
    }
    if( (dadv != kidv) && valid){
      //no need to check binary kidv -- it being a previous dad makes it valid in a smaller subtree
      if(binaryX[kidv] == 0){
        L = 0; U = xi[kidv].size() -1;
        nstar->rg(kidv, &L, &U); //[L,U] = range of cutoff for kidv at nstar
        valid = pathcheck(nstar, kidv, L, U); //is kidv valid in subtree after swap?
      }
    }
    return(valid);
  }

  /*Does cutoffs of v in the subtree rooted at n stay in correct range?*/
  /*[L, U] is the range of cutoff for v at the root n*/
  bool pathcheck(tree::tree_p n, size_t v, int L, int U){
    if(!(n->l)){ //if bottom node
      return true;
    } else { //if not bottom node
      size_t rv = n->v;
      if( rv == v){
        int cut;
        cut = (int)(n->c);
        if( (cut>=L) && (cut<=U) ){ //valid, move on to check children
          if( !pathcheck(n->l, v, L, cut-1) ){
            return false;
          } else if ( !pathcheck(n->r, v, cut+1, U) ){
            return false;
          } else {
            return true;
          }
        } else { //invalid, stop on path and return
          return false;
        }
      } else { // rv != v, move on to check children
        if( !pathcheck(n->l, v, L, U) ){
          return false;
        } else if ( !pathcheck(n->r, v, L, U) ){
          return false;
        } else {
          return true;
        }
      }//end rv != v
    }//end if not bottom node
  }//end function

  /*Calculate log prior log(p(T)) of a subtree rooted at n*/
  double lpiT(tree::tree_p n, xinfo& xi, pinfo& pi, Rcpp::NumericVector& binaryX){

    double res;
    std::vector<size_t> goodvars; //variables n can split on

    size_t dn = n->depth();
    double pgrow =  pi.alpha/pow(1.0 + dn,pi.beta); //prior prob of n being non terminal

    if(!(n->l)){ //bottom node
      res = log(1.0-pgrow);
    } else {
      res = log(pgrow);
      getgoodvars(n, xi, goodvars, binaryX);
      res -= log(goodvars.size());

      size_t v = (n->v);
      if(binaryX[v] == 0){
        int L,U;
        L=0; U = xi[v].size()-1;
        n->rg(v,&L,&U);
        res -= log(U-L+1); //U-L+1 is number of available split points
      }
      //when v is binary, p(choose cutoff) = 1 b.c. all cutoffs are the same
      //log(1) = 0

      res += lpiT(n->l,xi, pi, binaryX) + lpiT(n->r,xi, pi, binaryX);
    }

    return(res);
  }

  /* Calculate integrated likelihood of  the subtree */
  /*If upd == 1, update mu in bottom nodes of subtree rooted at n*/
  /*x is the whole tree, n is an internal node of x*/
  double lilT(std::vector<std::vector<double> >& X, tree::tree_p x,
              tree::tree_p n, xinfo& xi, pinfo& pi, dinfo& di, int upd){
    tree::npv bnv; //all the bottom nodes of subtree at n
    n->getbots(bnv);

    double* xx = new double[di.n_cov];
    double y;  //current y
    double yb, yb2, S, d, sig2, tau2, sum;
    double LT = 0;
    double a, b;

    double* bnn = new double[(int) (bnv.size()) ]();
    double* bnsy = new double[(int) (bnv.size()) ]();
    double* bnsy2 = new double[(int) (bnv.size()) ]();

    for(size_t i=0;i<di.n_samp;i++) {

      for(size_t j=0;j<di.n_cov; j++){
        xx[j] = X[i][j];
      }//j

      for(size_t k = 0; k != bnv.size(); k++){//loop over all bottom nodes of n
        if(bnv[k] == x->bn(xx,xi)) { //does the bottom node = xx's bottom node
          y = di.y[i];
          bnn[k]++;
          bnsy[k] += y;
          bnsy2[k] +=  y*y;
          break; //end loop over bnv early if found bn
        }
      }//k
    }//i

    for(size_t k = 0; k != bnv.size(); k++){
      if(bnn[k] == 0.0) return(10086.0);
    }//empty bottom node after action, invalid

    sig2 = pi.sigma * pi.sigma; // sigma^2
    tau2 = pi.tau * pi.tau;
    a= 1.0/tau2; //a = 1/tau^2

    
    for(size_t k = 0; k != bnv.size(); k++){//loop over all bottom nodes of n
      yb = bnsy[k]/ bnn[k];
      yb2 = yb*yb;
      S = bnsy2[k] - (bnn[k]*yb2);
      d = bnn[k]*tau2 + sig2;
      sum = S/sig2 + (bnn[k]*yb2)/d;
      LT += -(bnn[k]*LTPI/2.0) - (bnn[k]-1)*log(pi.sigma) - log(d)/2.0 - sum/2.0;

      //update the mu’s if upd==1
      if(upd == 1){
        b = bnn[k]/sig2; // b=n/sigma^2
        bnv[k]->mu = b*yb/(a+b) + norm_rand()/sqrt(a+b);
      }//upd

    }//k
    

    delete[] xx; xx = 0;
    delete[] bnn; bnn = 0;
    delete[] bnsy; bnsy = 0;
    delete[] bnsy2; bnsy2 = 0;
    return(LT);
  }


  /* Calculate integrated likelihood of the subtree */
  /*If upd == 1, update mu in bottom nodes of subtree rooted at n*/
  /*x is the whole tree, n is an internal node of x*/
  /*this is only different from lilT in finding y in bn using map as in allsuff, might be faster */
  double lilT1(std::vector<std::vector<double> >& X, size_t ncov, tree& x,
               tree::tree_p n, xinfo& xi, pinfo& pi, dinfo& di, tree::npv& nbnv, int upd){
    
    tree::tree_cp tbn; //the pointer to the bottom node for the current observations
    size_t ni;         //the  index into vector of the current bottom node
    
    tree::npv bnv; //all the bottom nodes of x
    x.getbots(bnv);
    
    typedef tree::npv::size_type bvsz;
    bvsz nb = bnv.size();
    
    std::map<tree::tree_cp,size_t> bnmap;
    for(size_t i=0;i!=(size_t) nb;i++) bnmap[bnv[i]]= i;
    
    
    double* xx = new double[ncov];
    double y;  //current y
    double yb, yb2, S, d, sig2, tau2, sum;
    double LT = 0;
    double a, b;
    
    double* bnn = new double[(int) nb]();
    double* bnsy = new double[(int) nb]();
    double* bnsy2 = new double[(int) nb]();
    
    for(size_t i=0;i<di.n_samp;i++) {
      
      for(size_t j=0;j<ncov; j++){
        xx[j] = X[i][j];
      }//j
      
      y=di.y[i];
      
      tbn = x.bn(xx,xi);
      ni = bnmap[tbn];
      
      ++(bnn[ni]);
      bnsy[ni] += y;
      bnsy2[ni] +=  y*y;
      
    }//i
    
    for(size_t k = 0; k != (size_t) nb; k++){
      if(bnn[k] == 0.0) return(10086.0);
    }//empty bottom node after action, invalid
    
    
    //update bn of the subtree rooted at n
    
    nbnv.clear();//all the bottom nodes of subtree n
    n->getbots(nbnv);
    
    sig2 = pi.sigma * pi.sigma; // sigma^2
    tau2 = pi.tau * pi.tau;
    a= 1.0/tau2; //a = 1/tau^2
    
    
    for(bvsz k = 0; k != nbnv.size(); k++){//loop over all bottom nodes of x
      tbn = nbnv[k];// pointer to the kth subtree bn
      ni = bnmap[tbn];//corresponding index in bnv, bnn, bnsy, bnsy2
      
      yb = bnsy[ni]/ bnn[ni];
      yb2 = yb*yb;
      S = bnsy2[ni] - (bnn[ni]*yb2);
      d = bnn[ni]*tau2 + sig2;
      sum = S/sig2 + (bnn[ni]*yb2)/d;
      
      LT += -(bnn[ni]*LTPI/2.0) - (bnn[ni]-1)*log(pi.sigma) - log(d)/2.0 - sum/2.0;
      
      //update the mu’s if upd==1
      if(upd == 1){
        b = bnn[ni]/sig2; // b=n/sigma^2
        bnv[ni]->mu = b*yb/(a+b) + norm_rand()/sqrt(a+b);
      }//upd
      
    }//k
    
    
    delete[] xx; xx = 0;
    delete[] bnn; bnn = 0;
    delete[] bnsy; bnsy = 0;
    delete[] bnsy2; bnsy2 = 0;
    return(LT);
  }

  
  double lilT1(Rcpp::NumericMatrix& X, tree& x,
              tree::tree_p n, xinfo& xi, pinfo& pi, dinfo& di, tree::npv& nbnv, int upd){

    tree::tree_cp tbn; //the pointer to the bottom node for the current observations
    size_t ni;         //the  index into vector of the current bottom node

    tree::npv bnv; //all the bottom nodes of x
    x.getbots(bnv);

    typedef tree::npv::size_type bvsz;
    bvsz nb = bnv.size();

    std::map<tree::tree_cp,size_t> bnmap;
    for(size_t i=0;i!=(size_t) nb;i++) bnmap[bnv[i]]= i;


    double* xx = new double[di.n_cov];
    double y;  //current y
    double yb, yb2, S, d, sig2, tau2, sum;
    double LT = 0;
    double a, b;

    double* bnn = new double[(int) nb]();
    double* bnsy = new double[(int) nb]();
    double* bnsy2 = new double[(int) nb]();

    for(size_t i=0;i<di.n_samp;i++) {

      for(size_t j=0;j<di.n_cov; j++){
        xx[j] = X(i,j);
      }//j

      y=di.y[i];

      tbn = x.bn(xx,xi);
      ni = bnmap[tbn];

      ++(bnn[ni]);
      bnsy[ni] += y;
      bnsy2[ni] +=  y*y;

    }//i

    for(size_t k = 0; k != (size_t) nb; k++){
      if(bnn[k] == 0.0) return(10086.0);
    }//empty bottom node after action, invalid


    //update bn of the subtree rooted at n

    nbnv.clear();//all the bottom nodes of subtree n
    n->getbots(nbnv);

    sig2 = pi.sigma * pi.sigma; // sigma^2
    tau2 = pi.tau * pi.tau;
    a= 1.0/tau2; //a = 1/tau^2

    
    for(bvsz k = 0; k != nbnv.size(); k++){//loop over all bottom nodes of x
      tbn = nbnv[k];// pointer to the kth subtree bn
      ni = bnmap[tbn];//corresponding index in bnv, bnn, bnsy, bnsy2

      yb = bnsy[ni]/ bnn[ni];
      yb2 = yb*yb;
      S = bnsy2[ni] - (bnn[ni]*yb2);
      d = bnn[ni]*tau2 + sig2;
      sum = S/sig2 + (bnn[ni]*yb2)/d;

      LT += -(bnn[ni]*LTPI/2.0) - (bnn[ni]-1)*log(pi.sigma) - log(d)/2.0 - sum/2.0;

      //update the mu’s if upd==1
      if(upd == 1){
        b = bnn[ni]/sig2; // b=n/sigma^2
        bnv[ni]->mu = b*yb/(a+b) + norm_rand()/sqrt(a+b);
      }//upd

    }//k
    

    delete[] xx; xx = 0;
    delete[] bnn; bnn = 0;
    delete[] bnsy; bnsy = 0;
    delete[] bnsy2; bnsy2 = 0;
    return(LT);
  }
