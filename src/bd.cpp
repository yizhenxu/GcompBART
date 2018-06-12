#include <Rcpp.h>
#include <iostream>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include "info.h"
#include "tree.h"
#include "bd.h"
#include "funs.h"

/*
notation: (as in old code): going from state x to state y (eg, incoming tree is x).
note: rather than have x and making a tree y
we just figure out what we need from x, the drawn bottom node,the drawn (v,c).
note sure what the right thing to do is.
Could make y (using a birth) and figure stuff out from y.
That is how the old code works.
*/

bool bd(Rcpp::NumericMatrix& X, tree& x, xinfo& xi, dinfo& di, pinfo& pi, size_t minobsnode, Rcpp::NumericVector& binaryX,
        std::vector<double>& L1,std::vector<double>& L2,std::vector<double>& L3,std::vector<double>& L4,
        std::vector<double>& L5,std::vector<double>& L6,std::vector<double>& L7,std::vector<double>& L8)
{
  //GetRNGstate();
  tree::npv goodbots;  //nodes we could birth at (split on)
  double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x
  double u;
  u = unif_rand();
  //Rprintf("u = %f, PBx = %f\n", u, PBx);
  if(u < PBx) { //do birth

    //--------------------------------------------------
    //draw proposal

    //invalid birth if no goodbots
    if(goodbots.size() == 0){
      L1.push_back(0);
      return false;
    }

    //draw bottom node, choose node index ni from list in goodbots
    size_t ni = floor(unif_rand()*goodbots.size());
    tree::tree_p nx = goodbots[ni]; //pointer to the bottom node we might birth at

    //draw v,  the variable
    std::vector<size_t> goodvars; //variables nx can split on
    getgoodvars(nx,xi,goodvars, binaryX);
    size_t vi = floor(unif_rand()*goodvars.size()); //index of chosen split variable
    size_t v = goodvars[vi];

    //draw c, the cutpoint
    int L,U;
    L=0; U = xi[v].size()-1;
    nx->rg(v,&L,&U);
    size_t c = L + floor(unif_rand()*(U-L+1)); //U-L+1 is number of available split points

    //--------------------------------------------------
    //compute things needed for metropolis ratio

    double Pbotx = 1.0/goodbots.size(); //proposal prob of choosing nx
    size_t dnx = nx->depth();
    double PGnx = pi.alpha/pow(1.0 + dnx,pi.beta); //prior prob of growing at nx

    double PGly, PGry; //prior probs of growing at new children (l and r) of proposal
    if(goodvars.size()>1) { //know there are variables we could split l and r on
      PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta); //depth of new nodes would be one more
      PGry = PGly;
    } else { //only had v to work with, if it is exhausted at either child need PG=0
      if((int)(c-1)<L) { //v exhausted in new left child l, new upper limit would be c-1
        PGly = 0.0;
      } else {
        PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
      }
      if(U < (int)(c+1)) { //v exhausted in new right child r, new lower limit would be c+1
        PGry = 0.0;
      } else {
        PGry = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
      }
    }



    double PDy; //prob of proposing death at y (prune l and r so nx is back to bn)
    if(goodbots.size()>1) { //can birth at y because splittable nodes left
      PDy = 1.0 - pi.pb;
    } else { //nx was the only node you could split on
      if((PGry==0) && (PGly==0)) { //cannot birth at y
        PDy=1.0;
      } else { //y can birth at either l or r
        PDy = 1.0 - pi.pb;
      }
    }

    double Pnogy; //death prob of choosing the nog node at y
    size_t nnogs = x.nnogs();
    tree::tree_cp nxp = nx->getp();//constant pointer to parent node of nx
    if(nxp==0) { //no parent, nx is the top and only node
      Pnogy=1.0;
    } else {
       if(nxp->isnog()) { //if parent is a nog, number of nogs same at x and y
        Pnogy = 1.0/nnogs;
      } else { //if parent is not a nog, y has one more nog.
        Pnogy = 1.0/(nnogs+1.0);
      }
    }

    //--------------------------------------------------
    //compute sufficient statistics
    sinfo sl,sr; //sl for left from nx and sr for right from nx (using rule (v,c))
    getsuff(X, x,nx,v,c,xi,di,sl,sr);
    //--------------------------------------------------
    //compute alpha

    double alpha=0.0,alpha1=0.0,alpha2=0.0;
    double lill=0.0,lilr=0.0,lilt=0.0;
    if((sl.n>=minobsnode) && (sr.n>=minobsnode)) {
      lill = lil(sl.n,sl.sy,sl.sy2,pi.sigma,pi.tau);
      lilr = lil(sr.n,sr.sy,sr.sy2,pi.sigma,pi.tau);
      lilt = lil(sl.n+sr.n,sl.sy+sr.sy,sl.sy2+sr.sy2,pi.sigma,pi.tau);

      alpha1 = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*PBx*Pbotx);
      alpha2 = alpha1*exp(lill+lilr-lilt);
      alpha = std::min(1.0,alpha2);
    } else {
      alpha=0.0;
    }


    //--------------------------------------------------
    //finally ready to try metrop
    double a,b,s2,yb;
    double mul,mur; //means for new bottom nodes, left and right
    if(unif_rand() < alpha) {
      //draw mul, mean for left node
      a= 1.0/(pi.tau*pi.tau); //a = 1/tau^2
      s2 = pi.sigma*pi.sigma; // sigma^2
      //left mean
      yb = sl.sy/sl.n;
      b = sl.n/s2; // b=n/sigma^2
      mul = b*yb/(a+b) + norm_rand()/sqrt(a+b);
      //draw mul, mean for left node
      yb = sr.sy/sr.n;
      b = sr.n/s2; // b=n/sigma^2
      mur = b*yb/(a+b) + norm_rand()/sqrt(a+b);
      x.birth(nx->nid(),v,c,mul,mur);
      
      L1.push_back(1);
      L2.push_back((double)nx->nid());
      L3.push_back((double) v);L3.push_back((double) c);L3.push_back(mul);L3.push_back(mur);
      
      return true;
    } else {
      L1.push_back(0);
      return false;
    }
  } else if(u < pi.pbd) { //do death
    //--------------------------------------------------
    //draw proposal

    //draw nog node, any nog node is a possibility
    tree::npv nognds; //nog nodes
    x.getnogs(nognds);

    //invalid death if no nog nodes
    if(nognds.size() == 0){
      L1.push_back(0);
      return false;
    }

    size_t ni = floor(unif_rand()*nognds.size());
    tree::tree_p nx = nognds[ni]; //the nog node we might kill children at

    //--------------------------------------------------
    //compute things needed for metropolis ratio

    double PGny; //prob the nog node grows
    size_t dny = nx->depth();
    PGny = pi.alpha/pow(1.0+dny,pi.beta);

    //better way to code these two?
    double PGlx = pgrow(nx->getl(),xi,pi);
    double PGrx = pgrow(nx->getr(),xi,pi);

    double PBy;  //prob of birth move at y
    if(!(nx->p)) { //is the nog node nx the top node
      PBy = 1.0;
    } else {
      PBy = pi.pb;
    }

    double Pboty;  //prob of choosing the nog as bot to split on when y
    int ngood = goodbots.size();
    if(cansplit(nx->getl(),xi)) --ngood; //if can split at left child, lose this one
    if(cansplit(nx->getr(),xi)) --ngood; //if can split at right child, lose this one
    ++ngood;  //know you can split at nx
    Pboty=1.0/ngood;

    double PDx = 1.0-PBx; //prob of a death step at x
    double Pnogx = 1.0/nognds.size();

    //--------------------------------------------------
    //compute sufficient statistics
    sinfo sl,sr; //sl for left from nx and sr for right from nx (using rule (v,c))
    getsuff(X, x,nx->getl(),nx->getr(),xi,di,sl,sr);
    //--------------------------------------------------
    //compute alpha

    double lill = lil(sl.n,sl.sy,sl.sy2,pi.sigma,pi.tau);
    double lilr = lil(sr.n,sr.sy,sr.sy2,pi.sigma,pi.tau);
    double lilt = lil(sl.n+sr.n,sl.sy+sr.sy,sl.sy2+sr.sy2,pi.sigma,pi.tau);

    double alpha1 = ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
    double alpha2 = alpha1*exp(lilt - lill - lilr);
    double alpha = std::min(1.0,alpha2);


    //--------------------------------------------------
    //finally ready to try metrop
    double a,b,s2,yb;
    double mu;
    size_t n;
    if(unif_rand()<alpha) {
      //draw mu for nog (which will be bot)
      n = sl.n + sr.n;
      a= 1.0/(pi.tau*pi.tau); //a = 1/tau^2
      s2 = pi.sigma*pi.sigma; // sigma^2
      yb = (sl.sy+sr.sy)/n;
      b = n/s2; // b=n/sigma^2
      mu = b*yb/(a+b) + norm_rand()/sqrt(a+b);
      //do death
      x.death(nx->nid(),mu);

      L1.push_back(2);
      L2.push_back((double)nx->nid());
      L4.push_back(mu);
      
      return true;
    } else {
      L1.push_back(0);
      return false;
    }
  } else if(u < (pi.pbd + pi.pswap)) { //do swap (pswap 0.1)
    //--------------------------------------------------
    //draw proposal

    //draw swappable parent node, any non-nog internal node is a possibility
    tree::npv dadnds; //vector of pointers to swappable parent nodes
    x.getswapdad(dadnds);
    if(dadnds.size()==0){ //no swappable pairs
      L1.push_back(0);
      return false;
    }
    size_t ni = floor(unif_rand()*dadnds.size());
    tree::tree_p nx = dadnds[ni]; //pointer to the parent node we might swap at

    //check which child/children to swap
    int lind = 0, rind = 0;
    swapkidInd(nx->getl(), nx->getr(), &lind, &rind); // update lind and rind

    //make a copy xstar of the x, swap xstar at node nxstar (same location as nx to x)
    //check if the swap is valid at nxstar
    bool checkrule;
    tree xstar;//to-swap copy of x
    tree::tree_p xp =  x.getptr(1); //pointer to x
    xstar.pcp(xp); //get a copy of x (by its pointer xp) to xstar
    tree::tree_p nxstar = xstar.getptr(nx->nid());//pointer to the swap parent node of xstar
    checkrule = validswap(xstar, nxstar, xi, lind, rind, binaryX);
    if(!checkrule){ //not a valid swap
      xstar.tonull();
      nxstar = 0;
      L1.push_back(0);
      return false;
    }

    //--------------------------------------------------
    //compute things needed for metropolis ratio

    double XLogPi, XLogL, YLogPi, YLogL;
    tree::npv nxbnv, nxstarbnv;
    YLogL = lilT1(X, xstar, nxstar, xi, pi, di, nxstarbnv, 1); // update bottom mu’s of the swapped tree copy
    //is there empty bottom node(s)?
    if(YLogL == 10086.0){ //empty bn from swap, invalid
      xstar.tonull();
      nxstar = 0;
      L1.push_back(0);
      return false;
    }
    YLogPi = lpiT(nxstar, xi, pi, binaryX);

    XLogPi = lpiT(nx, xi, pi, binaryX);
    XLogL = lilT1(X, x, nx, xi, pi, di, nxbnv, 0); //do not update bottom mu’s

    double alpha;
    alpha = std::min(1.0,exp(YLogPi+YLogL-XLogPi-XLogL));

    if(unif_rand() < alpha) {
      //update xp (x) with xstar
      nx->v = nxstar->v; nx->c = nxstar->c;
      if(lind == 1){
        (nx->l)->v = (nxstar->l)->v;(nx->l)->c = (nxstar->l)->c;
      }
      if(rind == 1){
        (nx->r)->v = (nxstar->r)->v;(nx->r)->c = (nxstar->r)->c;
      }

      for(size_t k = 0; k != nxbnv.size(); k++){
        nxbnv[k]->mu = nxstarbnv[k]->mu;
      }

      xstar.tonull();
      nxstar = 0;
      
      L1.push_back(3);
      L2.push_back((double)nx->nid());
      L5.push_back((double) lind);L5.push_back((double) rind);
      L5.push_back((double) nx->v);L5.push_back((double) nx->c);
      if(lind == 1){
        L5.push_back((double) (nx->l)->v);L5.push_back((double) (nx->l)->c);
      } else {
        L5.push_back((double) (nx->r)->v);L5.push_back((double) (nx->r)->c);
      }
      L5.push_back((double) nxbnv.size());
      for(size_t k = 0; k != nxbnv.size(); k++){
        L6.push_back(nxbnv[k]->mu);
      }
      
      return true;
    } else {
      xstar.tonull();
      nxstar = 0;
      L1.push_back(0);
      return false;
    }

  } else { //do change (pchange 0.4)
    //--------------------------------------------------
    //draw proposal
    //draw internal node for rule change
    tree::npv intnds; //vector of pointers to internal nodes
    x.getintn(intnds);
    if(intnds.size()==0){ //no changeable internal nodes
      L1.push_back(0);
      return false;
    }
    size_t ni = floor(unif_rand()*intnds.size());
    tree::tree_p nx = intnds[ni]; //pointer to the internal node we want to change rule

    //draw v,  the variable
    std::vector<size_t> goodvars; //variables nx can split on (based on root to nx)
    getgoodvars(nx,xi,goodvars, binaryX);
    size_t vi = floor(unif_rand()*goodvars.size()); //index of chosen split variable
    size_t v = goodvars[vi];

    if(binaryX[v] == 1){
      //if binary v appeared in nx's descendents, this change will cause empty nodes -- invalid
      //when nx->v is also v, it's guaranteed that v only appeared once (at nx) on its path
      if((v!=(nx->v)) && nx->nuse(v) > 0){
        L1.push_back(0);
        return false;
      } 
        
    }

    int L,U;
    L=0; U = xi[v].size()-1;
    if(binaryX[v] == 0){
      nx->rg(v, &L, &U);//range of v at nx -- from root to nx
      nx->downrg(v,&L,&U); //range of v at nx -- from nx to its bottom nodes
    }
    //it's meaningless looking for cutoff range for binary X...

    size_t c;
    if(U>=L){ //randomly pick cutoff (can happen c==nx->c, v==nx->v)
      c = L + floor(unif_rand()*(U-L+1));
    } else {
      L1.push_back(0);
      return false;
    }



    //make a copy xstar of the x, change xstar at node nxstar (same location as nx to x)
    tree xstar;
    tree::tree_p xp =  x.getptr(1); //pointer to x
    xstar.pcp(xp); //get a copy of x (by its pointer xp) to xstar
    tree::tree_p nxstar = xstar.getptr(nx->nid()); //point nstar to the same location in xstar as nx in xp
    nxstar->v = v;
    nxstar->c = c;

    //--------------------------------------------------
    //compute things needed for metropolis ratio

    double XLogPi, XLogL, YLogPi, YLogL;
    tree::npv nxbnv, nxstarbnv;
    YLogL = lilT1(X, xstar, nxstar, xi, pi, di, nxstarbnv, 1); // update bottom mu’s of the changed tree copy
    //is there empty bottom node(s)?
    if(YLogL == 10086.0){ //empty bn from rule change, invalid
      xstar.tonull();
      nxstar = 0;
      L1.push_back(0);
      return false;
    }
    YLogPi = lpiT(nxstar, xi, pi, binaryX);

    XLogPi = lpiT(nx, xi, pi, binaryX);
    XLogL = lilT1(X, x, nx, xi, pi, di, nxbnv, 0); //do not update bottom mu’s

    double alpha;
    alpha = std::min(1.0,exp(YLogPi+YLogL-XLogPi-XLogL));

    if(unif_rand() < alpha) {
      //replace xp (x) with xstar
      nx->v = nxstar->v; nx->c = nxstar->c;

      for(size_t k = 0; k != nxbnv.size(); k++){
        nxbnv[k]->mu = nxstarbnv[k]->mu;
      }

      xstar.tonull();
      nxstar = 0;
      
      L1.push_back(4);
      L2.push_back((double) nx->nid());
      L7.push_back((double) nx->v);L7.push_back((double) nx->c);
      L7.push_back((double) nxbnv.size());
      for(size_t k = 0; k != nxbnv.size(); k++){
        L8.push_back(nxbnv[k]->mu);
      }
      
      return true;
    } else {
      xstar.tonull();
      nxstar = 0;
      L1.push_back(0);
      return false;
    }

  } //End if-else
  //PutRNGstate();
}


//---------------------------------------------------------------


bool bd1(Rcpp::NumericMatrix& X, tree& x, xinfo& xi, dinfo& di, pinfo& pi, size_t minobsnode, Rcpp::NumericVector& binaryX,
        double* nLtD, double* pA, double* nN, double* nL, double* tD, std::vector<double>& iP,
        std::vector<double>& L1,std::vector<double>& L2,std::vector<double>& L3,std::vector<double>& L4,
        std::vector<double>& L5,std::vector<double>& L6,std::vector<double>& L7,std::vector<double>& L8)
{
  //GetRNGstate();
  tree::npv goodbots;  //nodes we could birth at (split on)
  double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x
  double u;
  u = unif_rand();
  //Rprintf("u = %f, PBx = %f\n", u, PBx);
  if(u < PBx) { //do birth

    //--------------------------------------------------
    //draw proposal

    //invalid birth if no goodbots
    if(goodbots.size() == 0){
      (*pA) = 0.0;
      L1.push_back(0);
      return false;
    }

    //draw bottom node, choose node index ni from list in goodbots
    size_t ni = floor(unif_rand()*goodbots.size());
    tree::tree_p nx = goodbots[ni]; //pointer to the bottom node we might birth at

    //draw v,  the variable
    std::vector<size_t> goodvars; //variables nx can split on
    getgoodvars(nx,xi,goodvars, binaryX);
    size_t vi = floor(unif_rand()*goodvars.size()); //index of chosen split variable
    size_t v = goodvars[vi];

    //draw c, the cutpoint
    int L,U;
    L=0; U = xi[v].size()-1;
    nx->rg(v,&L,&U);
    size_t c = L + floor(unif_rand()*(U-L+1)); //U-L+1 is number of available split points

    //--------------------------------------------------
    //compute things needed for metropolis ratio

    double Pbotx = 1.0/goodbots.size(); //proposal prob of choosing nx
    size_t dnx = nx->depth();
    double PGnx = pi.alpha/pow(1.0 + dnx,pi.beta); //prior prob of growing at nx

    double PGly, PGry; //prior probs of growing at new children (l and r) of proposal
    if(goodvars.size()>1) { //know there are variables we could split l and r on
      PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta); //depth of new nodes would be one more
      PGry = PGly;
    } else { //only had v to work with, if it is exhausted at either child need PG=0
      if((int)(c-1)<L) { //v exhausted in new left child l, new upper limit would be c-1
        PGly = 0.0;
      } else {
        PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
      }
      if(U < (int)(c+1)) { //v exhausted in new right child r, new lower limit would be c+1
        PGry = 0.0;
      } else {
        PGry = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
      }
    }



    double PDy; //prob of proposing death at y (prune l and r so nx is back to bn)
    if(goodbots.size()>1) { //can birth at y because splittable nodes left
      PDy = 1.0 - pi.pb;
    } else { //nx was the only node you could split on
      if((PGry==0) && (PGly==0)) { //cannot birth at y
        PDy=1.0;
      } else { //y can birth at either l or r
        PDy = 1.0 - pi.pb;
      }
    }

    double Pnogy; //death prob of choosing the nog node at y
    size_t nnogs = x.nnogs();
    tree::tree_cp nxp = nx->getp();//constant pointer to parent node of nx
    if(nxp==0) { //no parent, nx is the top and only node
      Pnogy=1.0;
    } else {
      if(nxp->isnog()) { //if parent is a nog, number of nogs same at x and y
        Pnogy = 1.0/nnogs;
      } else { //if parent is not a nog, y has one more nog.
        Pnogy = 1.0/(nnogs+1.0);
      }
    }

    //--------------------------------------------------
    //compute sufficient statistics
    sinfo sl,sr; //sl for left from nx and sr for right from nx (using rule (v,c))
    getsuff(X, x,nx,v,c,xi,di,sl,sr);
    //--------------------------------------------------
    //compute alpha

    double alpha=0.0,alpha1=0.0,alpha2=0.0;
    double lill=0.0,lilr=0.0,lilt=0.0;
    if((sl.n>=minobsnode) && (sr.n>=minobsnode)) {
      lill = lil(sl.n,sl.sy,sl.sy2,pi.sigma,pi.tau);
      lilr = lil(sr.n,sr.sy,sr.sy2,pi.sigma,pi.tau);
      lilt = lil(sl.n+sr.n,sl.sy+sr.sy,sl.sy2+sr.sy2,pi.sigma,pi.tau);

      alpha1 = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*PBx*Pbotx);
      alpha2 = alpha1*exp(lill+lilr-lilt);
      alpha = std::min(1.0,alpha2);
    } else {
      alpha=0.0;
    }


    //--------------------------------------------------
    //finally ready to try metrop
    double a,b,s2,yb;
    double mul,mur; //means for new bottom nodes, left and right
    if(unif_rand() < alpha) {
      //draw mul, mean for left node
      a= 1.0/(pi.tau*pi.tau); //a = 1/tau^2
      s2 = pi.sigma*pi.sigma; // sigma^2
      //left mean
      yb = sl.sy/sl.n;
      b = sl.n/s2; // b=n/sigma^2
      mul = b*yb/(a+b) + norm_rand()/sqrt(a+b);
      //draw mul, mean for left node
      yb = sr.sy/sr.n;
      b = sr.n/s2; // b=n/sigma^2
      mur = b*yb/(a+b) + norm_rand()/sqrt(a+b);
      x.birth(nx->nid(),v,c,mul,mur);

      L1.push_back(1);
      L2.push_back((double) nx->nid());
      L3.push_back((double) v);L3.push_back((double) c);L3.push_back(mul);L3.push_back(mur);
      
      (*pA) = 1.0;
      (*nN) += 2.0;
      (*nL)++;
      if((*tD) == (double)dnx){
        (*tD)++;
        (*nLtD) = 2.0;
      } else if ((*tD)-1 == (double)dnx) {
        (*nLtD) += 2.0;
      }
      iP[v]++;

      return true;
    } else {
      (*pA) = 0.0;
      L1.push_back(0);
      return false;
    }
  } else if(u < pi.pbd) { //do death
    //--------------------------------------------------
    //draw proposal

    //draw nog node, any nog node is a possibility
    tree::npv nognds; //nog nodes
    x.getnogs(nognds);

    //invalid death if no nog nodes
    if(nognds.size() == 0){
      (*pA) = 0.0;
      L1.push_back(0);
      return false;
    }

    size_t ni = floor(unif_rand()*nognds.size());
    tree::tree_p nx = nognds[ni]; //the nog node we might kill children at

    //--------------------------------------------------
    //compute things needed for metropolis ratio

    double PGny; //prob the nog node grows
    size_t dny = nx->depth();
    PGny = pi.alpha/pow(1.0+dny,pi.beta);

    //better way to code these two?
    double PGlx = pgrow(nx->getl(),xi,pi);
    double PGrx = pgrow(nx->getr(),xi,pi);

    double PBy;  //prob of birth move at y
    if(!(nx->p)) { //is the nog node nx the top node
      PBy = 1.0;
    } else {
      PBy = pi.pb;
    }

    double Pboty;  //prob of choosing the nog as bot to split on when y
    int ngood = goodbots.size();
    if(cansplit(nx->getl(),xi)) --ngood; //if can split at left child, lose this one
    if(cansplit(nx->getr(),xi)) --ngood; //if can split at right child, lose this one
    ++ngood;  //know you can split at nx
    Pboty=1.0/ngood;

    double PDx = 1.0-PBx; //prob of a death step at x
    double Pnogx = 1.0/nognds.size();

    //--------------------------------------------------
    //compute sufficient statistics
    sinfo sl,sr; //sl for left from nx and sr for right from nx (using rule (v,c))
    getsuff(X, x,nx->getl(),nx->getr(),xi,di,sl,sr);
    //--------------------------------------------------
    //compute alpha

    double lill = lil(sl.n,sl.sy,sl.sy2,pi.sigma,pi.tau);
    double lilr = lil(sr.n,sr.sy,sr.sy2,pi.sigma,pi.tau);
    double lilt = lil(sl.n+sr.n,sl.sy+sr.sy,sl.sy2+sr.sy2,pi.sigma,pi.tau);

    double alpha1 = ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
    double alpha2 = alpha1*exp(lilt - lill - lilr);
    double alpha = std::min(1.0,alpha2);


    //--------------------------------------------------
    //finally ready to try metrop
    double a,b,s2,yb;
    double mu;
    size_t n;
    if(unif_rand()<alpha) {
      //draw mu for nog (which will be bot)
      n = sl.n + sr.n;
      a= 1.0/(pi.tau*pi.tau); //a = 1/tau^2
      s2 = pi.sigma*pi.sigma; // sigma^2
      yb = (sl.sy+sr.sy)/n;
      b = n/s2; // b=n/sigma^2
      mu = b*yb/(a+b) + norm_rand()/sqrt(a+b);
      //do death
      x.death(nx->nid(),mu);

      L1.push_back(2);
      L2.push_back((double) nx->nid());
      L4.push_back(mu);
      
      (*pA) = 1.0;
      (*nN) -= 2.0;
      (*nL)--;
      if((*tD)-1 == (double)dny){
        if((*nLtD) > 2){//depth *tD unchanged
          (*nLtD) -= 2;
        } else {//*nLtD==2
          (*tD)--;
          (*nLtD) = x.Dbots(0, *tD);
        }
      }

      return true;
    } else {
      (*pA) = 0.0;
      L1.push_back(0);
      return false;
    }
  } else if(u < (pi.pbd + pi.pswap)) { //do swap (pswap 0.1)
    //--------------------------------------------------
    //draw proposal

    //draw swappable parent node, any non-nog internal node is a possibility
    tree::npv dadnds; //vector of pointers to swappable parent nodes
    x.getswapdad(dadnds);
    if(dadnds.size()==0){ //no swappable pairs
      (*pA) = 0.0;
      L1.push_back(0);
      return false;
    }
    size_t ni = floor(unif_rand()*dadnds.size());
    tree::tree_p nx = dadnds[ni]; //pointer to the parent node we might swap at

    //check which child/children to swap
    int lind = 0, rind = 0;
    swapkidInd(nx->getl(), nx->getr(), &lind, &rind); // update lind and rind

    //make a copy xstar of the x, swap xstar at node nxstar (same location as nx to x)
    //check if the swap is valid at nxstar
    bool checkrule;
    tree xstar;//to-swap copy of x
    tree::tree_p xp =  x.getptr(1); //pointer to x
    xstar.pcp(xp); //get a copy of x (by its pointer xp) to xstar
    tree::tree_p nxstar = xstar.getptr(nx->nid());//pointer to the swap parent node of xstar
    checkrule = validswap(xstar, nxstar, xi, lind, rind, binaryX);
    if(!checkrule){ //not a valid swap
      xstar.tonull();
      nxstar = 0;
      (*pA) = 0.0;
      L1.push_back(0);
      return false;
    }

    //--------------------------------------------------
    //compute things needed for metropolis ratio

    double XLogPi, XLogL, YLogPi, YLogL;
    tree::npv nxbnv, nxstarbnv;
    YLogL = lilT1(X, xstar, nxstar, xi, pi, di, nxstarbnv, 1); // update bottom mu’s of the swapped tree copy
    //is there empty bottom node(s)?
    if(YLogL == 10086.0){ //empty bn from swap, invalid
      xstar.tonull();
      nxstar = 0;
      (*pA) = 0.0;
      L1.push_back(0);
      return false;
    }
    YLogPi = lpiT(nxstar, xi, pi, binaryX);

    XLogPi = lpiT(nx, xi, pi, binaryX);
    XLogL = lilT1(X, x, nx, xi, pi, di, nxbnv, 0); //do not update bottom mu’s

    double alpha;
    alpha = std::min(1.0,exp(YLogPi+YLogL-XLogPi-XLogL));

    if(unif_rand() < alpha) {
      //update xp (x) with xstar
      nx->v = nxstar->v; nx->c = nxstar->c;
      if(lind == 1){
        (nx->l)->v = (nxstar->l)->v;(nx->l)->c = (nxstar->l)->c;
      }
      if(rind == 1){
        (nx->r)->v = (nxstar->r)->v;(nx->r)->c = (nxstar->r)->c;
      }

      for(size_t k = 0; k != nxbnv.size(); k++){
        nxbnv[k]->mu = nxstarbnv[k]->mu;
      }

      xstar.tonull();
      nxstar = 0;
      
      L1.push_back(3);
      L2.push_back((double) nx->nid());
      L5.push_back((double) lind);L5.push_back((double) rind);
      L5.push_back((double) nx->v);L5.push_back((double) nx->c);
      if(lind == 1){
        L5.push_back((double) (nx->l)->v);L5.push_back((double) (nx->l)->c);
      } else {
        L5.push_back((double) (nx->r)->v);L5.push_back((double) (nx->r)->c);
      }
      L5.push_back(nxbnv.size());
      for(size_t k = 0; k != nxbnv.size(); k++){
        L6.push_back(nxbnv[k]->mu);
      }
      
      (*pA) = 1.0;
      if(lind == 1 && rind == 1){
        iP[nx->v]--;
        iP[(nx->l)->v]++;
      }

      return true;
    } else {
      xstar.tonull();
      nxstar = 0;
      (*pA) = 0.0;
      L1.push_back(0);
      return false;
    }

  } else { //do change (pchange 0.4)
    //--------------------------------------------------
    //draw proposal
    //draw internal node for rule change
    tree::npv intnds; //vector of pointers to internal nodes
    x.getintn(intnds);
    if(intnds.size()==0){ //no changeable internal nodes
      (*pA) = 0.0;
      L1.push_back(0);
      return false;
    }
    size_t ni = floor(unif_rand()*intnds.size());
    tree::tree_p nx = intnds[ni]; //pointer to the internal node we want to change rule

    //draw v,  the variable
    std::vector<size_t> goodvars; //variables nx can split on (based on root to nx)
    getgoodvars(nx,xi,goodvars, binaryX);
    size_t vi = floor(unif_rand()*goodvars.size()); //index of chosen split variable
    size_t v = goodvars[vi];

    if(binaryX[v] == 1){
      //if binary v appeared in nx's descendents, this change will cause empty nodes -- invalid
      //when nx->v is also v, it's guaranteed that v only appeared once (at nx) on its path
      if((v!=(nx->v)) && nx->nuse(v) > 0){
        (*pA) = 0.0;
        L1.push_back(0);
        return false;
      }
    }

    int L,U;
    L=0; U = xi[v].size()-1;
    if(binaryX[v] == 0){
      nx->rg(v, &L, &U);//range of v at nx -- from root to nx
      nx->downrg(v,&L,&U); //range of v at nx -- from nx to its bottom nodes
    }
    //it's meaningless looking for cutoff range for binary X...

    size_t c;
    if(U>=L){ //randomly pick cutoff (can happen c==nx->c, v==nx->v)
      c = L + floor(unif_rand()*(U-L+1));
    } else {
      (*pA) = 0.0;
      L1.push_back(0);
      return false;
    }



    //make a copy xstar of the x, change xstar at node nxstar (same location as nx to x)
    tree xstar;
    tree::tree_p xp =  x.getptr(1); //pointer to x
    xstar.pcp(xp); //get a copy of x (by its pointer xp) to xstar
    tree::tree_p nxstar = xstar.getptr(nx->nid()); //point nstar to the same location in xstar as nx in xp
    nxstar->v = v;
    nxstar->c = c;

    //--------------------------------------------------
    //compute things needed for metropolis ratio

    double XLogPi, XLogL, YLogPi, YLogL;
    tree::npv nxbnv, nxstarbnv;
    YLogL = lilT1(X, xstar, nxstar, xi, pi, di, nxstarbnv, 1); // update bottom mu’s of the changed tree copy
    //is there empty bottom node(s)?
    if(YLogL == 10086.0){ //empty bn from rule change, invalid
      xstar.tonull();
      nxstar = 0;
      (*pA) = 0.0;
      L1.push_back(0);
      return false;
    }
    YLogPi = lpiT(nxstar, xi, pi, binaryX);

    XLogPi = lpiT(nx, xi, pi, binaryX);
    XLogL = lilT1(X, x, nx, xi, pi, di, nxbnv, 0); //do not update bottom mu’s

    double alpha;
    alpha = std::min(1.0,exp(YLogPi+YLogL-XLogPi-XLogL));

    if(unif_rand() < alpha) {

      (*pA) = 1.0;
      iP[nx->v]--;
      iP[nxstar->v]++;

      //replace xp (x) with xstar
      nx->v = nxstar->v; nx->c = nxstar->c;

      for(size_t k = 0; k != nxbnv.size(); k++){
        nxbnv[k]->mu = nxstarbnv[k]->mu;
      }

      xstar.tonull();
      nxstar = 0;
      
      L1.push_back(4);
      L2.push_back((double) nx->nid());
      L7.push_back((double) nx->v);L7.push_back((double) nx->c);
      L7.push_back((double) nxbnv.size());
      for(size_t k = 0; k != nxbnv.size(); k++){
        L8.push_back(nxbnv[k]->mu);
      }
      
      return true;
    } else {
      xstar.tonull();
      nxstar = 0;
      (*pA) = 0.0;
      L1.push_back(0);
      return false;
    }

  } //End if-else
  //PutRNGstate();
}
