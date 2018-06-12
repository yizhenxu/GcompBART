#ifndef GUARD_info_h
#define GUARD_info_h

//data
class dinfo {
public:
   dinfo() {n_cov=0;n_samp=0; n_dim = 0;y=0;}
   size_t n_samp;  //number of vars
   size_t n_cov;  //number of observations
   size_t n_dim;
   double *y;

};


//prior and mcmc
class pinfo
{
public:
   pinfo() {pbd=1.0;pb=.5;pswap=0.0;alpha=.95;beta=.5;tau=1.0;sigma=1.0;}
//mcmc info
   double pbd; //prob of birth/death
   double pb;  //prob of birth
   double pswap; //prob of swap
//prior info
   double alpha;
   double beta;
   double tau;
//sigma
   double sigma;
};

//sufficient statistics for 1 node
class sinfo
{
public:
   sinfo() {n=0;sy=0.0;sy2=0.0;}
   size_t n;
   double sy;
   double sy2;
};

#endif
