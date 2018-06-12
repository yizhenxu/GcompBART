#ifndef GUARD_tree_h
#define GUARD_tree_h

#include <iostream>
#include <cstddef>
#include <vector>

#include "info.h"


//xinfo xi, then xi[v][c] is the c^{th} cutpoint for variable v.
// left if x[v] < xi[v][c]
typedef std::vector<double> vec_d; //double vector
typedef std::vector<vec_d> xinfo; //vector of vectors, will be split rules

class tree {
public:
  //------------------------------


  friend bool bd(Rcpp::NumericMatrix& X, tree& x, xinfo& xi, dinfo& di, pinfo& pi, size_t minobsnode, Rcpp::NumericVector& binaryX,
                 std::vector<double>& L1,std::vector<double>& L2,std::vector<double>& L3,std::vector<double>& L4,
                 std::vector<double>& L5,std::vector<double>& L6,std::vector<double>& L7,std::vector<double>& L8);
  friend bool bd1(Rcpp::NumericMatrix& X, tree& x, xinfo& xi, dinfo& di, pinfo& pi, size_t minobsnode, Rcpp::NumericVector& binaryX, double* nLtD, double* pA, double* nN, double* nL, double* tD, std::vector<double>& iP,
                  std::vector<double>& L1,std::vector<double>& L2,std::vector<double>& L3,std::vector<double>& L4,
                  std::vector<double>& L5,std::vector<double>& L6,std::vector<double>& L7,std::vector<double>& L8);


  //------------------------------
  //typedefs
  typedef tree* tree_p;
  typedef const tree* tree_cp;
  typedef std::vector<tree_p> npv; //Node Pointer Vector
  typedef std::vector<tree_cp> cnpv; //const Node Pointer Vector

  //------------------------------
  //tree constructors, destructors
  tree();
  tree(const tree&);
  tree(double);
  ~tree() {tonull();}

  //------------------------------
  //operators
  tree& operator=(const tree&);

  //------------------------------
  //node access
  //you are freely allowed to change mu, v, c
  //set----------
  void setm(double mu) {this->mu=mu;}
  void setv(size_t v) {this->v = v;}
  void setc(size_t c) {this->c = c;}
  //get----------
  double getm() const {return mu;}
  size_t getv() const {return v;}
  size_t getc() const {return c;}
  tree_p getp() const {return p;}  //should this be tree_cp?
  tree_p getl() const {return l;}
  tree_p getr() const {return r;}

  //------------------------------
  //tree functions
  //stuff about the tree----------
  size_t treesize() const; //number of nodes in tree
  size_t nnogs() const;    //number of nog nodes (no grandchildren nodes)
  size_t nbots() const;    //number of bottom nodes
  double Dbots(int dnow, int td); //number of bottom nodes at depth td-dnow

  //birth death using nid----------
  bool birth(size_t nid,size_t v, size_t c, double ml, double mr);
  bool death(size_t nid, double mu);
  //vectors of node pointers----------
  void getbots(npv& bv);         //get bottom nodes
  void getnogs(npv& nv);         //get nog nodes (no granchildren)
  void getnodes(npv& v);         //get vector of all nodes
  void getnodes(cnpv& v) const;  //get all nodes
  void getswapdad(npv& sv); //get all swappable dad nodes
  void getintn(npv& iv); //get internal nodes (non-bn nodes) for rule change

  //find node from x and region for var----------
  tree_cp bn(double *x,xinfo& xi);  //find bottom node for x
  void rg(size_t v, int* L, int* U); //find region [L,U] for var v at this (top to this).
  size_t nuse(size_t v); //how many times var v is used in a rule.
  void downrg(size_t v, int* L, int* U); //find region for var v at this (top to bottom nodes)
  void DownMinMax(size_t vv, int* L, int* U); //used in downrg, adjust L & U by subtree paths of this

  //------------------------------
  //node functions
  size_t depth() const; //depth of a node
  size_t nid() const;   //node id
  char ntype() const;   //t:top;b:bottom;n:nog;i:interior, carefull a t can be bot
  tree_p getptr(size_t nid); //get node pointer from node id, 0 if not there.
  bool isnog() const;
  void tonull(); //like a "clear", null tree has just one node

  //------------------------------
  //parameter for node
  double mu;
  //------------------------------
  //rule: left if x[v] < xinfo[v][c]
  size_t v;
  size_t c;
  //------------------------------
  //tree structure
  tree_p p; //parent
  tree_p l; //left child
  tree_p r; //right child
  //------------------------------
  //utility functions
  void cp(tree_p n,  tree_cp o); //copy tree o to n
  void pcp(tree_cp o);//copy tree o to current pointer (assumed empty)

  void birthp(tree_p np,size_t v, size_t c, double ml, double mr);
  void deathp(tree_p nb, double mu); //kill children of nog node nb
};


#endif
