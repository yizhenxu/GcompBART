#include <Rcpp.h>

#include <string>
#include <vector>
#include <map>
#include "tree.h"
#include "info.h"

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
using std::string;
//--------------------------------------------------
// constructors
tree::tree(): mu(0.0),v(0),c(0),p(0),l(0),r(0) {}
tree::tree(double m): mu(m),v(0),c(0),p(0),l(0),r(0) {}
tree::tree(const tree& n): mu(0.0),v(0),c(0),p(0),l(0),r(0) { cp(this,&n); }
//--------------------------------------------------
//operators
tree& tree::operator=(const tree& rhs)
{
  if(&rhs != this) {
    tonull(); //kill left hand side (this)
    cp(this,&rhs); //copy right hand side to left hand side
  }
  return *this;
}
//--------------------------------------------------
//public functions
// find bottom node pointer given x
//--------------------
tree::tree_cp tree::bn(double *x,xinfo& xi)
{
  if(l==0) return this; //bottom node
  if(x[v] < xi[v][c]) {
    return l->bn(x,xi);
  } else {
    return r->bn(x,xi);
  }
}
//--------------------
//find region for a given variable from root to this
void tree::rg(size_t v, int* L, int* U)
{
  if(p==0)  { //no parent
    return;
  }
  if(p->v == v) { //does my parent use v?
    if(this == p->l) { //am I left or right child
      if((int)(p->c) <= (*U)) *U = (p->c)-1;
    } else {
      if((int)(p->c) >= *L) *L = (p->c)+1;
    }
  }
  p->rg(v,L,U);
}
//--------------------
//tree size
size_t tree::treesize() const
{
  if(l==0) return 1;  //if bottom node, tree size is 1
  else return (1+l->treesize()+r->treesize());
}
//--------------------
size_t tree::nnogs() const
{
  if(!l) return 0; //bottom node
  if(l->l || r->l) { //not a nog
    return (l->nnogs() + r->nnogs());
  } else { //is a nog
    return 1;
  }
}
size_t tree::nuse(size_t v)
{
  npv nds;
  this->getnodes(nds);
  size_t nu=0; //return value
  for(size_t i=0;i!=nds.size();i++) {
    if(nds[i]->l && nds[i]->v==v) nu+=1;
  }
  return nu;
}
//--------------------
size_t tree::nbots() const
{
  if(l==0) { //if a bottom node
    return 1;
  } else {
    return l->nbots() + r->nbots();
  }
}
//--------------------
//depth of node
size_t tree::depth() const
{
  if(!p) return 0; //no parents
  else return (1+p->depth());
}
//--------------------
// node id
size_t tree::nid() const
  //recursion up the tree
{
  if(!p) return 1; //if you don't have a parent, you are the top
  if(this==p->l) return 2*(p->nid()); //if you are a left child
  else return 2*(p->nid())+1; //else you are a right child
}
//--------------------
//node type
char tree::ntype() const
{
  //t:top, b:bottom, n:no grandchildren, i:internal
  if(!p) return 't';
  if(!l) return 'b';
  if(!(l->l) && !(r->l)) return 'n';
  return 'i';
}
//--------------------
//get bottom nodes
//recursion down the tree
void tree::getbots(npv& bv)
{
  if(l) { //have children
    l->getbots(bv);
    r->getbots(bv);
  } else {
    bv.push_back(this);
  }
}
//--------------------
//get nog nodes
//recursion down the tree
void tree::getnogs(npv& nv)
{
  if(l) { //have children
    if((l->l) || (r->l)) {  //have grandchildren
      if(l->l) l->getnogs(nv);
      if(r->l) r->getnogs(nv);
    } else {
      nv.push_back(this);
    }
  }
}
//--------------------
//get all nodes
//recursion down the tree
void tree::getnodes(npv& v)
{
  v.push_back(this);
  if(l) {
    l->getnodes(v);
    r->getnodes(v);
  }
}
void tree::getnodes(cnpv& v)  const
{
  v.push_back(this);
  if(l) {
    l->getnodes(v);
    r->getnodes(v);
  }
}
//--------------------
//add children to  bot node nid
bool tree::birth(size_t nid,size_t v, size_t c, double ml, double mr)
{
  tree_p np = getptr(nid);
  if(np==0) {
    Rcpp::stop("error in birth: bottom node not found\n");
    return false; //did not find note with that nid
  }
  if(np->l) {
    Rcpp::stop("error in birth: found node has children\n");
    return false; //node is not a bottom node
  }

  //add children to bottom node np
  tree_p l = new tree;
  l->mu=ml;
  tree_p r = new tree;
  r->mu=mr;
  np->l=l;
  np->r=r;
  np->v = v; np->c=c;
  l->p = np;
  r->p = np;

  return true;
}
//--------------------
//is the node a nog node
bool tree::isnog() const
{
  bool isnog=true;
  if(l) {
    if(l->l || r->l) isnog=false; //one of the children has children.
  } else {
    isnog=false; //no children
  }
  return isnog;
}
//--------------------
//kill children of  nog node nid
bool tree::death(size_t nid, double mu)
{
  tree_p nb = getptr(nid);

  if(nb==0) {
    Rcpp::stop("error in death, nid invalid\n");
    return false;
  }
  if(nb->isnog()) {
    delete nb->l;
    delete nb->r;
    nb->l=0;
    nb->r=0;
    nb->v=0;
    nb->c=0;
    nb->mu=mu;
    return true;
  } else {
    Rcpp::stop("error in death, node is not a nog node\n");
    return false;
  }
}
//--------------------
//add children to bot node *np
void tree::birthp(tree_p np,size_t v, size_t c, double ml, double mr)
{
  tree_p l = new tree;
  l->mu=ml;
  tree_p r = new tree;
  r->mu=mr;
  np->l=l;
  np->r=r;
  np->v = v; np->c=c;
  l->p = np;
  r->p = np;
}
//--------------------
//kill children of  nog node *nb
void tree::deathp(tree_p nb, double mu)
{
  delete nb->l;
  delete nb->r;
  nb->l=0;
  nb->r=0;
  nb->v=0;
  nb->c=0;
  nb->mu=mu;
}
//--------------------------------------------------
//private functions
//--------------------
//copy tree o to tree n
void tree::cp(tree_p n, tree_cp o)
  //assume n has no children (so we don't have to kill them)
  //recursion down
{
  if(n->l) {
    Rcpp::stop("cp:error node has children\n");
    return;
  }

  n->mu = o->mu;
  n->v = o->v;
  n->c = o->c;

  if(o->l) { //if o has children
    n->l = new tree;
    (n->l)->p = n;
    cp(n->l,o->l);
    n->r = new tree;
    (n->r)->p = n;
    cp(n->r,o->r);
  }
}

void tree::pcp(tree_cp o)
  //assume this has no children (so we don't have to kill them)
  //recursion down
{
  if(l) {
    Rcpp::stop("cp:error node has children\n");
    return;
  }

  mu = o->mu;
  v = o->v;
  c = o->c;

  if(o->l) { //if o has children
    l = new tree;
    l->p = this;
    l->pcp(o->l);
    r = new tree;
    r->p = this;
    r->pcp(o->r);
  }
}

//--------------------
//cut back to one node
void tree::tonull()
{
  size_t ts = treesize();
  while(ts>1) { //if false ts=1
    npv nv;
    getnogs(nv);
    for(size_t i=0;i<nv.size();i++) {
      delete nv[i]->l;
      delete nv[i]->r;
      nv[i]->l=0;
      nv[i]->r=0;
    }
    ts = treesize();
  }
  mu=0.0;
  v=0;c=0;
  p=0;l=0;r=0;
}


//--------------------
// get pointer for node from its nid
tree::tree_p tree::getptr(size_t nid)
{
  if(this->nid() == nid) return this; //found it
  if(l==0) return 0; //no children, did not find it
  tree_p lp = l->getptr(nid);
  if(lp) return lp; //found on left
  tree_p rp = r->getptr(nid);
  if(rp) return rp; //found on right
  return 0; //never found it
}

//--------------------
//For SWAP
//--------------------
//get parent nodes of swappable pairs
//internal nodes that have children AND grandchildren
void tree::getswapdad(npv& sv)
{
  if(l) { //have children
    if((l->l) || (r->l)) {  //have grandchildren
      sv.push_back(this);
      if(l->l) l-> getswapdad(sv);
      if(r->l) r-> getswapdad(sv);
    }
  }
}

//--------------------
//get internal nodes for rule change
void tree::getintn(npv& iv)
{
  if(l) { //have children
    iv.push_back(this);
    l-> getintn(iv);
    r-> getintn(iv);
  }
}

//--------------------
//find region for a given variable from this to bottom nodes
//in bd, range of L & U from root to this is given to downrg
void tree::downrg(size_t v, int* L, int* U)
{
  int lmin,lmax,rmin,rmax;
  lmin=(*U)+1;
  rmin=(*U)+1;
  lmax=(*L)-1;
  rmax=(*L)-1;

  (this->l)->DownMinMax(v,&lmin,&lmax);
  (this->r)->DownMinMax (v,&rmin,&rmax);

  *L = (int)std::max((*L), lmax+1);
  *U = (int)std::min((*U), rmin-1);
}

//--------------------
//adjust L & U from this to bottom nodes
void tree::DownMinMax(size_t vv, int* L, int* U)
{
  if(l) {
    if(vv == (this->v) ) {
      if( (int)(this->c) < (*L) )  (*L) =  (this->c);
      if((int)(this->c) > (*U))  (*U) =  (this->c);
    }
    (this->l)->DownMinMax(vv, L, U);
    (this->r)->DownMinMax(vv, L, U);
  }

}


//--------------------
//x.Dbots(0,tree_depth): to get the number of bottom nodes of x on the level of tree_depth
//x.Dbots(d1, d2):say the current x has depth d1, get the number of bn of x on the level d1+d2
double tree::Dbots(int dnow, int td)
{
  if(l==0) {
    if(dnow == td) {
      return 1.0;
    } else {
      return 0.0;
    }
  } else return( l->Dbots(dnow+1, td) + r->Dbots(dnow+1, td) );
}
