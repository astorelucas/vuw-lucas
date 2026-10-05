#include <Rcpp.h>
#include <vector>
#include <cmath>

using namespace Rcpp;

// =====================================================
// Union-Find
// =====================================================

class UnionFind {
public:
  std::vector<int> parent;
  std::vector<int> rank;

  UnionFind(int n) {
    parent.resize(n);
    rank.resize(n,0);
    for(int i=0;i<n;i++) parent[i]=i;
  }

  int find(int x){
    if(parent[x]!=x)
      parent[x]=find(parent[x]);
    return parent[x];
  }

  void unite(int a,int b){
    a=find(a);
    b=find(b);

    if(a==b) return;

    if(rank[a]<rank[b]){
      parent[a]=b;
    } else if(rank[a]>rank[b]){
      parent[b]=a;
    } else {
      parent[b]=a;
      rank[a]++;
    }
  }
};

// =====================================================
// Helpers
// =====================================================

inline int idx(int i,int j,int L){
  return i*L+j;
}

// =====================================================
// Energy
// =====================================================

// [[Rcpp::export]]
double ising_energy(IntegerMatrix spins, double J=1.0){

  int L = spins.nrow();
  double E = 0.0;

  for(int i=0;i<L;i++){
    for(int j=0;j<L;j++){

      int s = spins(i,j);

      int right = spins(i,(j+1)%L);
      int down  = spins((i+1)%L,j);

      E -= J * s * right;
      E -= J * s * down;
    }
  }

  return E;
}

// =====================================================
// Magnetization
// =====================================================

// [[Rcpp::export]]
double ising_magnetization(IntegerMatrix spins){

  int L = spins.nrow();
  double sum = 0.0;

  for(int i=0;i<L;i++){
    for(int j=0;j<L;j++){
      sum += spins(i,j);
    }
  }

  return sum/(L*L);
}

// =====================================================
// Single Metropolis sweep
// =====================================================

// [[Rcpp::export]]
IntegerMatrix metropolis_step(
    IntegerMatrix spins,
    double beta,
    double J=1.0){

  int L = spins.nrow();

  for(int n=0;n<L*L;n++){

    int i = floor(R::runif(0,L));
    int j = floor(R::runif(0,L));

    int s = spins(i,j);

    int nn =
      spins((i+1)%L,j) +
      spins((i-1+L)%L,j) +
      spins(i,(j+1)%L) +
      spins(i,(j-1+L)%L);

    double dE = 2.0 * J * s * nn;

    if(dE <= 0.0 ||
       R::runif(0,1) < std::exp(-beta*dE)){
      spins(i,j) = -s;
    }
  }

  return spins;
}

// =====================================================
// Heat-bath (Gibbs) sweep
// =====================================================

// [[Rcpp::export]]
IntegerMatrix gibbs_step(
    IntegerMatrix spins,
    double beta,
    double J=1.0){

  int L = spins.nrow();

  for(int i=0;i<L;i++){
    for(int j=0;j<L;j++){

      int h =
        spins((i+1)%L,j) +
        spins((i-1+L)%L,j) +
        spins(i,(j+1)%L) +
        spins(i,(j-1+L)%L);

      double p =
        1.0/(1.0 + std::exp(-2.0*beta*J*h));

      spins(i,j) =
        (R::runif(0,1) < p) ? 1 : -1;
    }
  }

  return spins;
}

// =====================================================
// Swendsen-Wang sweep
// =====================================================

// [[Rcpp::export]]
IntegerMatrix swendsen_wang_step(
    IntegerMatrix spins,
    double beta,
    double J=1.0){

  int L = spins.nrow();
  int N = L*L;

  UnionFind uf(N);

  double p = 1.0 - std::exp(-2.0*beta*J);

  for(int i=0;i<L;i++){
    for(int j=0;j<L;j++){

      int current = idx(i,j,L);

      int ir = i;
      int jr = (j+1)%L;

      int idr = idx(ir,jr,L);

      if(spins(i,j)==spins(ir,jr)){
        if(R::runif(0,1)<p){
          uf.unite(current,idr);
        }
      }

      int idown = idx((i+1)%L,j,L);

      if(spins(i,j)==spins((i+1)%L,j)){
        if(R::runif(0,1)<p){
          uf.unite(current,idown);
        }
      }
    }
  }

  std::vector<int> cluster_spin(N,0);
  std::vector<bool> assigned(N,false);

  for(int k=0;k<N;k++){

    int root = uf.find(k);

    if(!assigned[root]){
      cluster_spin[root] =
        (R::runif(0,1) < 0.5) ? 1 : -1;

      assigned[root]=true;
    }
  }

  IntegerMatrix out(L,L);

  for(int i=0;i<L;i++){
    for(int j=0;j<L;j++){

      int root =
        uf.find(idx(i,j,L));

      out(i,j) =
        cluster_spin[root];
    }
  }

  return out;
}

