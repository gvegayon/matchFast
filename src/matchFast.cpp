#include <Rcpp.h>
#include <utility>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector stl_partial_sort(NumericVector x, int n=4) {
  NumericVector y = clone(x);
  if (n > x.size()) n = x.size();
  std::nth_element(y.begin(), y.begin()+n, y.end());
  return y;
}


// [[Rcpp::export]]
double pnorm(const NumericVector & X0, const NumericVector & X1, double p=2.0) {
  return pow(sum(pow(X0 - X1, p)), 1/p);
}


// [[Rcpp::export]]
NumericMatrix PNORM(const NumericMatrix & X, double p=2.0) {

  int n=X.nrow();
  NumericMatrix ans(n,n);
  for (int i=0;i<n;i++)
    for (int j=(i+1);j<n;j++) {
      ans(j,i) = pnorm(X(i,_), X(j,_), p);
      ans(i,j) = ans(j,i);
    }

  return ans;
}

// Functions to find the nth closests individuals
typedef std::pair<int, double> mypair;
bool comparator ( const mypair& l, const mypair& r) {
  return l.second < r.second;
}


void print_vector(std::vector<mypair> & v) {
  for (std::vector<mypair>::iterator it=v.begin(); it!=v.end();it++)
    std::cout << it->first << "\t" << it->second << '\n';
  return;
}

// [[Rcpp::export]]
NumericMatrix sort_index(NumericVector x, int id, int n=4) {

  // Creating objects
  std::vector< mypair > xind(x.size());
  std::vector<int> ans_index;
  std::vector<double> ans_value;

  // Filling the data
  for (int i=0;i<x.size();i++)
    xind[i] = std::make_pair(i, x[i]);

  // Sorting
  std::partial_sort(xind.begin(), xind.begin()+n+1, xind.end(), comparator);

  int i = 0;
  int count = 1;
  for (std::vector<mypair>::iterator it=xind.begin(); it!=xind.end(); it++) {
    // Same id doesnt count
    if (count == n) break;
    if ((it->first) == id) continue;

    // Filling
    ans_index.push_back((it->first)+1);
    ans_value.push_back(it->second);

    // Should it continue?
    if ((i>0) && (it->second != ans_value[i-1])) ++count;

    i++;

  }

  // Creating the matrix
  NumericMatrix ans(ans_index.size(),2);
  for (int i=0;i<ans_index.size();i++)
    ans(i,0) = ans_index[i], ans(i,1) = ans_value[i];

  return ans;
}

// [[Rcpp::export]]
List matchFast(const NumericMatrix & X, int m=2, double p=2.0) {
  int n = X.nrow();
  int k = X.nrow();
  List ans(n);

  // Computing distances
  NumericMatrix D = PNORM(X, p);

  // Matching
  for (int i=0;i<n;i++)
    ans[i] = sort_index(D(_,i), i, m);

  return ans;
}


/***R
set.seed(1231231)
n <- 1e3
k <- 10
X <- matrix(runif(n*k), ncol=k)
X <- rbind(X, X[1:2,], X[1:2,])
Y <- runif(nrow(X))

# D <- matrix(ncol=n, nrow=n)
# for (i in 1:(n-1))
#   for (j in (i+1):n)
#     D[i,j] = pnorm(X[i,], X[j,])
# D

# dist(X)
# PNORM(X)

matches <- matchFast(X, 100)
sapply(matches, nrow)
Xmatch <- sapply(seq_len(nrow(X)), function(i) {
  sum(matches[[i]][,2]*Y[matches[[i]][,1]])/sum(matches[[i]][,2])
})

plot(cbind(Xmatch, Y))
*/
