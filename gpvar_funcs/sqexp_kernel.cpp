// load Rcpp
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// ================================================================================================
// Listed below are functions to create different kernel matrices

// [[Rcpp::export]]
arma::mat LogGKcpp(arma::mat Xin, arma::mat Xstar, double sc = 1){
    int i,j;
    double ns = Xstar.n_cols;;
    double ni = Xin.n_cols;
    mat K = zeros<mat>(ns,ni);
    for (i = 0; i<ns; i++){
        for(j = 0; j<ni; j++){
            if(i==j){
                break;
            }
            else{
                K(i,j) = (-1/(2*sc)*sum(pow(Xstar.col(i)-Xin.col(j),2)));
            }
        }
    }
    return K + trans(K);
}


// [[Rcpp::export]]
arma::mat GKcpp(arma::mat Xin, arma::mat Xstar,double h = 1, double sigf = 1, double sc = 1){
  int i,j;
  double ns = Xstar.n_cols;
  double ni = Xin.n_cols;
  mat K = zeros<mat>(ns,ni);
  for (i = 0; i<ns; i++){
    for(j = 0; j<ni; j++){
      K(i,j) = sigf*exp(-h/(2*sc)*sum(pow(Xstar.col(i)-Xin.col(j),2)));
    }
  }
  return K;
}





