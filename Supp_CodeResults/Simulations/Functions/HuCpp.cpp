#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::export]]
arma::mat evalGrad(arma::mat Yinner, arma::mat X, arma::mat beta, arma::mat missingindicator){
  
  int n = Yinner.n_rows; 
  int p = beta.n_rows; 
  int q = beta.n_cols;
  vec onemat; onemat.ones(p); 
  mat grad(p, q); grad.zeros(); 
  
  for(int i=0; i<n; ++i){
    grad += - (X(span(i,i), span::all).t()*(Yinner(span(i,i), span::all) - X(span(i,i), span::all)*beta))%(onemat*missingindicator(span(i,i), span::all));
  }
  
  return grad; 
}

//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::export]]
double evalObj(arma::mat Yinner, arma::mat X, arma::mat beta, arma::mat missingindicator){
  
  int n = Yinner.n_rows; 
  int p = beta.n_rows; 
  int q = beta.n_cols;
  double obj = 0.0; 
  
  for(int i=0; i<n; ++i){
    obj += accu((pow(Yinner(span(i,i), span::all) - X(span(i,i), span::all)*beta, 2))%(missingindicator(span(i,i), span::all))); 
  }
  
  return .5*obj; 
}



