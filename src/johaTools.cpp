// JohansenTools

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include "tools.h"

using namespace std;
using namespace arma;

// function returns the log-likelihood given...
double logLik(mat Z0, mat Z1, vec Z2, mat Psi, mat alpha, mat beta, mat Omega){
  int N = Z0.n_cols;
  mat Pi = alpha*beta.t();
  mat Omega1 = inv(Omega);
  mat tmp = Z0-Pi*Z1-Psi*Z2;
  // sumZ = sum of (Z0-ab'Z1-PsiZ2)'Om1(Z0-ab'Z1-PsiZ2) for t=1:T
  double sumZ = trace(tmp.t()*Omega1*tmp);
  double output = -0.5*N*log(det(Omega))-0.5*sumZ;
  return(output);
}

// Calculate M_{ij} matrices (Johansen procedure), input are two of Z0,Z1 or Z2
mat M(mat X,mat Y){
  int N = X.n_cols;
  int x = X.n_rows;
  int y = Y.n_rows;
  mat tmp = zeros<mat>(x,y);
  for(int n=0;n<N;n++){
    tmp += X.col(n)*Y.col(n).t();
  }

  return(tmp/N);
}

// Calculate S_{ij} matrices (Johansen procedure), input are two of R0 or R1
mat S(mat X, mat Y){
  int N = X.n_cols;
  int x = X.n_rows;
  int y = Y.n_rows;

  mat tmp = zeros<mat>(x,y);
  for(int n=0;n<N;n++){
    tmp += X.col(n)*Y.col(n).t();
  }
  return(tmp/N);
}


// Calculate the orthongonal complement to M
mat M_perp(mat M){
    int p = M.n_rows;
    int m = M.n_cols;
    mat Q,R,out;
    // QR decomposition of M
    qr(Q,R,M);
    if(p != m){
      // Orthogonal complent is p times p-m
      out = Q.cols(m,p-1);
    } else{
      out = zeros<mat>(1,1);
    }
    return(out);
}


// Calculate the normalization matrix of A
mat M_bar(mat M){
  mat out = M*inv(M.t()*M);
  return out;
}
