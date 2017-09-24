// Perform Johansen procedure for alpha/beta-restricted VECM models.
// Model is dX = (Pi*X[t-1]+Phi*D[t])dt + Sigma*dW[t], 3-dimensional... for now.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include "johaTools.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


arma::mat vecm(arma::mat X, int r, arma::mat A, arma::mat B, double dt){
  // Find Z0, Z1 and Z2 for input data matrix X
  int N = X.n_rows-1;   // Number of observations
  int p = X.n_cols;     // Dimension of the system
  //int m = A.n_cols;     // Restrictions for alpha matrix (unused)
  int s = B.n_cols;     // Restirctions for beta matrix

  X = X.t();            // Easier to work with column vectors
  mat Z0 = zeros<mat>(p,N);
  mat Z1 = zeros<mat>(p,N);
  mat Z2 = ones<mat>(1,N);

  for(int n=0;n<N;n++){
    Z0.col(n) = X.col(n+1)-X.col(n);
    Z1.col(n) = X.col(n);
  }


  // Find M_{ij} matrices
  mat M00,M11,M22,M01,M02,M12,M22_1;
  M00 = M(Z0,Z0);
  M11 = M(Z1,Z1);
  M22 = M(Z2,Z2);
  M01 = M(Z0,Z1);
  M02 = M(Z0,Z2);
  M12 = M(Z1,Z2);


  M22_1 = inv(M22);

  // Find residuals R0,R1
  mat R0,R1;
  R0 = Z0-M02*M22_1*Z2;
  R1 = Z1-M12*M22_1*Z2;

  // Find matrices S_{ij}
  mat S00,S11,S01,S10,S00_1;
  S00 = S(R0,R0);
  S11 = S(R1,R1);
  S01 = S(R0,R1);
  S10 = S(R1,R0);
  S00_1 = inv(S00);

  mat solveMat;
  mat cholS;

  // For restricted A: correct residuals
  mat Ap = M_perp(A);
  mat Ab = M_bar(A);
  mat Rt0,Rt1,St00,St11,St01,St10,St00_1;

  if(Ap.size() != 1){ // A.perp is not degenerate!
    Rt0 = R0-S00*Ap*inv(Ap.t()*S00*Ap)*Ap.t()*R0;
    Rt1 = R1-S10*Ap*inv(Ap.t()*S00*Ap)*Ap.t()*R0;

    St00 = S(Rt0,Rt0);
    St11 = S(Rt1,Rt1);
    St01 = S(Rt0,Rt1);
    St10 = S(Rt1,Rt0);
    St00_1 = inv(Ab.t()*S00*Ab);

    solveMat = inv(B.t()*St11*B)*B.t()*St10*Ab*St00_1*Ab.t()*St01*B;

    cholS = chol(B.t()*St11*B,"lower");

  } else{
    // Solve for the eigenvalues (and vectors)
    solveMat = inv(B.t()*S11*B)*B.t()*S10*S00_1*S01*B;
    cholS = chol(B.t()*S11*B,"lower");
  }

  cx_vec eigval_cx;
  cx_mat eigvec_cx;

  eig_gen(eigval_cx, eigvec_cx, solveMat);

  // C++ function returns complex vectors/matrices, so extract real parts
  vec eigval = real(eigval_cx);
  mat eigvec = real(eigvec_cx);

  // Sort by eigenvalues, descending (sort vectors first!)
  eigvec = eigvec.cols(sort_index(eigval,"descend"));
  eigval = eigval(sort_index(eigval,"descend"));

  // Normalize eigenvectors
  for(int i=0;i<s;i++){

    if(Ap.size() != 1){ // A.perp is not degenerate!
      double nf = as_scalar(sqrt(eigvec.col(i).t()*(B.t()*St11*B)*eigvec.col(i)));
      eigvec.col(i) = eigvec.col(i)/sqrt(nf);
    } else{
      double nf = as_scalar(sqrt(eigvec.col(i).t()*(B.t()*S11*B)*eigvec.col(i)));
      eigvec.col(i) = eigvec.col(i)/sqrt(nf);
    }
  }

  // To use cumsum for the teststats, the eigenvalues must be sorted "in reverse", this is ascending...
  vec testStat = -N*cumsum(log(1-eigval(sort_index(eigval,"ascend"))));
  //cout << testStat << endl;
  testStat = sort(testStat,"descend");
  //uvec idx = linspace<uvec>(0, r-1, 1);

  // Normalize beta using c {p x r} = (I{r x r},0{()p-r) x r})
  mat b_hat = B*eigvec.cols(0,r-1);
  mat Ir = eye<mat>(r,r);
  mat c = join_cols(Ir,zeros<mat>(p-r,r));
  //b_hat = b_hat/b_hat(0,0);
  b_hat = b_hat*inv(c.t()*b_hat);

  mat BS11B_1 = inv(b_hat.t()*S11*b_hat);

  mat a_hat = A*Ab.t()*S01*b_hat*BS11B_1;
  mat Psi_hat = (M02*M22_1-a_hat*b_hat.t()*M12*M22_1);
  mat Omega_hat = S00-S01*b_hat*BS11B_1*b_hat.t()*S10;

  // Calculate residuals
  mat Pi_hat = a_hat*b_hat.t();
  mat res = zeros<mat>(p,N);
  for(int n=0;n<N;n++){
    res.col(n) = Z0.col(n)-Pi_hat*(Z1.col(n))-Psi_hat;
  }

  int outRows = a_hat.n_cols+b_hat.n_cols+1+p+2+N;
  int outCols = p;
  mat out = zeros<mat>(outRows,outCols);

  // Insert alpha estimate in output
  for(int i=0; i<a_hat.n_cols;i++){
    out.row(i) = a_hat.col(i).t()/dt;
  }

  // Insert beta estimate in output
  for(int i=0;i<b_hat.n_cols;i++){
    int j = a_hat.n_cols;
    out.row(j+i) = b_hat.col(i).t();
  }
  // Insert Psi estimate in output
  int j = a_hat.n_cols+b_hat.n_cols;
  out.row(j) = Psi_hat.t()/dt;

  // Insert Omega estimate in output
  for(int i=0;i<p;i++){
    j = a_hat.n_cols+b_hat.n_cols+1;
    out.row(j+i) = Omega_hat.col(i).t()/dt;
  }

  // Insert test statistic in output
  // Append zeros to match dimensions...
  testStat = join_cols(testStat,zeros<mat>(p-s,1));
  eigval = join_cols(eigval,zeros<mat>(p-s,1));

  // cout << "Test " << testStat << endl;
  // cout << "Eigval " << eigval << endl;

  j = a_hat.n_cols+b_hat.n_cols+1+p;
  out.row(j) = testStat.t();
  j = a_hat.n_cols+b_hat.n_cols+1+p+1;
  out.row(j) = eigval.t();

  // Insert residuals in output
  for(int i=0;i<N;i++){
    j = a_hat.n_cols+b_hat.n_cols+1+p+2;
    out.row(j+i) = res.col(i).t();
  }
  return out;
}


// Estimation for a VAR model (if r=0...)
arma::mat var(arma::mat X, double dt){
  int N = X.n_rows-1;
  int p = X.n_cols;

  X = X.t();

  mat Z0 = zeros<mat>(p,N);
  mat Psi_hat = zeros<mat>(p,1);

  // Estimation using Least Squares
  // Formula: Psi_hat = T^-1*sum_{t=1}^T dX_t (this is 3x1 dimensional)
  for(int n=0;n<N;n++){
    Z0.col(n) = X.col(n+1)-X.col(n);
    Psi_hat += Z0.col(n);
  }
  Psi_hat = Psi_hat/N;

  // Calculate residuals and covariance estimator (see Lütkepohls book, p. 75)
  mat res = zeros<mat>(p,N);
  mat Omega = zeros<mat>(p,p);
  for(int n=0;n<N;n++){
    res.col(n) = Z0.col(n)-Psi_hat;
    Omega += res.col(n)*res.col(n).t();
  }
  Omega = (Omega/N)*(N/(N-1));


  // Calculate r-hypotheses statistics
  int r_tmp = 1;
  mat tmp = vecm(X.t(),r_tmp, eye<mat>(p,p), eye<mat>(p,p), dt);
  mat test  = tmp.rows(2*r_tmp+p+1,2*r_tmp+p+1);
  mat eigs  = tmp.rows(2*r_tmp+p+2,2*r_tmp+p+2);
  mat joha = join_cols(test,eigs);

  // Model estimates
  mat est = join_cols(Psi_hat.t()/dt,Omega/dt);

  // Add test statistics for r- hypotheses
  mat est_test = join_cols(est,joha);

  // Output estimates and residuals
  mat out  = join_cols(est_test,res.t());

  // Estimation using MLE
  // Not implemented yet.

  return out;
}



// Johansen estimation procedure, returns eigenvalues, eigenvectors
// [[Rcpp::export]]
arma::mat johansenCpp(arma::mat X, int r, arma::mat A, arma::mat B, double dt=1){
  int N = X.n_rows-1;   // Number of observations
  int p = X.n_cols;     // Dimension of the system

  mat out = zeros<mat>(N,p);
  if(r > 0){
    out = vecm(X,r,A,B,dt);
  } else{
    out = var(X,dt);
  }
  return out;
}
