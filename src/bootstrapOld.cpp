// Bootstrapping for critical values for Johansens rank test for a VECM
// Model is dX = (Pi*X[t-1]+Phi*D[t])dt + Sigma*dW[t], 3-dimensional...

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include "tools.h"
#include "johansen.h"

using namespace Rcpp;
using namespace arma;
using namespace std;



// Bootstrapping loop, using the Johansen procedure for estimation of statistic-distribution (only for p=3 at the moment...)
arma::mat bootLoop(arma::mat X, int B, double dt){

  int N = X.n_rows-1;
  int p = X.n_cols;
  int r;
  // mat Omega = Minit.rows( 2*r+1     , 2*r+d     );
  // mat test  = Minit.rows( 2*r+d+1   , 2*r+d+1   );
  // mat eigs  = Minit.rows( 2*r+d+2   , 2*r+d+2   );

  arma::mat Minit, res, alpha, beta, Pi_0, Pi_1, Pi_2, Psi_0, Psi_1, Psi_2, S2_0, S2_1, S2_2, res_0, res_1, res_2, X_0, X_1, X_2, Xboot_0, Xboot_1, Xboot_2;
  // Get residuals from model estimation with r=0, r=1 and r=2
  r = 0;

  Minit   = var(X, dt);
  res     = Minit.rows(p+2+1, p+2+N);

  Pi_0    = zeros<mat>(p,p);
  Psi_0   = Minit.rows(0,0);
  Psi_0   = Psi_0.t();
  S2_0    = Minit.rows( 1     , p   );
  res_0   = res.t();

  X_0     = X.t();
  Xboot_0 = zeros<mat>(p,N);



  r = 1;
  Minit   = vecm(X,r,eye<mat>(p,p),eye<mat>(p,p), dt);
  alpha   = Minit.rows( 0         , r-1       );
  beta    = Minit.rows( r         , 2*r-1     );
  res     = Minit.rows( 2*r+p+2+1 , 2*r+p+2+N );

  Pi_1    = alpha.t()*beta;
  Psi_1   = Minit.rows( 2*r       , 2*r       );
  Psi_1   = Psi_1.t();
  S2_1    = Minit.rows( 2*r+1     , 2*r+p   );
  res_1   = res.t();

  X_1     = X.t();
  Xboot_1 = zeros<mat>(p,N);



  r = 2;
  Minit   = vecm(X,r,eye<mat>(p,p),eye<mat>(p,p), dt);
  alpha   = Minit.rows( 0         , r-1       );
  beta    = Minit.rows( r         , 2*r-1     );
  res     = Minit.rows( 2*r+p+2+1 , 2*r+p+2+N );

  Pi_2    = alpha.t()*beta;
  Psi_2   = Minit.rows( 2*r       , 2*r       );
  Psi_2   = Psi_2.t();
  S2_2    = Minit.rows( 2*r+1     , 2*r+p   );
  res_2   = res.t();


  X_2     = X.t();
  Xboot_2 = zeros<mat>(p,N);

  mat resBoot_0,resBoot_1,resBoot_2,Mboot_0,Mboot_1,Mboot_2;
  mat out = zeros<mat>(B,p);

  r=1;

  for(int b=0;b<B;b++){
    resBoot_0 = (res_0-repmat(mean(res_0,1),1,res_0.n_cols))%randn(p,N);
    resBoot_1 = (res_1-repmat(mean(res_1,1),1,res_1.n_cols))%randn(p,N);
    resBoot_2 = (res_2-repmat(mean(res_2,1),1,res_2.n_cols))%randn(p,N);
    for(int n=1;n<N;n++){
      Xboot_0.col(n) = Xboot_0.col(n-1) + (Pi_0*X_0.col(n-1) + Psi_0)*dt + chol(S2_0).t()*resBoot_0.col(n);
      Xboot_1.col(n) = Xboot_1.col(n-1) + (Pi_1*X_1.col(n-1) + Psi_1)*dt + chol(S2_1).t()*resBoot_1.col(n);
      Xboot_2.col(n) = Xboot_2.col(n-1) + (Pi_2*X_2.col(n-1) + Psi_2)*dt + chol(S2_2).t()*resBoot_2.col(n);
    }

    Mboot_0 = vecm(Xboot_0.t(),1,eye<mat>(p,p),eye<mat>(p,p), dt);
    Mboot_1 = vecm(Xboot_1.t(),1,eye<mat>(p,p),eye<mat>(p,p), dt);
    Mboot_2 = vecm(Xboot_2.t(),1,eye<mat>(p,p),eye<mat>(p,p), dt);
    // cout << Mboot_1 << endl;
    // cout << Mboot_1(2*r+d+1,0) << " " << Mboot_1(2*r+d+1,1) << " " << Mboot_1(2*r+d+1,2) <<endl;
    out(b,0) =  Mboot_0(2*r+p+1,0);
    out(b,1) =  Mboot_1(2*r+p+1,1);
    out(b,2) =  Mboot_2(2*r+p+1,2);
  }

  return out;
}



// The the function returns ...
// [[Rcpp::export]]
arma::mat bootstrapCppOld(arma::mat X, int B, double dt){

  // Include support for restricted alpha/beta matrices!
  //int N = X.n_rows;
  int p = X.n_cols;
  mat est,boot,out;

  // if(r>0){
  //   est = vecm(X,r,eye<mat>(p,p),eye<mat>(p,p), dt);
  // } else {
  //   est = var(X, dt);
  // }
  if(B > 0){
    // Run bootstrap method to determine test-statistic distributions!
    boot = bootLoop(X,B,dt);
    //mat boot = zeros<mat>(B,d);
    out = boot;//join_cols(est,boot);
  } else{
    out = zeros<mat>(1,p);//est;
  }


  return out;
}
