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
arma::mat bootLoop(arma::mat X, int r, int B, double dt){
  int N = X.n_rows-1;
  int p = X.n_cols;
  // mat Omega = Minit.rows( 2*r+1     , 2*r+d     );
  // mat test  = Minit.rows( 2*r+d+1   , 2*r+d+1   );
  // mat eigs  = Minit.rows( 2*r+d+2   , 2*r+d+2   );

  arma::mat Minit, res, alpha, beta, Pi, Psi, S2, Xboot;
  // Get residuals from model estimation with r=0, r=1 and r=2
  if(r == 0){
    Minit   = var(X, dt);
    Pi      = zeros<mat>(p,p);
    // res     = Minit.rows(p+2+1, p+2+N);
    // Psi     = Minit.rows(0,0);
    // Psi     = Psi.t();
    // S2      = Minit.rows( 1     , p   );
  } else{
    Minit   = vecm(X,r,eye<mat>(p,p),eye<mat>(p,p), dt);
    alpha   = Minit.rows( 0         , r-1       );
    beta    = Minit.rows( r         , 2*r-1     );
    Pi      = alpha.t()*beta;
  }

  res     = Minit.rows( 2*r+p+2+1 , 2*r+p+2+N );
  Psi     = Minit.rows( 2*r       , 2*r       );
  Psi     = Psi.t();
  S2      = Minit.rows( 2*r+1     , 2*r+p   );
  res     = res.t();

  X       = X.t();
  Xboot   = zeros<mat>(p,N);

  mat resBoot, Mboot;
  mat out = zeros<mat>(B,1);

  for(int b=0;b<B;b++){
    resBoot = (res-repmat(mean(res,1),1,res.n_cols))%randn(p,N);
    for(int n=1;n<N;n++){
      Xboot.col(n) = Xboot.col(n-1) + (Pi*X.col(n-1) + Psi)*dt + chol(S2).t()*resBoot.col(n);
    }

    Mboot = vecm(Xboot.t(),1,eye<mat>(p,p),eye<mat>(p,p), dt);

    // cout << Mboot_1 << endl;
    // cout << Mboot_1(2*r+d+1,0) << " " << Mboot_1(2*r+d+1,1) << " " << Mboot_1(2*r+d+1,2) <<endl;
    out(b,0) =  Mboot(2*1+p+1,r);
  }
  // cout << Mboot << endl;
  return out;
}



// The the function returns ...
// [[Rcpp::export]]
arma::mat bootstrapCpp(arma::mat X, int B, double dt){

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
        out = zeros<mat>(B,p);
        for(int r=0;r<p;r++){
          boot = bootLoop(X,r,B,dt);
          //mat boot = zeros<mat>(B,d);
          out.col(r) = boot; //join_cols(est,boot);
        }
    } else{
      out = zeros<mat>(1,p);//est;
    }
    return out;
}
