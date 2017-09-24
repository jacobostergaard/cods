#include <RcppArmadillo.h>


// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace arma;



// [[Rcpp::export]]
arma::mat testCpp(arma::mat A) {
  cout << "Success!" << endl;
  double dt = 0.1;
  //arma::vec B = arma::vec(A.col(0));
  // arma::mat B(3,3, fill::randn);
  // B = B*dt;
  ///arma:mat B = sqrt(dt)*randn(3,4);
  // vec B = regspace< vec>(0,10);
  // cout << B << endl;

  //cout << "The sum of 3 and 4 is: " << add(3, 4) << endl;


  return A + dt * (A);
}

