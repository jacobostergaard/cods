// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef integralApprox_H
#define integralApprox_H

using namespace arma;
// This is the content of the .h file, which is where the declarations go
unsigned long long  factorialCpp(unsigned long long  x, unsigned long long  result = 1);
arma::mat kronSum(arma::mat A);
arma::mat matPow(arma::mat A, int k);
arma::vec integralApprox(arma::mat P, double t, bool verbose=false);
arma::mat Omega(arma::mat P,double t,double s);

// This is the end of the header guard
#endif
