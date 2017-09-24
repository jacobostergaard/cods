#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include <Rcpp.h>

using namespace std;
using namespace arma;

// Return factorial of input x
// [[Rcpp::export]]
unsigned long long  factorialCpp(unsigned long long  x, unsigned long long result = 1) {
  if (x == 1) return result; else return factorialCpp(x - 1, x * result);
}


// Return Kronecker SUM of matrix A
// [[Rcpp::export]]
arma::mat kronSum(arma::mat A){
  int d = A.n_rows;               // Dimension of the matrix
  arma::mat I(d,d,fill::eye);     // Create identity matrix of size d
  arma::mat out = arma::kron(A,I)+arma::kron(I,A);  // Kronecker sum definition
  return out;
}

// Return A^k = AA...A and not A = {a_ij}^k
// [[Rcpp::export]]
arma::mat matPow(arma::mat A, int k){
  int d = A.n_rows;                 // Dimension of the matrix
  arma::mat out(d,d,fill::eye);     // Create identity matrix of size d
  for(int n=0; n<k; n++){
    out = out*A;
  }

  // Use tmp*A/k instead... then A^k/k! = (A/1)*(A/2)*...*(A/k)
  return out;
}


// Function returns approximate integral from 0 to t of exp[(t-s)*A] wrt s, where A is an arbitrary quadratic matrix.
// [[Rcpp::export]]
arma::mat integralApprox(arma::mat P, double t, bool verbose=false){
  double tol=1e-7;
  int d = P.n_rows;               // Dimension of matrix
  int maxIter = 50;              // ensure the loop doesn't run wild...!
  maxIter = 10; // c++ is imprecise with large powers...!
  arma::mat out(d,d,fill::zeros);
  //arma::mat A = kronSum(P);
  arma::mat pre = out;
  int k = 0;
  double chg = 2*tol; // Initiate with a change larger than the tolerance level.

  while((chg > tol) && (k < maxIter)){
    out = pre+matPow(P,k)/factorialCpp(k+1)*pow(t,k+1);
    chg = arma::norm(abs(out-pre),"inf"); // Check against the maximum of the differences

    if(verbose){
      cout << "-------------------------------------------------------------------" << k << endl;
      cout << "k=" << k << endl;
      cout << "chg= " << chg << endl;
      cout << "t=" << t << endl;
      cout << "factor is " << factorialCpp(k+1) << endl;
      cout << round(1e4*matPow(P*t,k))/1e4 << endl;
      cout << round(1e4*out)/1e4 << endl;
    }

    pre = out;
    k += 1;
  }

  return out;
}

// [[Rcpp::export]]
arma::mat Omega(arma::mat P,double t,double s){
  int d = P.n_rows;
  arma::mat I3(3, 3, fill::eye);
  arma::mat ZZ = arma::vectorise(s*I3);

  arma::mat out(d,d,fill::zeros);
  arma::mat intVal = integralApprox(kronSum(P),t)*ZZ; // This is a 9x1 vector
  out.col(0) = intVal.rows(0,2);
  out.col(1) = intVal.rows(3,5);
  out.col(2) = intVal.rows(6,8);

  // Convert back from vector to matrix


  return out;
}
