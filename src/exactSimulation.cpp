/* Main function winfreeOscillator integrates a stochastically oscillating system,
* based on linearly cointegrated phases.
*
  * This version :  07.07.2016
*
  * Auxillary functions used (these are found in other files):
  *  phase()       : return the phase of input of a set of (x,y) coordinates
*  amplitude()   : return the amplitude of input of a set of (x,y) coordinates
*  polarToXT()   : return (x,y) coordinates from polar input (gam,phi)
*  f()           : return the deterministic trend (including the cointegration relations) for the phase process
*  g()           : return the deterministic trend for the amplitude process
*/


  // [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "tools.h"
#include "f.h"
#include "g.h"
#include "integralApprox.h"

using namespace std;
using namespace arma;



// [[Rcpp::export]]
arma::mat exactSimulation(int N,              // Number of simulations
                          double dt,           // Discretization timestep for gamma process
                          arma::vec z0,       // Initial values of process
                          arma::vec freq,     // Intrinsic frequency
                          arma::vec lvl,      // Level for gamma process
                          double s_phi,    // Noise for phi
                          double s_gam,    // Noise for gamma
                          const char* model
                          ){

  // Idenfity the number of oscillators
  int p = 3;
  double t = 0;
  arma::vec xi_0, xi_1, xi_2, xi_3, mu;
  arma::mat Om_0, Om_1, Om_2, Om_3;
  arma::mat Pi_0, Pi_1, Pi_2, Pi_3;
  arma::mat I3(3, 3, fill::eye);
  // arma::mat ZZ = arma::vectorise(s_phi*I3);

  // Define Pi_x matrices:
  Pi_0 = zeros<mat>(3,3);
    arma::mat a,b;
    a = b = zeros(3,1);
    a(0,0) = -0.5;
    b(0,0) = 1;
    b(1,0) = -1;
  Pi_1 = a*b.t();
    a(1,0) = 0.5;
  Pi_2 = a*b.t();
    a = b = zeros(3,2);
    a(0,0) = a(1,1) = -0.5;
    a(0,1) = a(1,0) = a(2,0) = a(2,1) = 0.25;
    b(0,0) = b(1,1) = 1;
    b(2,0) = b(2,1) = -1;
  Pi_3 = a*b.t();

  // Initialize processes, note that the gamma process is shared!
  // amplitude process
  arma::mat gam(p,N+1, fill::zeros);
  gam.col(0)  = amplitude(z0);

  // phase process
  arma::mat phi0(p,N+1, fill::zeros);
  arma::mat phi1(p,N+1, fill::zeros);
  arma::mat phi2(p,N+1, fill::zeros);
  arma::mat phi3(p,N+1, fill::zeros);
  phi0.col(0)  = phase(z0);
  phi1.col(0)  = phase(z0);
  phi2.col(0)  = phase(z0);
  phi3.col(0)  = phase(z0);

  // observed process
  arma::mat Z0(2*p,N+1, fill::zeros);
  arma::mat Z1(2*p,N+1, fill::zeros);
  arma::mat Z2(2*p,N+1, fill::zeros);
  arma::mat Z3(2*p,N+1, fill::zeros);
  Z0.col(0)    = z0;
  Z1.col(0)    = z0;
  Z2.col(0)    = z0;
  Z3.col(0)    = z0;

  // Draw noise processes for phi and gamma and scale appropriately
  arma::mat dW_phi = randn(p,N+1);
  arma::mat dW_gam = sqrt(dt)*randn(p,N+1);

  // Main loop for iterations
  for(int n=1; n<N+1; n++){

    // Iterate gamma process
    gam.col(n) = gam.col(n-1)+g(gam.col(n-1), lvl, model)*dt+s_gam*dW_gam.col(n);
    // Fix possible negative gamma values
    gam.col(n) = abs(gam.col(n));

    // Exact phi process, given gam(n).

    mu = freq+gam.col(n);
    t = n*dt;
    xi_0 = phi0.col(0)+mu*t;
    xi_1 = phi1.col(0)+mu*t-2*(exp(-t/2)-1)*Pi_1*phi1.col(0) + 4*(exp(-t/2)+t/2-1)*Pi_1*mu;
    xi_2 = phi2.col(0)+mu*t-(exp(-t)-1)*Pi_2*phi2.col(0) + (exp(-t)+t-1)*Pi_2*mu;
    xi_3 = phi3.col(0)+mu*t-4/3*(exp(-3*t/4)-1)*Pi_3*phi3.col(0) + pow(4/3,2)*(exp(-3*t/4)+3*t/4-1)*Pi_3*mu;
    // Om_0 = s_phi*t*I3;
    // Om_1 = s_phi*t*I3;//+2*pow(s_phi,2)*Pi_1*(exp(-t)+t-1);
    // Om_2 = s_phi*t*I3;//+pow(s_phi,2)/2*Pi_2*(exp(-2*t+2*t-1))*Pi_2;
    // Om_3 = s_phi*t*I3;//+8/9*pow(s_phi,2)*Pi_3*(exp(-3*t/2)+3*t/2-1);

    // cout << "t=" << t << endl;
    // cout << "-----------" << endl;
    Om_0 = Omega(Pi_0, t, s_phi);
    // cout << round(1e4*Om_0)/1e4 << endl;
    // cout << chol(Om_0) << endl;
    // cout << "-----------" << endl;
    // Om_1 = Omega(Pi_1, t, s_phi);
    Om_1 = Omega(Pi_1, t, s_phi);
    // cout << round(1e4*Om_1)/1e4 << endl;
    // cout << chol(Om_1) << endl;
    // cout << "-----------" << endl;
    //Om_2 = Omega(Pi_2, t, s_phi);
    Om_2 = Omega(Pi_2, t, s_phi);
    // cout << round(1e4*Om_2)/1e4 << endl;
    // cout << chol(Om_2) << endl;
    // cout << "-----------" << endl;
    // Om_3 = Omega(Pi_3, t, s_phi);
    Om_3 = Omega(Pi_3, t, s_phi);
    // cout << round(1e4*Om_3)/1e4 << endl;
    // cout << chol(Om_3) << endl;
    // cout << "-----------" << endl;

    // cout << "t=" << t << endl;
    // cout << "-----------" << endl;
    // cout << Om_0 << endl;
    // cout << "-----------" << endl;
    // cout << Om_1 << endl;
    // cout << "-----------" << endl;
    // cout << Om_2 << endl;
    // cout << "-----------" << endl;
    // cout << Om_3 << endl;
    // cout << "-----------" << endl;


    // Iterate phi process
    phi0.col(n) = xi_0+chol(Om_0).t()*dW_phi.col(n);
    phi1.col(n) = xi_1+chol(Om_1).t()*dW_phi.col(n);
    phi2.col(n) = xi_2+chol(Om_2).t()*dW_phi.col(n);
    phi3.col(n) = xi_3+chol(Om_3).t()*dW_phi.col(n);
    // Find observed Z-process
    Z0.col(n) = polarToXY(gam.col(n), phi0.col(n));
    Z1.col(n) = polarToXY(gam.col(n), phi1.col(n));
    Z2.col(n) = polarToXY(gam.col(n), phi2.col(n));
    Z3.col(n) = polarToXY(gam.col(n), phi3.col(n));
  }

  // Output is Z process + timesteps, 2*p+1 rows, N+1 columns
  vec t_out = regspace<vec>(0,N);
  t_out = dt*t_out;
  arma::mat Z;
  Z = join_cols(Z0,Z1);
  Z = join_cols(Z,Z2);
  Z = join_cols(Z,Z3);
  Z.insert_rows(0,t_out.t());

  return Z.t();
}

// THE END
