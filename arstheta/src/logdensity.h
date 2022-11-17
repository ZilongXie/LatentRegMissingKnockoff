#ifndef __LOGDENS_INCLUDED__
#define __LOGDENS_INCLUDED__

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

class LogDensity {
  // Abstract base class for Log Densities
public:
  virtual NumericVector h(NumericVector x) = 0;
  virtual NumericVector h_prime(NumericVector x) = 0;
  virtual double h(double x) = 0;
  virtual double h_prime(double x) = 0;
};

class RLogDensity : public LogDensity {
  // Log Density class if functions are passed from R
public:
  Function R_h;
  Function R_h_prime;

  RLogDensity(Function R_h, Function R_h_prime);
  NumericVector h(NumericVector x);
  NumericVector h_prime(NumericVector x);
  double h(double x);
  double h_prime(double x);
};


class ThetaLogPost : public LogDensity {
  // Log Density of theta posterior
public:
  NumericVector Y;
  NumericVector avec;
  NumericVector d1vec;
  NumericVector d2vec;
  double mu;
  double var;
  
  ThetaLogPost(NumericVector Y, NumericVector avec, NumericVector d1vec, NumericVector d2vec, double mu, double var);
  
  NumericVector h(NumericVector theta);
  NumericVector h_prime(NumericVector theta);
  double h(double theta);
  double h_prime(double theta);
};

#endif
