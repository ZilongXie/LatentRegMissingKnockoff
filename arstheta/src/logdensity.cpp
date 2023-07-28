#include <Rcpp.h>
#include "logdensity.h"
using namespace Rcpp;
using namespace std;


// RLogDensity implementation
// -> Log Density class if functions are passed from R
RLogDensity::RLogDensity(Function R_h, Function R_h_prime) :
    R_h(R_h), R_h_prime(R_h_prime) {}

NumericVector RLogDensity::h(NumericVector x) {
  NumericVector rv = R_h(x);
  return rv;
}

NumericVector RLogDensity::h_prime(NumericVector x) {
  NumericVector rv = R_h_prime(x);
  return rv;
}

double RLogDensity::h(double x) {
  NumericVector rv = R_h(x);
  return rv.at(0);
}

double RLogDensity::h_prime(double x) {
  NumericVector rv = R_h_prime(x);
  return rv.at(0);
}


// Theta posterior implementation
ThetaLogPost::ThetaLogPost(NumericVector Y, NumericVector avec, NumericVector d1vec, NumericVector d2vec, double mu, double var) :
  Y(Y), avec(avec), d1vec(d1vec), d2vec(d2vec), mu(mu), var(var) {}

NumericVector ThetaLogPost::h(NumericVector theta) {
  // Compute the logdensity, taking vector as input, return vector
  NumericVector rv = - 0.5*pow((theta-mu), 2)/var;
  for (int i = 0; i < Y.size(); ++i) {
    if (!NumericVector::is_na(Y.at(i))) {
      if (NumericVector::is_na(d2vec.at(i))) {
        rv += Y.at(i)*(avec.at(i)*theta + d1vec.at(i)) - log(1 + exp(avec.at(i)*theta + d1vec.at(i)));  
      }
      if (!NumericVector::is_na(d2vec.at(i))) {
        rv += -log(1 + exp(avec.at(i)*theta + d1vec.at(i)) + exp(2*avec.at(i)*theta + d1vec.at(i) + d2vec.at(i)));
        if (Y.at(i) >= 1) {
          rv += avec.at(i)*theta + d1vec.at(i);
        }
        if (Y.at(i) >= 2) {
          rv += avec.at(i)*theta + d2vec.at(i);
        }
      }
    }
  }
  return rv;
}

NumericVector ThetaLogPost::h_prime(NumericVector theta) {
  // Compute the detivative of logdensity, taking vector as input, return vector
  NumericVector rv = - (theta - mu)/var;
  for (int i = 0; i < Y.size(); ++i) {
    if (!NumericVector::is_na(Y.at(i))) {
      if (NumericVector::is_na(d2vec.at(i))) {
        rv += Y.at(i)*avec.at(i) - avec.at(i)/(1+exp(-(avec.at(i)*theta+d1vec.at(i))));
      }
      if (!NumericVector::is_na(d2vec.at(i))) {
        rv += Y.at(i)*avec.at(i) - avec.at(i)*(
          (exp(avec.at(i)*theta+d1vec.at(i)) + 2*exp(2*avec.at(i)*theta+d1vec.at(i)+d2vec.at(i)))
        /(1 + exp(avec.at(i)*theta+d1vec.at(i)) + exp(2*avec.at(i)*theta+d1vec.at(i)+d2vec.at(i))));
      }
    }
  }
  return rv;
}

double ThetaLogPost::h(double theta) {
  // Compute the logdensity, taking scalar as input, return scalar
  NumericVector rv = h(wrap(theta));
  return rv.at(0);
}

double ThetaLogPost::h_prime(double theta) {
  // Compute the detivative of logdensity, taking scalar as input, return scalar
  NumericVector rv = h_prime(wrap(theta));
  return rv.at(0);
}
