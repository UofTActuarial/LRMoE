// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <Rcpp.h>
#include <RcppNumerical.h>
using namespace Rcpp;
using namespace Numer;

// For accessing constant isinf
#define _USE_MATH_DEFINES
#include <math.h>


// //' @export
// [[Rcpp::export]]
double sumBinomialY(double n, double p, double lower_, double upper_)
{
  double result = 0;
  double lower = fmax(floor(lower_), 0.0);
  double upper = fmin(ceil(upper_), n);

  if(!isinf(upper)){
    for(int j=0; j<upper-lower+1; j++){
      result = result + (lower+j)* exp((lower+j)*log(p) + (n-(lower+j))*log(1-p)) * tgamma(n+1) / (tgamma(lower+j+1) * tgamma(n-(lower+j)+1));
    }
  }else{
    if(lower==upper){
      result = (lower)* exp((lower)*log(p) + (n-(lower))*log(1-p)) * tgamma(n+1) / (tgamma(lower+1) * tgamma(n-(lower)+1));
    }else{
      double temp = 0;
      if(lower<=1){
        temp = 0.0;
      }else{
        for(int j=0; j<lower; j++){
          temp = temp + (j)* exp((j)*log(p) + (n-(j))*log(1-p)) * tgamma(n+1) / (tgamma(j+1) * tgamma(n-(j)+1));
        }
      }
      result = n*p - temp;
    }
  }

  return result;
}

////' @export
// [[Rcpp::export]]
SEXP sumBinomialYObs(SEXP n_, SEXP p_, SEXP lower_, SEXP upper_)
{
  NumericVector n(n_);
  NumericVector p(p_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    result(j) = sumBinomialY(n(0), p(0), lower(j), upper(j));
  }

  return result;
}

////' @export
// [[Rcpp::export]]
SEXP sumBinomialYLat(SEXP n_, SEXP p_, SEXP lower_, SEXP upper_)
{
  NumericVector n(n_);
  NumericVector p(p_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    double temp1 = sumBinomialY(n(0), p(0), 0, lower(j)-1);
    double temp2 = sumBinomialY(n(0), p(0), upper(j)+1, n(0));
    result(j) = temp1 + temp2;
  }

  return result;
}
