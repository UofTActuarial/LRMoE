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
double sumZTPoissonY(double mu, double lower_, double upper_)
{
  double result = 0;
  double lower = fmax(floor(lower_), 0.0);
  double upper = ceil(upper_);

  if(!isinf(upper)){
    for(int j=0; j<upper-lower+1; j++){
      result = result + (lower+j)* exp(-mu + (lower+j)*log(mu) - lgamma(lower+j+1)) / (1-exp(-mu));
    }
  }else{
    if(lower==upper){
      result = (lower)* exp(-mu + (lower)*log(mu) - lgamma(lower+1)) / (1-exp(-mu));
    }else{
      double temp = 0;
      if(lower<=1){
        temp = 0.0;
      }else{
        for(int j=0; j<lower; j++){
          temp = temp + (j)* exp(-mu + (j)*log(mu) - lgamma(j+1)) / (1-exp(-mu));
        }
      }
      result = mu / (1-exp(-mu)) - temp;
    }
  }

  return result;
}

////' @export
// [[Rcpp::export]]
SEXP sumZTPoissonYObs(SEXP mu_, SEXP lower_, SEXP upper_)
{
  NumericVector mu(mu_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    result(j) = sumZTPoissonY(mu(0), lower(j), upper(j));
  }

  return result;
}

////' @export
// [[Rcpp::export]]
SEXP sumZTPoissonYLat(SEXP mu_, SEXP lower_, SEXP upper_)
{
  NumericVector mu(mu_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    double temp1 = sumZTPoissonY(mu(0), 0, lower(j)-1);
    double temp2 = sumZTPoissonY(mu(0), upper(j)+1, R_PosInf);
    result(j) = temp1 + temp2;
  }

  return result;
}
