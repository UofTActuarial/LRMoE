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
double sumPoissonY(double mu, double lower_, double upper_)
{
  double result = 0;
  double lower = floor(lower_);
  double upper = ceil(upper_);

  if(!isinf(upper)){
    for(int j=0; j<upper-lower+1; j++){
      result = result + (lower+j)* exp(-mu + (lower+j)*log(mu)) / tgamma(lower+j+1);
    }
  }else{
    if(lower_==upper_){
      result = (lower)* exp(-mu + (lower)*log(mu)) / tgamma(lower+1);
    }else{
      double temp = 0;
      if(lower<=1){
        temp = 0.0;
      }else{
        for(int j=0; j<lower; j++){
          temp = temp + (j)* exp(-mu + (j)*log(mu)) / tgamma(j+1);
        }
      }
      result = mu - temp;
    }
  }

  return result;
}

////' @export
// [[Rcpp::export]]
SEXP sumPoissonYObs(SEXP mu_, SEXP lower_, SEXP upper_)
{
  NumericVector mu(mu_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
      result(j) = sumPoissonY(mu(0), lower(j), upper(j));
  }

  return result;
}

////' @export
// [[Rcpp::export]]
SEXP sumPoissonYLat(SEXP mu_, SEXP lower_, SEXP upper_)
{
  NumericVector mu(mu_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    double temp1 = sumPoissonY(mu(0), 0, lower(j)-1);
    double temp2 = sumPoissonY(mu(0), upper(j)+1, R_PosInf);
    result(j) = temp1 + temp2;
  }

  return result;
}
