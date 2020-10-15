// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <Rcpp.h>
#include <RcppNumerical.h>
using namespace Rcpp;
using namespace Numer;

// For accessing constant isinf
#define _USE_MATH_DEFINES
// #include <math.h>
#include <cmath>

// For qnbinom
#include <Rmath.h>


// //' @export
// [[Rcpp::export]]
double sumNegativeBinomialY(double n, double p, double lower_, double upper_)
{
  double result = 0;
  double lower = fmax(floor(lower_), 0.0);
  double upper = ceil(upper_);

  if(!std::isinf(upper)){
    for(int j=0; j<upper-lower+1; j++){
      result = result + (lower+j)* exp( (n)*log(p) + (lower+j)*log(1-p) + lgamma(n+(lower+j)) - lgamma(n) - lgamma((lower+j)+1) );
        // result + (lower+j)* exp((n)*log(p) + (lower+j)*log(1-p)) * tgamma(n+(lower+j)) / (tgamma(n) * tgamma((lower+j)+1));
    }
  }else{
    if(lower==upper){
      result = (lower)* exp((n)*log(p) + (lower)*log(1-p) + lgamma(n+(lower)) - lgamma(n) - lgamma((lower)+1) );
        // (lower)* exp((n)*log(p) + (lower)*log(1-p)) * tgamma(n+(lower)) / (tgamma(n) * tgamma((lower)+1));
    }else{
      double temp = 0;
      if(lower<=1){
        temp = 0.0;
      }else{
        for(int j=0; j<lower; j++){
          temp = temp + (j)* exp((n)*log(p) + (j)*log(1-p) + lgamma(n+(j)) - lgamma(n) - lgamma((j)+1) );
            // temp + (j)* exp((n)*log(p) + (j)*log(1-p)) * tgamma(n+(j)) / (tgamma(n) * tgamma((j)+1));
        }
      }
      result = n*(1-p)/p - temp;
    }
  }

  // if(isnan(result)){
  //   result = 0.0;
  // }
  return result;
}

////' @export
// [[Rcpp::export]]
SEXP sumNegativeBinomialYObs(SEXP n_, SEXP p_, SEXP lower_, SEXP upper_)
{
  NumericVector n(n_);
  NumericVector p(p_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    result(j) = sumNegativeBinomialY(n(0), p(0), lower(j), upper(j));
  }

  return result;
}

////' @export
// [[Rcpp::export]]
SEXP sumNegativeBinomialYLat(SEXP n_, SEXP p_, SEXP lower_, SEXP upper_)
{
  NumericVector n(n_);
  NumericVector p(p_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    double temp1 = sumNegativeBinomialY(n(0), p(0), lower(j), upper(j));
      // sumNegativeBinomialY(n(0), p(0), 0, lower(j)-1);
    // double temp2 = sumNegativeBinomialY(n(0), p(0), upper(j)+1, n(0));
    result(j) = n(0) * (1-p(0)) / p(0) - temp1; // temp1 + temp2;
  }

  return result;
}

///////////////////////////////////////////////////////////////////////////////////

// //' @export
// [[Rcpp::export]]
double sumNegativeBinomialLfacY(double n, double p, double nn, double lower_, double upper_)
{
  double result = 0;
  double lower = fmax(floor(lower_), 0.0);
  double upper = ceil(upper_);

  // if(!isinf(upper)){
    for(int j=0; j<upper-lower+1; j++){
      result = result + lgamma(lower+j+nn)* exp( (n)*log(p) + (lower+j)*log(1-p) + lgamma(n+(lower+j)) - lgamma(n) - lgamma((lower+j)+1) );
    }
  // }else{
  //   if(lower==upper){
  //     result = lgamma(lower+nn)* exp((n)*log(p) + (lower)*log(1-p) + lgamma(n+(lower)) - lgamma(n) - lgamma((lower)+1) );
  //   }else{
  //     double temp = 0;
  //     if(lower<=1){
  //       temp = 0.0;
  //     }else{
  //       for(int j=0; j<lower; j++){
  //         temp = temp + lgamma(j+nn)* exp((n)*log(p) + (j)*log(1-p) + lgamma(n+(j)) - lgamma(n) - lgamma((j)+1) );
  //       }
  //     }
  //     result = n*(1-p)/p - temp;
  //   }
  // }


  return result;
}

////' @export
// [[Rcpp::export]]
SEXP sumNegativeBinomialLfacYObs(SEXP n_, SEXP p_, SEXP nn_, SEXP lower_, SEXP upper_)
{
  NumericVector n(n_);
  NumericVector p(p_);
  NumericVector nn(nn_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector uupper(upper.length());
  uupper = pmin(upper, R::qnbinom(0.9999999999, n(0), p(0), 1, 0));

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    result(j) = sumNegativeBinomialLfacY(n(0), p(0), nn(0), lower(j), uupper(j));
  }

  return result;
}

////' @export
// [[Rcpp::export]]
SEXP sumNegativeBinomialLfacYLat(SEXP n_, SEXP p_, SEXP nn_, SEXP lower_, SEXP upper_)
{
  NumericVector n(n_);
  NumericVector p(p_);
  NumericVector nn(nn_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector uupper(upper.length());
  uupper = pmax(upper+1, R::qnbinom(0.9999999999, n(0), p(0), 1, 0));

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    double temp1 = sumNegativeBinomialLfacY(n(0), p(0), nn(0), 0, lower(j)-1);
    double temp2 = sumNegativeBinomialLfacY(n(0), p(0), nn(0), upper(j)+1, uupper(j));
    result(j) = temp1 + temp2;
  }

  return result;
}
