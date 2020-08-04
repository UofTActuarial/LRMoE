// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <Rcpp.h>
#include <RcppNumerical.h>
using namespace Rcpp;
using namespace Numer;

// For accessing constant PI
#define _USE_MATH_DEFINES
#include <math.h>

class BurrLogY: public Func
{
private:
  double k;
  double c;
  double lambda;
public:
  BurrLogY(double k_, double c_, double lambda_) : k(k_), c(c_), lambda(lambda_) {}

  double operator()(const double& u) const
  {
    double temp;
    if(isinf(u)){
      temp = 0.0;
    }else{
      temp = u* exp( log(c*k/lambda) + (c-1)*log(exp(u)/lambda) + (-k-1)*log(1+pow(exp(u)/lambda, c)) + u );
        // log(u) * c*k/lambda * pow((u/lambda), c-1) * pow(1+pow(u/lambda, c), -k-1);
    }
    if(isnan(temp) || isinf(temp)){
      return 0.0;
    }else{
      return temp;
    }
  }
};

Rcpp::List intBurrLogY(double k, double c, double lambda, double lower, double upper)
{
  BurrLogY f(k, c, lambda);
  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code,
                               250, 1e-8, 1e-8, Integrator<double>::GaussKronrod71);
  return Rcpp::List::create(
    // Rcpp::Named("true") = true_val,
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}

////' @export
// [[Rcpp::export]]
SEXP intBurrLogYObs(SEXP k_, SEXP c_, SEXP lambda_, SEXP lower_, SEXP upper_)
{
  NumericVector k(k_);
  NumericVector c(c_);
  NumericVector lambda(lambda_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    result(j) = intBurrLogY(k(0), c(0), lambda(0), lower(j), upper(j))[0];
  }

  return result;
}

////' @export
// [[Rcpp::export]]
SEXP intBurrLogYLat(SEXP k_, SEXP c_, SEXP lambda_, SEXP lower_, SEXP upper_)
{
  NumericVector k(k_);
  NumericVector c(c_);
  NumericVector lambda(lambda_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    double temp1 = intBurrLogY(k(0), c(0), lambda(0), R_NegInf, lower(j))[0];
    double temp2 = intBurrLogY(k(0), c(0), lambda(0), upper(j), R_PosInf)[0];
    result(j) = temp1 + temp2;

  }

  return result;
}

/////////////////////////////////////////////////////////////////////////////////////////

class BurrPolY: public Func
{
private:
  double k;
  double c;
  double lambda;
  double cc;
  double ll;
public:
  BurrPolY(double k_, double c_, double lambda_, double cc_, double ll_) : k(k_), c(c_), lambda(lambda_), cc(cc_), ll(ll_) {}

  double operator()(const double& u) const
  {
    double temp;
    if(isinf(u)){
      temp = 0.0;
    }else{
      temp = log1p( exp(cc* (u-log(ll)) ) ) * exp( log(c*k/lambda) + (c-1)*log(exp(u)/lambda) + (-k-1)*log(1+pow(exp(u)/lambda, c)) + u );
        // log(1+pow(exp(u)/ll, cc))* exp( log(c*k/lambda) + (c-1)*log(exp(u)/lambda) + (-k-1)*log(1+pow(exp(u)/lambda, c)) + u );
      // log(u) * c*k/lambda * pow((u/lambda), c-1) * pow(1+pow(u/lambda, c), -k-1);
    }
    if(isnan(temp) || isinf(temp)){
      return 0.0;
    }else{
      return temp;
    }
  }
};

Rcpp::List intBurrPolY(double k, double c, double lambda, double cc, double ll, double lower, double upper)
{
  BurrPolY f(k, c, lambda, cc, ll);
  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code,
                               250, 1e-8, 1e-8, Integrator<double>::GaussKronrod71);
  return Rcpp::List::create(
    // Rcpp::Named("true") = true_val,
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}

////' @export
// [[Rcpp::export]]
SEXP intBurrPolYObs(SEXP k_, SEXP c_, SEXP lambda_, SEXP cc_, SEXP ll_, SEXP lower_, SEXP upper_)
{
  NumericVector k(k_);
  NumericVector c(c_);
  NumericVector lambda(lambda_);
  NumericVector cc(cc_);
  NumericVector ll(ll_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    result(j) = intBurrPolY(k(0), c(0), lambda(0), cc(0), ll(0), lower(j), upper(j))[0];
  }

  return result;
}

////' @export
// [[Rcpp::export]]
SEXP intBurrPolYLat(SEXP k_, SEXP c_, SEXP lambda_, SEXP cc_, SEXP ll_, SEXP lower_, SEXP upper_)
{
  NumericVector k(k_);
  NumericVector c(c_);
  NumericVector lambda(lambda_);
  NumericVector cc(cc_);
  NumericVector ll(ll_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    double temp1 = intBurrPolY(k(0), c(0), lambda(0), cc(0), ll(0), R_NegInf, lower(j))[0];
    double temp2 = intBurrPolY(k(0), c(0), lambda(0), cc(0), ll(0), upper(j), R_PosInf)[0];
    result(j) = temp1 + temp2;

  }

  return result;
}

