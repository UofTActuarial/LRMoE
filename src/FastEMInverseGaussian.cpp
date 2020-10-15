// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <Rcpp.h>
#include <RcppNumerical.h>
using namespace Rcpp;
using namespace Numer;

// For accessing constant PI
#define _USE_MATH_DEFINES
// #include <math.h>
#include <cmath>

class InvGaussLogY: public Func
{
private:
  double mu;
  double lambda;
public:
  InvGaussLogY(double mu_, double lambda_) : mu(mu_), lambda(lambda_) {}

  double operator()(const double& u) const
  {
    double temp = u * exp( 0.5*log(lambda/(2*M_PI*exp(3*u))) - (lambda*(exp(u)-mu)*(exp(u)-mu))/(2*mu*mu*exp(u)) + u );
    if(std::isnan(temp) || std::isinf(temp)){
      return 0.0;
    }else{
      return temp;
    }
  }
};

Rcpp::List intInvGaussLogY(double mu, double lambda, double lower, double upper)
{
  InvGaussLogY f(mu, lambda);
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

// //' @export
// [[Rcpp::export]]
SEXP intInvGaussLogYObs(SEXP mu_, SEXP lambda_, SEXP lower_, SEXP upper_)
{
  NumericVector mu(mu_);
  NumericVector lambda(lambda_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    result(j) = intInvGaussLogY(mu(0), lambda(0), lower(j), upper(j))[0];
  }

  return result;
}

// //' @export
// [[Rcpp::export]]
SEXP intInvGaussLogYLat(SEXP mu_, SEXP lambda_, SEXP lower_, SEXP upper_)
{
  NumericVector mu(mu_);
  NumericVector lambda(lambda_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    double temp1 = intInvGaussLogY(mu(0), lambda(0), R_NegInf, lower(j))[0];
    double temp2 = intInvGaussLogY(mu(0), lambda(0), upper(j), R_PosInf)[0];
    result(j) = temp1 + temp2;

  }

  return result;
}


//////////////////////////////////////////////////////////////////////////////////



class InvGaussY: public Func
{
private:
  double mu;
  double lambda;
public:
  InvGaussY(double mu_, double lambda_) : mu(mu_), lambda(lambda_) {}

  double operator()(const double& u) const
  {
    double temp = exp(u) * exp( 0.5*log(lambda/(2*M_PI*exp(3*u))) - (lambda*(exp(u)-mu)*(exp(u)-mu))/(2*mu*mu*exp(u)) + u );
    if(std::isnan(temp) || std::isinf(temp)){
      return 0.0;
    }else{
      return temp;
    }
  }
};

Rcpp::List intInvGaussY(double mu, double lambda, double lower, double upper)
{
  InvGaussY f(mu, lambda);
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

// //' @export
// [[Rcpp::export]]
SEXP intInvGaussYObs(SEXP mu_, SEXP lambda_, SEXP lower_, SEXP upper_)
{
  NumericVector mu(mu_);
  NumericVector lambda(lambda_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    result(j) = intInvGaussY(mu(0), lambda(0), lower(j), upper(j))[0];
  }

  return result;
}

// //' @export
// [[Rcpp::export]]
SEXP intInvGaussYLat(SEXP mu_, SEXP lambda_, SEXP lower_, SEXP upper_)
{
  NumericVector mu(mu_);
  NumericVector lambda(lambda_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    double temp1 = intInvGaussY(mu(0), lambda(0), R_NegInf, lower(j))[0];
    double temp2 = intInvGaussY(mu(0), lambda(0), upper(j), R_PosInf)[0];
    result(j) = temp1 + temp2;

  }

  return result;
}


//////////////////////////////////////////////////////////////////////////////////



class InvGaussInvY: public Func
{
private:
  double mu;
  double lambda;
public:
  InvGaussInvY(double mu_, double lambda_) : mu(mu_), lambda(lambda_) {}

  double operator()(const double& u) const
  {
    double temp = exp(-u) * exp( 0.5*log(lambda/(2*M_PI*exp(3*u))) - (lambda*(exp(u)-mu)*(exp(u)-mu))/(2*mu*mu*exp(u)) + u );
    if(std::isnan(temp) || std::isinf(temp)){
      return 0.0;
    }else{
      return temp;
    }
  }
};

Rcpp::List intInvGaussInvY(double mu, double lambda, double lower, double upper)
{
  InvGaussInvY f(mu, lambda);
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

// //' @export
// [[Rcpp::export]]
SEXP intInvGaussInvYObs(SEXP mu_, SEXP lambda_, SEXP lower_, SEXP upper_)
{
  NumericVector mu(mu_);
  NumericVector lambda(lambda_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    result(j) = intInvGaussInvY(mu(0), lambda(0), lower(j), upper(j))[0];
  }

  return result;
}

// //' @export
// [[Rcpp::export]]
SEXP intInvGaussInvYLat(SEXP mu_, SEXP lambda_, SEXP lower_, SEXP upper_)
{
  NumericVector mu(mu_);
  NumericVector lambda(lambda_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    double temp1 = intInvGaussInvY(mu(0), lambda(0), R_NegInf, lower(j))[0];
    double temp2 = intInvGaussInvY(mu(0), lambda(0), upper(j), R_PosInf)[0];
    result(j) = temp1 + temp2;

  }

  return result;
}
