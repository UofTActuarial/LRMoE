// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <Rcpp.h>
#include <RcppNumerical.h>
using namespace Rcpp;
using namespace Numer;

class GammaLogY: public Func
{
private:
  double m;
  double theta;
public:
  GammaLogY(double m_, double theta_) : m(m_), theta(theta_) {}

  double operator()(const double& u) const
  {
    return u * exp(- m*log(theta) - lgamma(m) + (m-1)*u - exp(u)/theta + u );
  }
};

Rcpp::List intGammaLogY(double m, double theta, double lower, double upper)
{
  // const double a = 3, b = 10;
  // const double lower = 0.3, upper = 0.8;
  // const double true_val = R::pbeta(upper, a, b, 1, 0) -
  //   R::pbeta(lower, a, b, 1, 0);

  // NumericVector m(m_);
  // NumericVector theta(theta_);
  // NumericVector lower(lower_);
  // NumericVector upper(upper_);

  // NumericVector result(lower.length());

  // double mm = m(0);
  // double tt = theta(0);
  // double ll = lower(0);
  // double uu = upper(0);


  GammaLogY f(m, theta);
  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code);
  return Rcpp::List::create(
    // Rcpp::Named("true") = true_val,
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );

  // for(int j=0; j<lower.length(); j++){
  //   ll = lower(j);
  //   uu = upper(j);
  //   double res = integrate(f, ll, uu, err_est, err_code);
  //   result(j) = res;
  // }
  //
  // return result;
}

// //' @export
// [[Rcpp::export]]
SEXP intGammaLogYObs(SEXP m_, SEXP theta_, SEXP lower_, SEXP upper_)
{
  NumericVector m(m_);
  NumericVector theta(theta_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    // Rcpp::List res = intGammaLogY(m(0), theta(0), lower(j), upper(j));
    // if(lower(j)!=upper(j)){
      result(j) = intGammaLogY(m(0), theta(0), lower(j), upper(j))[0];
    // }
  }

  return result;
}

// //' @export
// [[Rcpp::export]]
SEXP intGammaLogYLat(SEXP m_, SEXP theta_, SEXP lower_, SEXP upper_)
{
  NumericVector m(m_);
  NumericVector theta(theta_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    // Rcpp::List res = intGammaLogY(m(0), theta(0), lower(j), upper(j));
    // if(lower(j)!=R_NegInf){
    //   result(j) = intGammaLogY(m(0), theta(0), R_NegInf, lower(j))[0];
    // }
    // if(upper(j)!=R_PosInf){
    //   result(j) = result(j) + intGammaLogY(m(0), theta(0), upper(j), R_PosInf)[0];
    // }
    double temp1 = intGammaLogY(m(0), theta(0), R_NegInf, lower(j))[0];
    double temp2 = intGammaLogY(m(0), theta(0), upper(j), R_PosInf)[0];
    result(j) = temp1 + temp2;

  }

  return result;
}
