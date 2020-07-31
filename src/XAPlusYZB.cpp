#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
SEXP XAPlusYZB(SEXP x, SEXP a, SEXP y, SEXP z, SEXP b) {
  NumericVector xx(x);
  NumericVector yy(y);
  NumericVector zz(z);
  NumericVector aa(a);
  NumericVector bb(b);
  NumericVector result(xx.length());

  result = xx*aa + yy*zz*bb;
  return result;
}
