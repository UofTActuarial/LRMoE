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
SEXP XPlusYZ(SEXP x, SEXP y, SEXP z) {
  NumericVector xx(x);
  NumericVector yy(y);
  NumericVector zz(z);
  NumericVector result(xx.length());

  // for(int j=0; j<xx.ncol(); j++){
  //   result(_,j) = xx(_,j) + yy(_,j)*zz;
  // }

  result = xx + yy*zz;
  return result;
}

