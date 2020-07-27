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
SEXP XColMinusY(SEXP x, SEXP y) {
  NumericMatrix xx(x);
  NumericVector yy(y);
  NumericMatrix result(xx.nrow(), xx.ncol());

  for(int j=0; j<xx.ncol(); j++){
    result(_,j) = xx(_,j) - yy;
  }
  return result;
}


