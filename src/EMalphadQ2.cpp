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

// // [[Rcpp::export]]
// SEXP EMalphadQ2left(SEXP x, SEXP z, SEXP p, SEXP q) {
//  NumericMatrix xx(x);
//  NumericVector zz(z);
//  NumericVector pp(p);
//  NumericVector qq(q);
//  NumericMatrix result(xx.nrow(), xx.ncol());
//
//  for(int j=0; j<xx.ncol(); j++){
//    result(_,j) = xx(_,j) * (zz*pp*qq);
//  }
//  return result;
//
// }

// [[Rcpp::export]]
SEXP EMalphadQ2(SEXP x, SEXP z, SEXP p, SEXP q) {
  NumericMatrix xx(x);
  NumericVector zz(z);
  NumericVector pp(p);
  NumericVector qq(q);
  NumericMatrix temp(xx.nrow(), xx.ncol());
  NumericMatrix result(xx.ncol(), xx.ncol());

  // for(int j=0; j<xx.ncol(); j++){
  //   temp(_,j) = xx(_,j) * (zz*pp*qq);
  // }

  for(int i=0; i<result.nrow(); i++){
    temp(_,i) = xx(_,i) * (zz*pp*qq);
    for(int j=0; j<result.ncol(); j++){
      result(i,j) = -1.0* sum(temp(_,i)*xx(_,j));
    }
  }

  return result;
}
