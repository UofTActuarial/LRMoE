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
SEXP EMalphadQ(SEXP x, SEXP zj, SEXP z, SEXP p) {
  NumericMatrix xx(x);
  NumericVector zzj(zj);
  NumericVector zz(z);
  NumericVector pp(p);
  NumericVector result(xx.ncol());

  for(int j=0; j<xx.ncol(); j++){
    result(j) = sum(xx(_,j)*(zzj - zz*pp));
  }
  return result;

}
