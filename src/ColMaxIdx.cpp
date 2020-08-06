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

////' @export
// [[Rcpp::export]]
SEXP ColMaxIdx(SEXP w) {
  NumericMatrix ww(w);
  NumericVector result(ww.nrow());

  for(int i=0; i<ww.nrow(); i++){
    result(i) = 1;
    double mx = ww(i,0);
    for(int j=1; j<ww.ncol(); j++){
      if(ww(i,j)>mx){
        result(i) = j+1;
        mx = ww(i,j);
      }
    }
  }

  return result;
}
