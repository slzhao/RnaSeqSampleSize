#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]                                                             
NumericVector cumsumBorder(NumericVector x, double border) {
  NumericVector res(1);
  double acc = 0;
  for (int i=(x.size()-1); i>-1; --i) {
    acc += x[i];
    if (acc > border)  {
        res=x.size()-(i+2);
        return res;
    };
  }
  return 0;
}