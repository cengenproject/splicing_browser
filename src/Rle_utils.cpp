#include <Rcpp.h>
using namespace Rcpp;


// Parallel mean of a list of vectors the same length

// [[Rcpp::export]]
NumericVector cpmean(List x) {
  int nb=x.size();
  NumericVector x0=x[0];
  int n=x0.size();
  
  NumericMatrix mat(nb,n);
  
  //fill in matrix
  for(int i=0;i<nb;i++){
    NumericVector cur_x = x[i];
    
    for(int j=0;j<n;j++){
      mat(i,j) = cur_x[j];
    }
  }
  
  //Apply mean to each column
  NumericVector out(n);
  for(int j=0; j<n; j++){
    out[j] = 0;
    for(int i=0; i<nb; i++){
      out[j] = out[j] + mat(i,j);
    }
    out[j] = out[j] / nb;
  }
  return out;
}



/*** R
timesTwo(42)
*/
