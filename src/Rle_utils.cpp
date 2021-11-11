#include <Rcpp.h>
using namespace Rcpp;


/*** R
# For basic testing
test_list <- list(1:5, 2:6, 7:3, 8:4, 6:10)
stopifnot(length(unique(sapply(test_list, length))) == 1)

p_apply <- function(my_list, .f, ...){
  stopifnot(dplyr::n_distinct(sapply(my_list, length)) == 1)
  
  mat <- matrix(unlist(my_list), byrow = TRUE, nrow = length(my_list))
  
  apply(mat, 2, .f, ...)
}
*/




// Parallel mean of a list of vectors the same length (for tests)
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
stopifnot(all.equal(cpmean(test_list),
                    Reduce(`+`, test_list)/length(test_list)))
stopifnot(all.equal(cpmean(test_list),
                    p_apply(test_list, mean)))
*/






// Implement parallel median
// [[Rcpp::export]]
NumericVector cpmedian(List x) {
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
  
  //Apply median to each column
  NumericVector out(n);
  for(int j=0; j<n; j++){
    NumericVector cur_col = mat(_,j);
    out[j] = median(cur_col);
  }
  return out;
}

/*** R
stopifnot(all.equal(cpmedian(test_list),
                    p_apply(test_list, median)))
*/



// Implement parallel percentiles

// first, percentile for 1 vector
double percentile(NumericVector x, int nth){
  
  std::nth_element (x.begin(), x.begin()+nth-1, x.end());
  
  return x[nth-1];
}

// Use it for list, first with median to double check
// [[Rcpp::export]]
NumericVector cpmedian2(List x) {
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
  
  //Apply median to each column
  
  int nth = 1+nb/2;
  NumericVector out(n);
  for(int j=0; j<n; j++){
    NumericVector cur_col = mat(_,j);
    out[j] = percentile(cur_col, nth);
  }
  return out;
}

/*** R
stopifnot(all.equal(cpmedian(test_list),
                    cpmedian2(test_list)))
test_list2 <- list(c( 1, 1, 2, 2, 2),
                   c(12,12, 1, 7, 8),
                   c( 2, 2, 3, 3, 3),
                   c(10,10,12,12,12),
                   c(10,11,12,13,14))

stopifnot(all.equal(cpmedian(test_list2),
                    cpmedian2(test_list2)))
*/


// then lower value
// [[Rcpp::export]]
NumericVector c_low_percentile(List x) {
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
  
  //not actual percentile (as peak only in 1 neur type
  // would have ~3 samples, out of 140 that gets missed)
  
  int nth = 3;
  NumericVector out(n);
  for(int j=0; j<n; j++){
    NumericVector cur_col = mat(_,j);
    out[j] = percentile(cur_col, nth);
  }
  return out;
}

/*** R
# direct_quantile <- function(x, i) sort(x)[1+floor(length(x)*i)]
direct_rank <- function(x, i) sort(x)[i]
stopifnot(all.equal(c_low_percentile(test_list),
                    p_apply(test_list, direct_rank, 3)))
test_list2 <- list(c( 1, 1, 2, 2, 2),
                   c(12,12, 1, 7, 8),
                   c( 2, 2, 3, 3, 3),
                   c(10,10,12,12,12),
                   c(10,11,12,13,14))

stopifnot(all.equal(c_low_percentile(test_list2),
                    p_apply(test_list2, direct_rank, 3)))
*/




// and 90th percentile
// [[Rcpp::export]]
NumericVector c_high_percentile(List x) {
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
  
  //not actual percentile (as peak only in 1 neur type
  // would have ~3 samples, out of 140 that gets missed)
  
  int nth = nb-3;
  NumericVector out(n);
  for(int j=0; j<n; j++){
    NumericVector cur_col = mat(_,j);
    out[j] = percentile(cur_col, nth);
  }
  return out;
}

/*** R
stopifnot(all.equal(c_high_percentile(test_list),
                    p_apply(test_list, direct_rank, length(test_list) -3)))
test_list2 <- list(c( 1, 1, 2, 2, 2),
                   c(12,12, 1, 7, 8),
                   c( 2, 2, 3, 3, 3),
                   c(10,10,12,12,12),
                   c(10,11,12,13,14))

stopifnot(all.equal(c_high_percentile(test_list2),
                    p_apply(test_list2, direct_rank, length(test_list2) -3)))
*/



