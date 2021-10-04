

# pmean: parallel mean of a group of Rle.
# The main function is implemented in C++
Rcpp::sourceCpp(cpmean)
pmean <- function(...) cpmean(list(...))



# Support functions to call on every chromosome in a simpleRleList
# Copied from S4Vectors (v 0.30.0)

setMethod("pmean", "Rle", function(...)
  S4Vectors:::.psummary.Rle(pmean, ...))


setMethod("pmean", "RleList", function(...)
  mendoapply(pmean, ...))

