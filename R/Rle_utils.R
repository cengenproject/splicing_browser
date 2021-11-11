# Start ----



# For tests
test_list_of_SimpleRleList <- list(
  RleList(c(I = Rle(c(1,1,1,1,8,8,8,2,2,2,2)),
            II = Rle(c(2,2,2,5,5,5,2,2,2,2,3))), compress = FALSE),
  RleList(c(I = Rle(c(2,2,2,5,5,5,2,2,2,2,3)),
            II = Rle(c(1,1,1,1,8,8,8,2,2,2,2))), compress = FALSE),
  RleList(c(I = Rle(c(3,3,2,4,4,5,2,6,6,6,6)),
            II = Rle(c(1,5,8,8,8,8,8,8,2,2,2))), compress = FALSE),
  RleList(c(I = Rle(c(12,12,12,12,12,12,3,3,3,3,3)),
            II = Rle(c(9,9,9,5,5,5,5,5,2,2,2))), compress = FALSE)
)

# direct_quantile <- function(x, i) sort(x)[1+floor(length(x)*i)]
direct_rank <- function(x, i) sort(x)[i]


# Expect a list of SimpleRleList, where each has one Rle per chromosome.
# Function .f must accept a vector of integer values and return a single int.
# Note this is very slow and is only to use for testing small lists.
apply_per_position <- function(list_of_rlelist, .f, ...){
  # check number of chromosomes
  nb_chrom <- purrr::map_int(list_of_rlelist, length)
  nb_chrom <- unique(nb_chrom)
  if(length(nb_chrom) != 1) stop("Not all RleLists have the same number of chromosomes")
  
  out <- vector("list", length = nb_chrom)
  for(cur_chrom in seq_len(nb_chrom)){
    expanded_vecs <- purrr::map(list_of_rlelist,
                         ~ rep(.x[[cur_chrom]]@values, .x[[cur_chrom]]@lengths))
    
    chrom_length <- purrr::map_int(expanded_vecs, length)
    chrom_length <- unique(chrom_length)
    if(length(chrom_length) != 1) stop("The chromosome lengths differ.")
    
    
    chrom_as_mat <- matrix(unlist(expanded_vecs),
                           byrow = TRUE,
                           nrow = length(list_of_rlelist),
                           ncol = chrom_length)
    out[[cur_chrom]] <- Rle(apply(chrom_as_mat, 2, .f, ...))
  }
  out <- RleList(out)
  names(out) <- names(list_of_rlelist[[1]])
  out
}


Rcpp::sourceCpp("src/Rle_utils.cpp", echo = FALSE)




# Note: pmean is as or more efficient as R implementation,
# but good for tests

# # pmean: parallel mean of a group of Rle.
# # The main function is implemented in C++
# 
# pmean <- function(...) cpmean(list(...))
# 
# 
# # Support functions to call on every chromosome in a simpleRleList
# # Copied from S4Vectors (v 0.30.0)
# 
# setMethod("pmean", "Rle", function(...)
#   S4Vectors:::.psummary.Rle(pmean, ...))
# 
# 
# setMethod("pmean", "RleList", function(...)
#   mendoapply(pmean, ...))




pmean <- function(...) {
  ll <- list(...)
  Reduce(`+`, ll)/length(ll)
}

# Support functions to call on every chromosome in a simpleRleList
# Copied from S4Vectors (v 0.30.0)

setMethod("pmean", "Rle", function(...)
  S4Vectors:::.psummary.Rle(pmean, ...))


setMethod("pmean", "RleList", function(...)
  mendoapply(pmean, ...))

# Test!
stopifnot(all.equal(unlist(do.call(pmean,test_list_of_SimpleRleList)),
                    suppressWarnings(unlist(apply_per_position(test_list_of_SimpleRleList, mean)))))




# pmedian ----
pmedian <- function(...) cpmedian(list(...))
setMethod("pmedian", "Rle", function(...)
  S4Vectors:::.psummary.Rle(pmedian, ...))


setMethod("pmedian", "RleList", function(...)
  mendoapply(pmedian, ...))



# Test!
stopifnot(all.equal(unlist(do.call(pmedian,test_list_of_SimpleRleList)),
                    suppressWarnings(unlist(apply_per_position(test_list_of_SimpleRleList, median)))))




# p lower percentile ----
plower <- function(...) c_low_percentile(list(...))
setMethod("plower", "Rle", function(...)
  S4Vectors:::.psummary.Rle(plower, ...))


setMethod("plower", "RleList", function(...)
  mendoapply(plower, ...))




# Test!
stopifnot(all.equal(unlist(do.call(plower,
                                   test_list_of_SimpleRleList)),
                    suppressWarnings(unlist(apply_per_position(test_list_of_SimpleRleList,
                                                               direct_rank, 3)))))





# p higher percentile ----
phigher <- function(...) c_high_percentile(list(...))
setMethod("phigher", "Rle", function(...)
  S4Vectors:::.psummary.Rle(phigher, ...))


setMethod("phigher", "RleList", function(...)
  mendoapply(phigher, ...))




# Test!
stopifnot(all.equal(unlist(do.call(phigher,
                                   test_list_of_SimpleRleList)),
                    suppressWarnings(unlist(apply_per_position(test_list_of_SimpleRleList,
                                                               direct_rank, length(test_list_of_SimpleRleList)-3)))))

# bench::mark(do.call(phigher,
#                     test_list_of_SimpleRleList),
# apply_per_position(test_list_of_SimpleRleList,
#                    direct_quantile, .9),
# check = FALSE)


rm(test_list)
rm(test_list2)
rm(test_list_of_SimpleRleList)