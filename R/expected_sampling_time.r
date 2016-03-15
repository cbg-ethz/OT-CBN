#' average_sampling_time_for_genotypes
#' @export
average_sampling_time_for_genotypes <- function(genotype_list, unique_genotypes, S_list) {
  unique_geno_str = unlist(lapply(unique_genotypes, function(x){ paste(x, collapse='') }))
  geno_str = unlist(lapply(genotype_list, function(x){ paste(x, collapse='') }))

  unique_geno_str[unique_geno_str == ""] = "0"
  geno_str[geno_str==""] = "0"
  
  sapply(unique_geno_str, function(gen) { 
    mean(S_list[which(geno_str == gen)] )})
}


#' expected_sampling_time_for_genotype_list
#' @export
expected_sampling_time_for_genotype_list <- function(poset, lambdas, training_S_list, test_genotype_list, nr_sample_points) {
  f_est = density(training_S_list, n = nr_sample_points)
  S = f_est$x
  PS = f_est$y

  S.est = foreach( i=1:length(test_genotype_list), .combine=c ) %dopar% { 
    genotype = test_genotype_list[[i]]
    expected_sampling_time_for_genotype(poset, lambdas, genotype, S, PS)
  }
  S.est
}

#' expected_sampling_time_for_genotype
#' @export
expected_sampling_time_for_genotype <- function(poset, lambdas, genotype, S, PS) {
  P_g_S =   sapply(S, function(x) {probability_of_compatible_genoptype(poset, lambdas, genotype, x) })
  my_integrate(S, S * PS * P_g_S) / my_integrate(S, PS * P_g_S)  
}

# equi-distant x
my_integrate <- function(x, y)  {
  #print(y)
  n = length(x)-1
  h <- (max(x) - min(x))/n
  s <- h * (y[1]/2 + sum(y[2:n]) + y[n+1]/2)
  s
}

  
