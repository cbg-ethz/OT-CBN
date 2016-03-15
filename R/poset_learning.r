


#' loglikelihood_cbn
#' @export
loglikelihood_cbn <- function(poset, lambdas, sampling_times, obs_events){
  p = ncol(poset)
  lattice_size = orderIdeals(poset)$lattice_size
  
  compatible_geno = compatible_genotypes(obs_events, poset)
  nr_compatible_genotypes = length(compatible_geno$compatible_indexes)
  nr_incompatible_genotypes = nrow(obs_events) - length(compatible_geno$compatible_indexes)
  
  # fraction of the data that are compatible with the given poset
  alpha = nr_compatible_genotypes /nrow(obs_events)
  
  q_e = 1/(2^p - lattice_size  )
  
  incompatible_ll = 0.0
  if(nr_incompatible_genotypes > 0)
    incompatible_ll = nr_incompatible_genotypes* (log(1-alpha) + log(q_e))
  
  
  
  compatible_ll = 0.0
  if( nr_compatible_genotypes > 0 )
  {
    genotypes = obs_events[compatible_geno$compatible_indexes,]
    sampling_times = sampling_times[compatible_geno$compatible_indexes]
    
    compatible_ll = 0.0
    for(i in 1:nr_compatible_genotypes){
      p = probability_of_compatible_genoptype(poset, lambdas, which(genotypes[i,] == 1), sampling_times[i] )
      #     if(p == 0)
      #       p = 10^-9
      compatible_ll = compatible_ll + log(p)
    }
    
    compatible_ll = compatible_ll + nr_compatible_genotypes * log(alpha)
  }
  
  list(ll=incompatible_ll + compatible_ll, alpha=alpha)
  
}


#' violation_freqs
#' @export
violation_freqs <- function(obs_events) {
  p = ncol(obs_events)
  N = nrow(obs_events)
  violations <- matrix(0, p, p)
  for(i in 1:p) {
    for(j in 1:p) {
      violations[i, j] = sum(obs_events[, i] == 0 & obs_events[, j] == 1)
    }
  }
  diag(violations) = N+1
  violations/N
}


#' maximal_poset
#' @export
maximal_poset <- function(violations, eps) {
  p = ncol(violations)
  poset = matrix(0, p, p)
  for(i in 1:p) {
    for(j in 1:p) {
      if(violations[i, j] <= eps & poset[j, i] == 0) {

        parents_i = which(poset[, i] == 1)

            # check if adding (i,j) introduces a cycle                          
        if(all(poset[j, parents_i] == 0) ) 
        {
          poset[parents_i, j] = 1
          poset[i, j] = 1  
          
        } 
      }
    }
  }

  diag(poset) = 0
  poset
}
