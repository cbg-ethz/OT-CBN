
all_chains <- function(poset) {
  all_chains_nested <- function() {
    first_nodes <- which(apply(poset, 2, function(x) { sum(x[visited == FALSE]) }) == 0)
    first_nodes = first_nodes[which(visited[first_nodes] == FALSE)]
    
    for(f_node in first_nodes) {
      
      visited[f_node] <<- TRUE
      cur_chain[f_node] <<- sum(visited)
      
      if(all(visited)) {
        chain_index <<- chain_index + 1
        chains[[chain_index]] <<- order(cur_chain)
        visited[f_node] <<- FALSE
        return()
      } else {
        all_chains_nested()
      }
      visited[f_node] <<- FALSE
    }
  }
  
  n = ncol(poset)
  visited <- rep(FALSE, n)
  cur_chain <- rep(-1, n)
  chains <- list()
  chain_index <- 0
  
  all_chains_nested()
  return(chains)
}


possible_chains_for_genotype <- function(chains_m, genotype) {
  if(length(genotype) == 0)
    return(chains_m)
  
  number_of_mutations <- length(genotype)
  if(number_of_mutations == 1)
    m = matrix(chains_m[, 1:number_of_mutations], nrow = nrow(chains_m), byrow=TRUE)
  else
    m = chains_m[, 1:number_of_mutations]
  
  as.matrix(chains_m[apply(m, 1,  function(x) {all(genotype %in% x) } ) , ])
}


compute_lambda_exit <- function(poset, lambda, genotype) {
  n = ncol(poset)
  c_geno = setdiff(1:n, genotype)
  
  lambda_exit <- 0
  
  for(mut in c_geno) {
    parents <- which(poset[ ,mut] == 1)
    if(all(parents %in% genotype))
      lambda_exit = lambda_exit + lambda[mut]
  }
  
  return(lambda_exit)
}

time_components <- function( poset, chain, mut)
{ 
  parents <- which(poset[, mut] == 1)
  
  u_list = mut
  index = which(chain == mut)
  if(index > 1) {
    for(m_index in (index-1):1) { 
      
      if(chain[m_index] %in% parents)
        break
      u_list = c(u_list, chain[m_index])
    }
  }
  unlist(u_list)
}

lambda_exit_for_chain <- function(poset, lambda, chain)
{
  n = ncol(poset)
  genotype <- c()
  lambda_exit <- rep(0, n)
  for(i in 1:n) {
    
    mut = unlist(chain)[i]
    lambda_exit[i] = compute_lambda_exit(poset, lambda, genotype)
    genotype = c(genotype, mut)
  }
  lambda_exit
}

################################# 

compute_integral_over_chain <- function(poset, chain, genotype, lambda, S) {
  mutation_nr = ncol(poset)
  
  lambda_exit = lambda_exit_for_chain(poset, lambda, chain)
  
  last_obs_event = length(genotype)
  
  U <- rep(0, mutation_nr)
  for(k in 1:mutation_nr) {
    # E_function computes the integral: integ(u_k f(u)du)
    U[k] = E_function(lambda_exit, S, mutation_nr, k, last_obs_event)
  }
  U0 = E_function(lambda_exit, S, mutation_nr, 0, last_obs_event)
  
  U = U[order(unlist(chain))]
  
  E <- rep(0, mutation_nr+1)
  for(k in 1:mutation_nr) {
    u_list = time_components(poset, chain, k)
    E[k] = sum(U[u_list])
  }
  
  E[mutation_nr+1] = U0
  
  return(E)
}


integral_over_chain_v2 <- function(poset, chain, genotype, lambda, S, usePrec, prec) {
  mutation_nr = ncol(poset)
  
  lambda_exit = lambda_exit_for_chain(poset, lambda, chain)
  
  last_obs_event = length(genotype)
  if(usePrec == TRUE) {
    tmpE = fast_E_prec(mutation_nr, last_obs_event, S, lambda_exit, prec)
  } else {
    tmpE = exp(fast_E(mutation_nr, last_obs_event, S, lambda_exit)$finalE)
  }
  
  U = tmpE[2:(mutation_nr+1)]
  U0 = tmpE[1]
  
  U = U[order(unlist(chain))]
  
  #E <- tmpE[mutation_nr+1, , 1]  # initialize the 
  E <- tmpE  # initialize the 
  for(k in 1:mutation_nr) {
    u_list = time_components(poset, chain, k)
    E[k] = sum(U[u_list])
  }
  
  E[mutation_nr+1] = U0
  
  return(E)
}

expected_time_for_geno_v2 <- function(poset, lambda, genotype, S, all_possible_chains, usePrec, prec) {
  n = ncol(poset)
  E <- rep(0, n+1)
  
  geno_chains <- possible_chains_for_genotype(all_possible_chains, genotype)
  
   for(row in 1:nrow(geno_chains) ) {
     chain = geno_chains[row, ]
     E = E + integral_over_chain_v2(poset, chain, genotype, lambda, S, usePrec, prec)
     #print(row)
     #print(length(integral_over_chain_v2(poset, chain, genotype, lambda, S, usePrec, prec)))
     #print(length(E))
     #print("---")
   }
  
#   E2 = t(apply(geno_chains, 1, function(chain) {
#     integral_over_chain_v2(poset, chain, genotype, lambda, S, usePrec, prec)
#   }) )
  
#   print("**********")
#   print(dim(E2))
   
#     E = apply(E2, 2, sum)
#   print("**********")
#   print(E3)
#   
#    print("**********")
#    print(E3-E)
  
  # E[n+1] is equal to P
  E[1:n] = E[1:n]/E[n+1]
  
  if(any(E[1:n] < 0.0))
  {
	  print( paste("Error, poset size=", ncol(poset) ) )
  }
  return(E[1:n])
}


m_step <- function(genotype_list, S_list, old_lambda, poset, chains_m, usePrec, prec) {
  geno_nr = length(genotype_list)
  
  mutation_nr = ncol(poset)
  
  
  if(usePrec == FALSE) {
    expected_obs_time = matrix(0, geno_nr, mutation_nr)  
  } else {
    # 10 is not the ultimate precision. The precision of inner function will be written over this.
    expected_obs_time <- array(mpfr(0, prec=10), dim=c(geno_nr, mutation_nr))
    #old_lambda = mpfr(old_lambda, precBits=prec)
  }
  
#   print("m-step-begin")
#   print(Sys.time()-t1)
#   
  
  expected_obs_time <- foreach(i=1:geno_nr,.combine=rbind) %dopar%{
    expected_time_for_geno_v2(poset, old_lambda, genotype_list[[i]], S_list[i], chains_m, usePrec, prec)  
  }
  
  
#   expected_obs_time = t(sapply(1:geno_nr, function(i) { expected_time_for_geno_v2(poset, old_lambda, genotype_list[[i]], S_list[i], chains_m, usePrec, prec)}))
# 
#   for(i in 1:geno_nr) {
#   expected_obs_time[i, ] = expected_time_for_geno_v2(poset, old_lambda, genotype_list[[i]], S_list[i], chains_m, usePrec, prec)  
#   }
  
#   print("m-step-end")
#   print(Sys.time()-t1)
  
  #print(expected_obs_time)
  #print(apply(expected_obs_time, 2, mean))
  
  new_lambda = 1 / apply(expected_obs_time, 2, mean)
  new_lambda
}

do_em <- function(poset, genotype_list, S_list, maxIter, usePrec, prec, verbose) {
  mutation_nr = ncol(poset)
  
  old_lambda = lambda = rep(1/mean(S_list), mutation_nr)
    #runif(mutation_nr, 0.1, 0.8)
  
  if(usePrec == TRUE) {
    lambda = mpfr(lambda, precBits=prec)
    old_lambda = mpfr(old_lambda, precBits=prec)
  }
  
  chains <- all_chains(poset)
  chains_m <- data.frame(matrix(unlist(chains), nrow=length(chains), byrow=T))
  
  
  iter = 1
  while( ( (iter < 3) || (sum(abs(lambda-old_lambda)) > 10^(-3) * sum(abs(lambda))) ) && (iter < maxIter) ) {
    if(verbose) {
      print(paste("EM step ", iter, ":     ", paste(lambda, collapse=' ') , sep='') )
    }
    old_lambda = lambda
    lambda = m_step(genotype_list, S_list, old_lambda, poset, chains_m, usePrec, prec)
    iter = iter + 1
  }
  lambda
}


otcbn_on_samples <- function(poset, samples, maxIter = 100, usePrec=FALSE, precision=NA) {
  genotype_list <- samples$genotype_list
  S_list = samples$T_sampling
  
  do_em(poset, genotype_list, S_list, maxIter, usePrec, precision)
}

#' otcbn
#' @export
otcbn <- function(poset, genotype_list, S_list, maxIter = 100, verbose=FALSE) {
  if(verbose) {
    print(paste("Start EM on ", Sys.time()) )
  }
  do_em(poset, genotype_list, S_list, maxIter, FALSE, NA, verbose)
}


#' otcbn_prec
#' @export
otcbn_prec <- function(poset, genotype_list, S_list, maxIter = 100, prec=100, verbose=FALSE) {
  if(verbose) {
    print(paste("Start EM on ", Sys.time()) )
  }
  do_em(poset, genotype_list, S_list, maxIter, TRUE, prec, verbose)
}


#' probability_of_compatible_genoptype
#' @export
probability_of_compatible_genoptype <- function(poset, lambda, genotype, S) {
  if(S <= 0) {
    if(length(genotype) == 0) {
      return(1.0)
    } else {
      return(0.0)  
    }
  }
  n = ncol(poset)
  E <- rep(0, n+1)
  
  chains <- all_chains(poset)
  all_possible_chains <- data.frame(matrix(unlist(chains), nrow=length(chains), byrow=T))
  

  geno_chains <- possible_chains_for_genotype(all_possible_chains, genotype)
  
#   
#    for(row in 1:nrow(geno_chains) ) {
#      chain = geno_chains[row, ]
#      print(row)
#      E = E + integral_over_chain_v2(poset, chain, genotype, lambda, S, usePrec, prec)
#      print(integral_over_chain_v2(poset, chain, genotype, lambda, S, usePrec, prec))
#    }
  usePrec = FALSE
  prec = NA
  E2 = t(apply(geno_chains, 1, function(chain) {
    integral_over_chain_v2(poset, chain, genotype, lambda, S, usePrec, prec)
  }) )
  

   
  E = apply(E2, 2, sum)

  
  # E[n+1] is equal to P
  P = E[n+1] * prod(lambda)
  return(P)
}


#' get_compatible_observations
#' @export
get_compatible_observations <- function(poset, obs_events) {
  
  chains <- all_chains(poset)
  chains_m <- data.frame(matrix(unlist(chains), nrow=length(chains), byrow=T))
  
  genotype_list = apply(obs_events, 1, function(x) {
    which(x == 1)
  })
  
  indexes = c()
  i = 1
  for(genotype in genotype_list) {    
    if(nrow(possible_chains_for_genotype(chains_m, genotype) ) > 0 )
      indexes = c(indexes, i)
    i = i + 1 
  }
  
  indexes
}
