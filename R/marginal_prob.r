
 MargP_function <- function(lambdas, S, n, l) {
    
    if(n == l) {
      orders = order(abs(lambdas))
      
      if(all(orders == c(n:1)) || all(orders == c(1:n)) ) {
      } else{
        #print("Error, lambda vales should be sorted")
        #print(orders)
      }

      lambdas = lambdas[orders]
    }
    
    lambda_n = lambdas[n]
    #   print("----")
    #   print(lambdas)
    if(max(abs(lambdas)) < 10^(-16) & n == l) { # thr
      return( (S^n)/factorial(n))
    }    
  lambda_n = lambdas[n]
  
  ExpSL = exp(-lambda_n*S)
  
  if( n == 1) {
    tmp = ExpSL/lambda_n
    return( (1/lambda_n) - tmp )
  }
  
  
  if(n > l) {
    return(MargP_function(lambdas[1:(n-1)], S, n-1, l)/lambda_n)
  }
  
  # from here n = l
  Tmp = ExpSL * MargP_function(lambdas[1:(n-1)]-lambda_n, S, n-1, n-1)/lambda_n
  
  return( (MargP_function(lambdas[1:(n-1)], S, n-1, n-1)/lambda_n) - Tmp )
  
}

marg_integral_over_chain <- function(poset, chain, lambda, S) {
  mutation_nr = ncol(poset)
  
  lambda_exit = lambda_exit_for_chain(poset, lambda, chain)
  

  U <- rep(0, mutation_nr)
  
  for(k in 1:mutation_nr) {
    U[k] = MargP_function(lambda_exit, S, mutation_nr, k)
  }
  
   U = U[order(unlist(chain))]
#   E <- rep(0, mutation_nr)
#   for(k in 1:mutation_nr) {
#     u_list = time_components(poset, chain, k)
#     E[k] = sum(U[u_list])
#   }
#   
  return(U)
}

#' marginal_probability_for_mutations
#' @export
marginal_probability_for_mutations <- function(poset, lambda , s) {
  chains <- all_chains(poset)
  chains_m <- data.frame(matrix(unlist(chains), nrow=length(chains), byrow=T))
  
  
  n = ncol(poset)
  Marg_P <- rep(0, n)
  
  for(row in 1:nrow(chains_m) ) {
    chain = chains_m[row, ]
    Marg_P = Marg_P + marg_integral_over_chain(poset, chain, lambda, s)
  }
  
  return(prod(lambda) * Marg_P)
  
}

