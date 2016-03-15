
fast_E_prec <- function(n, l, S, lambdas, prec = 100, verbose=FALSE) {
  
  # k = 1 : u_k = 1 for prob. K > 1 means for u_{k-1}
  
  # E(p, k, q)   :  p is the number of total events (so far). k : the place of u_k, q: lambda vector is substracted with lambda_q
  # k = n + 1 means probability (u_{n+1} = 1)
  l = l + 1
  n = n + 1
  
  E <- array(mpfr(1.0, precBits=prec), dim=c(n, n, n))
  #E[1, , ] = 1.0
  
  #lambdas = mpfr(c(0, lambdas), precBits=prec)
  lambdas = c(mpfr(c(0), precBits=prec),  lambdas)
  
  for(p in   from_to_values(2, min(l+1, n))) {
    S_pl = ifelse(p > l, 1, -1)
    I_pl = ifelse(p <= l, 1, 0)
    
    #for(q in c(1, from_to_values(p+1, n)) ){
    q_indexes = c(1, from_to_values(p+1, n))
    lambda_vector = lambdas[p]-lambdas[q_indexes]
    zero_q_indexes = q_indexes[which(abs(lambda_vector) <= exp(-30))]
    non_zero_q_indexes = q_indexes[which(abs(lambda_vector) > exp(-30))]
    non_zero_lambdas = lambda_vector[which(abs(lambda_vector) > exp(-30))]
    
    tmpExp_vector = exp(-non_zero_lambdas * S)  / non_zero_lambdas
    k_indexes <- from_to_values(1, p-1)
    
    lambdas_inv_matrix = t(new("mpfrMatrix", rep(I_pl/non_zero_lambdas, length(k_indexes)), Dim=c(length(non_zero_q_indexes), length(k_indexes))))
                           
    tmpExp_matrix = t(new("mpfrMatrix", rep(S_pl  * tmpExp_vector, length(k_indexes)), Dim = c(length(non_zero_q_indexes), length(k_indexes))))
    
    E_result_for_zero = array(mpfr(c((S^(p-1) )/factorial(p-1), rep((S^p )/factorial(p), p-1)), precBits=prec), dim=c(length(k_indexes)+1, length(zero_q_indexes)))
    E[p, c(k_indexes, p), zero_q_indexes] = E_result_for_zero
    
    E[p, k_indexes, non_zero_q_indexes] = (lambdas_inv_matrix) * E[p-1, k_indexes, non_zero_q_indexes] +
    tmpExp_matrix * new("mpfrMatrix", rep(E[p-1, k_indexes, p], length(non_zero_q_indexes)), Dim=c(length(k_indexes), length(non_zero_q_indexes)))
    
    
    k = p
    tsum = mpfr(0.0, precBits=prec)
    if(p > 2)
    tsum = sum(E[p-1, from_to_values(2, p-1), p])
    
    E[p, k, non_zero_q_indexes] = (I_pl * E[p-1, 1, non_zero_q_indexes] / (non_zero_lambdas^2) ) + 
      (S_pl  * tmpExp_vector) * ( (S + 1/non_zero_lambdas)* E[p-1, 1, p] - tsum)
    #       }
    
    #}
  }
  
  # The term q is always one when p > l + 2
  q = 1
  
  
  for( p in  from_to_values(l+2, n) )
  {
    lambda_last = lambdas[p]-lambdas[q]
    for(k in from_to_values(1, p-1)) {
      E[p, k, q] = E[p-1, k, q] / lambda_last
    }
    
    k = p
    E[p, k, q] = E[p-1, 1, q] / (lambda_last^2)
  }
  if(verbose) {
    print(E)  
  }
  

  if(any(E) < 0 )
    print("SSS")
  return(list(finalE=E[n, , 1], totalE=E) )
}


