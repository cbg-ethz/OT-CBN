

logsumexp <- function(logs, signs) {
  indexes = which(signs != 0)
  signs = signs[indexes] 
  logs = logs[indexes] 
  
  m = max(logs)
  S = sum(signs*exp(logs-m))
  logR = m + log(abs(S))
  list(logR=logR, sign=sign(S))
}

myExp <- function(myNum) {
  myNum$sign * exp(myNum$logR)
}


fast_E <- function(n, l, S, lambdas, verbose=FALSE) {
  error = FALSE
  # k = 1 : u_k = 1 for prob. K > 1 means for u_{k-1}
  
  # E(p, k, q)   :  p is the number of total events (so far). k : the place of u_k, q: lambda vector is substracted with lambda_q
  # k = n + 1 means probability (u_{n+1} = 1)
  l = l + 1
  n = n + 1
  
  # This array represents the expected values of mutation occurrence time
  E <- array(NA, dim=c(n, n, n))
  E[1, , ] = 0.0
  
  # According to the prop. 2 we first sort the first 'l' rate parameters.
  lambdas = c(0,  lambdas)
  lambda_for_observed_events = lambdas[from_to_values(2, l)]
  sorted_lambda_for_observed_events = sort(lambda_for_observed_events)
  
  # we keep the original ranking. We need this when we want to report back the expected values for each mutation in the 
  # original order.
  original_orders = c()
  lambda_for_observed_events
  for(tmpLambda in lambda_for_observed_events){
    original_orders = c(original_orders, which(tmpLambda == sorted_lambda_for_observed_events))
  }
  original_orders = unique(original_orders)
  
  lambdas = c( 0, sorted_lambda_for_observed_events, lambdas[from_to_values(l+1, n)] )
  
  # We use a dynamic programming approach to fill the three-dimensional matrix E. The values of  
  # the matrix E for p=1 is zero (in log scale). We can compute the values for E[p, , ] using the 
  # the E[p-1, ,]. Therefore we iterate over p values.
  
  # The recursive formula is given in the prop 2. Following for loop compute the E matrix when p is associated to the observed mutations 
  # and the one immediately after that (l+1 mutations in the chain).
  for(p in   from_to_values(2, min(l+1, n))) {
    S_pl = ifelse(p > l, 1, -1)
    I_pl = ifelse(p <= l, 1, 0)
    
    for(q in c(1, from_to_values(p+1, min(n, l+1))) ){
      lambda_last = lambdas[p] - lambdas[q]
      
      if(abs(lambda_last) < exp(-30)) {
        
        ###### lambda_last zero
        E[p, 1, q] = (p-1) * log(S) - sum(log( 1:(p-1)))
        #(S^(p-1) )/factorial(p-1)
        
        for(k in  from_to_values(2, p)) {
          E[p, k, q] = p * log(S) - sum(log( 1:p))
          # (S^(p))/factorial(p)
        }
      } else  {
        
        
#         tmpExp = exp(-lambda_last*S)/lambda_last
        signTmpExp = sign(lambda_last)
        logTmpExp = -lambda_last*S - log(abs(lambda_last))
        
        
        for(k in  from_to_values(1, p-1)) {
           tmp_log = logsumexp( c(-log(abs(lambda_last) ) + E[p-1, k, q], logTmpExp + E[p-1, k, p] ), 
                                  c(sign(lambda_last)*I_pl, signTmpExp*S_pl))
           
           E[p, k, q] = tmp_log$logR
          #(I_pl / lambda_last) * E[p-1, k, q] +  (S_pl  * tmpExp) * E[p-1, k, p]
          if(tmp_log$sign == -1)
           error = TRUE
        }
        
        k = p
        #  tsum = 0.0
        #  if(p > 2)
        #    tsum = sum(E[p-1, from_to_values(2, p-1), p])
        #print("H1")
        
#            # tmp for log(abs(S + 1/lambda_last))
#          tmp_sum = logsumexp(c(log(S), -log(abs(lambda_last))), c(1, sign(lambda_last)) )$logR
#          sign_tmp_sum = logsumexp(c(log(S), -log(abs(lambda_last))), c(1, sign(lambda_last)) )$sign
#          
#          logs = c(E[p-1, 1, q]- 2*log(abs(lambda_last)), logTmpExp + tmp_sum + E[p-1, 1, p],
#                   logTmpExp + E[p-1, from_to_values(2, p-1), p] )
#          
#          signs = c(I_pl, S_pl * sign(lambda_last) * sign_tmp_sum, S_pl* signTmpExp * rep(-1, length(from_to_values(2, p-1))) )
#          
#          E[p, k, q] = logsumexp( logs ,signs )$logR 
         
        # (I_pl * E[p-1, 1, q] / (lambda_last^2) ) + 
        # (S_pl  * tmpExp) * ( (S + 1/lambda_last)* E[p-1, 1, p] - tsum)

        logs = c(E[p-1, 1, q], log(abs(lambda_last)) + log(S) -lambda_last*S + E[p-1, 1, p],
                 -lambda_last*S + E[p-1, 1, p], log(abs(lambda_last)) -lambda_last*S  + E[p-1, from_to_values(2, p-1), p] )
        
        signs = c(I_pl, S_pl * sign(lambda_last),  S_pl, S_pl* sign(lambda_last) * rep(-1, length(from_to_values(2, p-1))) )
        
        tmp_log = logsumexp( logs ,signs )
        E[p, k, q] = tmp_log$logR - 2*log(abs(lambda_last))
        
        if(tmp_log$sign == -1)
          error = TRUE
        
      }
    }
  }
  
  # The term q is always one when p > l + 2
  q = 1
  
  
  for( p in  from_to_values(l+2, n) )
  {
    lambda_last = lambdas[p]-lambdas[q]
    for(k in from_to_values(1, p-1)) {
      # no need to abs here
      E[p, k, q] = E[p-1, k, q] - log(lambda_last)
    }
    
    k = p
    E[p, k, q] = E[p-1, 1, q] - 2*log(lambda_last)
  }
  
  tmpE = E[n, 1:n, 1]
  
   P = tmpE[1]
  if(verbose == TRUE){
    print(E)  
  }
  
  
  E_unobserved = tmpE[from_to_values(l+1, n)]
  
  # observed E. reorder them
  E_observed = tmpE[from_to_values(2, l)]
  E_observed = E_observed[original_orders]
  
  if(length(E_observed) > 0 & sum(exp(E_observed-P)) > S) {
    if(verbose == TRUE){
      print("Numerical Error-")
    }
    E_observed = rep(log(S/length(E_observed)) + P, length(E_observed))
  }
  
  
  #return( c(exp(c(P, E_observed, E_unobserved)) ) )
  return( list(finalE=c(P, E_observed, E_unobserved), totalE=E, error = error ) ) 
}
