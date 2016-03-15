
#' make_linear_poset
#' @param p the number of events
#' @export
make_linear_poset <- function(p) {
  Poset <- matrix(0, p, p)
  if(p > 1){
    for(i in 1:(p-1))
      Poset[i, i+1] = 1
  }
  Poset
}

#' make_empty_poset
#' @param p the number of events
#' @export
make_empty_poset <- function(p) {
  matrix(0, p, p)
}
#' make_random_poset
#' @param p the number of events
#' @export
make_random_poset <- function(p) {
  poset <- matrix(0, p, p)
  if(p == 4){
    poset[1, 3] = poset[2, 3] = poset[2, 4] = 1
  }
  
  if(p == 6) {
    poset[1, 4] = poset[2, 4] = poset[3, 4] = 1
    poset[3, 5] = 1
    poset[4, 6] = poset[5, 6] = 1
  }
  poset
}


#' sampling.expo
#' @export
sampling.expo <- function(n, mean_s) {
  rexp(n, 1/mean_s)
}

#' sampling.norm
#' @export
sampling.norm <- function(n, mean_s) {
  #rnorm(n, mean_s, 0.1*(abs(mean_s)+1) )
  if(mean_s <= 0) {
    print("Error, mean_s should be positive!")
  }
  rtruncnorm(n, a=0, b=Inf, mean=mean_s, sd=0.1*mean_s )
}

#' sampling.const 
#' @export
sampling.const <- function(n, mean_s) {
  rep(mean_s, n)
}

#' sample_timed_genotypes
#' @export
sample_timed_genotypes <- function(n, poset,  mean_s, lambdas, sampling_fn) {
  T_sampling <- sampling_fn(n, mean_s)
  p = length(lambdas)
  T_events <- matrix(0, n, p)
  
  for(i in 1:p) {
    T_events[, i] <-  rexp(n, lambdas[i])
  }
  
  T_sum_events <- matrix(0, n, p)
  
  #topo_path = topological.sort(graph.adjacency(poset))
  topo_path = my.topological.sort(poset)
   
  for(e in topo_path) {
    parents <- which(poset[, e] == 1)

    if(length(parents) == 0) {
      T_sum_events[, e] = T_events[, e]
    } else if(length(parents) == 1) {
      T_sum_events[, e] = T_events[, e] + T_sum_events[, parents]
    } else {
      T_sum_events[, e] = T_events[, e] + apply(T_sum_events[, parents], 1, max)
    }
    
  }
  
  obs_events <- matrix(0, n, p)
  for(i in 1:p) {
    obs_events[, i] <- as.numeric(T_sum_events[, i] < T_sampling)
  }
  
  genotype_list = apply(obs_events, 1, function(x) { which(x == 1)})
  
  list(n=n, p = p, T_sampling=T_sampling, obs_events=obs_events, genotype_list=genotype_list, mean_s=mean_s, lambdas=lambdas, T_events=T_events, T_sum_events=T_sum_events)  
}
