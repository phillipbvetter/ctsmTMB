#######################################################
# MAIN SIMULATE FUNCTION THAT CALL OTHERS
#######################################################
perform_simulation <- function(self, private, use.cpp, n.sims){
  
  
  if(private$method=="laplace"){
    stop("Simulations are not available for the Laplace method")
  }
  
  if(private$method=="laplace.thygesen"){
    stop("Simulations are not available for the Laplace-Thygesen method")
  }
  
  # time prediction
  comptime <- system.time(
    {
      # Predict with Rcpp implementation
      if(use.cpp){
        
        if(!private$silent) message("Simulating with C++...")
        
        lkf_ekf_ukf_simulate_rcpp(private$pars, self, private, n.sims)
        
        # Predict with R implementation
      } else {
        
        if(!private$silent) message("Simulating with R...")
        
        ekf_lkf_ukf_simulate_r(private$pars, self, private, n.sims)
        
      }
      
    }, gcFirst = FALSE)
  
  private$timer_simulation <- comptime
  
  return(invisible(self))
}

#######################################################
# MAIN RETURN PREDICTIONS FUNCTION THAT CALL OTHERS
#######################################################
create_return_simulation <- function(return.k.ahead, n.sims, self, private){
  
  if(!private$silent) message("Returning results...")
  
  list.out = vector("list",length=private$number.of.states)
  names(list.out) = private$state.names
  
  # Compute the prediction times for each horizon
  ran = 0:(private$last.pred.index-1)
  t.j = private$data$t[rep(ran,each=private$n.ahead+1)+1+rep(0:private$n.ahead,private$last.pred.index)]
  t.j.splitlist = split(t.j, ceiling(seq_along(t.j)/(private$n.ahead+1)))
  list.of.time.vectors = lapply(t.j.splitlist, function(x) data.frame(t.j=x))
  
  for(i in seq_along(list.out)){
    list.out[[i]] = stats::setNames(
      lapply(private$simulation, function(ls.outer){
        t(do.call(cbind, lapply(ls.outer, function(ls.inner) ls.inner[,i])))
      }),
      paste0("i",ran)
    )
  }
  
  for(i in seq_along(list.out)){
    for(j in seq_along(list.out[[i]])){
      list.out[[i]][[j]] = data.frame(i = j-1, 
                                      j = (j-1):(j+private$n.ahead-1), 
                                      t.i = rep(private$data$t[i],private$n.ahead+1),
                                      t.j = list.of.time.vectors[[j]][,"t.j"], 
                                      k.ahead = 0:private$n.ahead,
                                      list.out[[i]][[j]]
      )
      nams = paste0(private$state.names,1:n.sims)
      names(list.out[[i]][[j]]) = c("i","j","t.i","t.j","k.ahead",nams)
    }
  }
  
  # Observations
  # eval(parse(text=private$rekf.function.strings$h))
  
  # We should take the state vector, and pass it through the observation vector, right?
  # The inputs should come from the inputMat for the j'th row...
  # There are so many though...
  
  private$simulation = list( states = list.out, observations = list() )
  
  return(invisible(self))
}
