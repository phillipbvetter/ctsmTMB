#' Stochastic Euler-Maruyama simulation function that calls the underlying Rcpp simulation function
#' 
#' @param self model object
#' @param private model object private fields
#' @param n.sims An integer giving the number of stochastic simulations to be
#' performed
rcpp_simulation = function(self, private, n.sims){
  
  # observation/input matrix
  obsMat = as.matrix(private$data[private$obs.names])
  inputMat = as.matrix(private$data[private$input.names])
  
  # non-na observation matrix
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)
  
  # number of non-na observations
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)
  
  # Call C++ function to perform simulation
  mylist = execute_ekf_simulation(private$Rcppfunction_f,
                                  private$Rcppfunction_g,
                                  private$Rcppfunction_dfdx,
                                  private$Rcppfunction_h,
                                  private$Rcppfunction_dhdx,
                                  private$Rcppfunction_hvar,
                                  obsMat,
                                  inputMat,
                                  private$pars,
                                  private$pred.initial.state$p0,
                                  private$pred.initial.state$x0,
                                  private$ode.timestep.size,
                                  private$ode.timesteps,
                                  private$simulation.timestep.size,
                                  private$simulation.timesteps,
                                  numeric_is_not_na_obsMat,
                                  number_of_available_obs,
                                  private$number.of.states,
                                  private$number.of.observations,
                                  private$number.of.diffusions,
                                  private$last.pred.index,
                                  private$n.ahead,
                                  private$ode.solver,
                                  n.sims)
  
  private$simulation = mylist
  
  return(invisible(NULL))
}

#' Generates a user-friendly data.frame of prediction results from private$prediction
#' 
#' @param return.k.ahead a vector of integers indicating which k.ahead prediction
#' that should be returned out of the 1:k.ahead that were calculated.
#' @param n.sims An integer number indicating the number of stochastic 
#' simulations
#' @param private model object private fields
#' @param self model object
create_return_simulation = function(return.k.ahead, n.sims, self, private){
  
  
  list.out = vector("list",length=private$number.of.states)
  names(list.out) = private$state.names
  
  setRownames = function(obj, nm){rownames(obj) = nm; return(obj)}
  setColnames = function(obj, nm){colnames(obj) = nm; return(obj)}
  
  # Compute the prediction times for each horizon
  ran = 0:(private$last.pred.index-1)
  t.j = private$data$t[rep(ran,each=private$n.ahead+1)+1+rep(0:private$n.ahead,private$last.pred.index)]
  t.j.splitlist = split(t.j, ceiling(seq_along(t.j)/(private$n.ahead+1)))
  list.of.time.vectors = lapply(t.j.splitlist, function(x) data.frame(t.j=x))
  
  for(i in seq_along(list.out)){
    list.out[[i]] = stats::setNames(
      lapply(private$simulation, function(ls.outer){
        # setRownames(
        t(do.call(cbind, lapply(ls.outer, function(ls.inner) ls.inner[,i])))
        # paste0("k.ahead", 0:private$n.ahead)
        # )
      }),
      # paste0("t", head(data$t, private$last.pred.index))
      paste0("i",ran)
    )
  }
  
  for(i in seq_along(list.out)){
    for(j in seq_along(list.out[[i]])){
      list.out[[i]][[j]] = data.frame(i=j-1, 
                                      j=(j-1):(j+private$n.ahead-1), 
                                      t.i=rep(private$data$t[i],private$n.ahead+1),
                                      t.j=list.of.time.vectors[[j]][,"t.j"], 
                                      k.ahead = 0:private$n.ahead,
                                      list.out[[i]][[j]]
      )
      nams = paste0(private$state.names,1:n.sims)
      names(list.out[[i]][[j]]) = c("i","j","t.i","t.j","k.ahead",nams)
      
      # if(!is.null(return.k.ahead)){
        
        # bool = list.out[[i]][[j]][["k.ahead"]] %in% return.k.ahead
        # list.out[[i]][[j]] = list.out[[i]][[j]][bool, ]
      # }
    }
  }
  
  private$simulation = list( states = list.out, observations = list() )
  
  return(invisible(NULL))
}
