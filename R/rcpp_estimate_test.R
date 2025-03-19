# Estimate function that calls the underlying Rcpp estimate ekf function
# rcpp_ekf_estimate = function(self, private){
#   
#   message("Estimating using the ekf Rcpp script")
#   
#   # observation/input matrix
#   obsMat = as.matrix(private$data[private$obs.names])
#   inputMat = as.matrix(private$data[private$input.names])
#   
#   # non-na observation matrix
#   numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
#   if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)
#   
#   # number of non-na observations
#   number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)
#   
#   ################################################
#   # Parameters
#   ################################################
#   
#   parameters = sapply(private$parameters, function(x) x[["initial"]]) # Initial parameter values
#   
#   
#   ################################################
#   # Create negative log-likelihood function
#   ################################################
#   nll <- list()
#   nll$par = parameters
#   nll$fn <- function(parVec){
#     ekf_rcpp_likelihood(private$Rcppfunction_f,
#                         private$Rcppfunction_g,
#                         private$Rcppfunction_dfdx,
#                         private$Rcppfunction_h,
#                         private$Rcppfunction_dhdx,
#                         private$Rcppfunction_hvar,
#                         obsMat,
#                         inputMat,
#                         # below is the parameter function argument
#                         parVec,
#                         # above is the parameter function argument
#                         private$initial.state$p0,
#                         private$initial.state$x0,
#                         private$ode.timestep.size,
#                         private$ode.timesteps,
#                         numeric_is_not_na_obsMat,
#                         number_of_available_obs,
#                         private$number.of.states,
#                         private$number.of.observations,
#                         private$ode.solver)$nll
#   }
#   private$nll <- nll
#   
#   ################################################
#   # Optimization
#   ################################################
#   
#   # Parameter Bounds
#   lower.parameter.bound = unlist(lapply(private$free.pars, function(par) par$lower))
#   upper.parameter.bound = unlist(lapply(private$free.pars, function(par) par$upper))
#   if(private$unconstrained.optim){
#     lower.parameter.bound = -Inf
#     upper.parameter.bound = Inf
#   }
#   
#   stats::nlminb(start=nll$par, 
#                 objective=nll$fn, 
#                 lower=lower.parameter.bound,
#                 upper=upper.parameter.bound, 
#                 control=list(trace=1)
#                 )
# 
#   return(invisible(NULL))
# }
