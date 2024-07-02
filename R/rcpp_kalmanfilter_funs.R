

create_rcpp_statespace_functions = function(self, private){
  
  # Create substitution translation list
  obsList = lapply(seq_along(private$obs.names), function(id) substitute(obsVec(i),list(i=as.numeric(id-1))))
  parList = lapply(seq_along(private$parameter.names), function(id) substitute(parVec(i),list(i=as.numeric(id-1))))
  stateList = lapply(seq_along(private$state.names), function(id) substitute(stateVec(i),list(i=as.numeric(id-1))))
  inputList = lapply(seq_along(private$input.names), function(id) substitute(inputVec(i),list(i=as.numeric(id-1))))
  names(obsList) = private$obs.names
  names(parList) = private$parameter.names
  names(stateList) = private$state.names
  names(inputList) = private$input.names
  subsList = c(obsList, parList, stateList, inputList)
  
  
  ##################################################
  # drift
  ##################################################
  
  # Perform substitution of parameters, inputs and states
  f = sapply( seq_along(private$state.names),
              function(i){
                drift.term = hat2pow(private$diff.terms[[i]]$dt)
                new.drift.term = do.call(substitute, list(drift.term, subsList))
                sprintf("f(%i) = %s;",i-1,deparse1(new.drift.term))
              })
  
  code = sprintf("SEXP f(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::VectorXd f(%s); 
                 %s
                 return Rcpp::wrap(f);
                 }",private$number.of.states, paste(f,collapse=""))
  private$Rcppfunction_f <- RcppXPtrUtils::cppXPtr(code,
                                                   depends=c("RcppEigen","sdeTMB"))
  
  ##################################################
  # drift jacobian
  ##################################################
  
  # calculate all the terms and substitute variables
  dfdx = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$state.names)){
      term = hat2pow(Deriv::Deriv(private$diff.terms[[i]]$dt,x=private$state.names[j], cache.exp=F))
      new.term = do.call(substitute, list(term, subsList))
      dfdx = c(dfdx, sprintf("dfdx(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  
  code = sprintf("SEXP dfdx(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd dfdx(%s,%s);
                 %s
                 return Rcpp::wrap(dfdx);
                 }",private$number.of.states, private$number.of.states, paste(dfdx,collapse=""))
  
  private$Rcppfunction_dfdx <- RcppXPtrUtils::cppXPtr(code, depends="RcppEigen")
  
  ##################################################
  # diffusion
  ##################################################
  
  # calculate all the terms and substitute variables
  g = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$diff.processes[-1])){
      term = hat2pow(private$diff.terms[[i]][[j+1]])
      new.term = do.call(substitute, list(term, subsList))
      g = c(g, sprintf("g(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  
  code = sprintf("SEXP g(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd g(%s,%s);
                 %s
                 return Rcpp::wrap(g);
                 }",private$number.of.states, private$number.of.diffusions, paste(g,collapse=""))
  
  private$Rcppfunction_g <- RcppXPtrUtils::cppXPtr(code, depends=c("RcppEigen"))
  
  ##################################################
  # observation
  ##################################################
  
  # calculate all the terms and substitute variables
  h = sapply(seq_along(private$obs.names), 
             function(i){
               term = hat2pow(private$obs.eqs.trans[[i]]$rhs)
               new.term = do.call(substitute, list(term, subsList))
               sprintf("h(%s) = %s;",i-1, deparse1(new.term))
             }) 
  
  code = sprintf("SEXP h(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::VectorXd h(%s);
                 %s
                 return Rcpp::wrap(h);
                 }",private$number.of.observations, paste(h,collapse=""))
  
  private$Rcppfunction_h <- RcppXPtrUtils::cppXPtr(code, depends="RcppEigen")
  
  ##################################################
  # observation jacobian
  ##################################################
  
  # calculate all the terms and substitute variables
  dhdx = c()
  for(i in seq_along(private$obs.names)){
    for(j in seq_along(private$state.names)){
      term = hat2pow(private$diff.terms.obs[[i]][[j]])
      new.term = do.call(substitute, list(term, subsList))
      dhdx = c(dhdx, sprintf("dhdx(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  
  code = sprintf("SEXP dhdx(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd dhdx(%s,%s);
                 %s
                 return Rcpp::wrap(dhdx);
                 }", private$number.of.observations, private$number.of.states, paste(dhdx,collapse=""))
  private$Rcppfunction_dhdx <- RcppXPtrUtils::cppXPtr(code, 
                                                      depends="RcppEigen")
  
  
  ##################################################
  # observation variance
  ##################################################
  
  hvar = lapply(seq_along(private$obs.var.trans), 
                function(i) {
                  term = hat2pow(private$obs.var.trans[[i]]$rhs)
                  new.term = do.call(substitute, list(term, subsList))
                  sprintf("hvar(%s,%s) = %s;", i-1, i-1, deparse1(new.term))
                })
  
  code = sprintf("SEXP hvar(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd hvar(%s,%s);
                 hvar.setZero();
                 %s
                 return Rcpp::wrap(hvar);
                 }", private$number.of.observations, private$number.of.observations, paste(hvar,collapse=""))
  
  private$Rcppfunction_hvar <- RcppXPtrUtils::cppXPtr(code, depends="RcppEigen")
  
}



perform_rcpp_ekf_prediction = function(self, private, pars){
  
  # observation matrix
  obsMat = as.matrix(private$data[private$obs.names])
  
  # non-na observation matrix
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)
  
  # number of non-na observations
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)
  
  # predict using c++ function
  mylist = execute_ekf_prediction(private$Rcppfunction_f,
                                  private$Rcppfunction_g,
                                  private$Rcppfunction_dfdx,
                                  private$Rcppfunction_h,
                                  private$Rcppfunction_dhdx,
                                  private$Rcppfunction_hvar,
                                  as.matrix(private$data[private$obs.names]),
                                  as.matrix(private$data[private$input.names]),
                                  pars,
                                  private$pred.initial.state$p0,
                                  private$pred.initial.state$x0,
                                  private$ode.timestep.size,
                                  private$ode.timesteps,
                                  numeric_is_not_na_obsMat,
                                  number_of_available_obs,
                                  private$number.of.states,
                                  private$number.of.observations,
                                  private$last.pred.index,
                                  private$n.ahead,
                                  private$ode.solver)
  
  return(mylist)
}

perform_rcpp_ekf_simulation = function(self, private, pars, initial.state, n.sims){
  
  # Calculate NA-vectors needed for update step in Kalman filter
  obsMat = as.matrix(private$data[private$obs.names])
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)
  
  # Call C++ function to perform simulation
  mylist = execute_ekf_simulation(private$Rcppfunction_f,
                                  private$Rcppfunction_g,
                                  private$Rcppfunction_dfdx,
                                  private$Rcppfunction_h,
                                  private$Rcppfunction_dhdx,
                                  private$Rcppfunction_hvar,
                                  as.matrix(private$data[private$obs.names]),
                                  as.matrix(private$data[private$input.names]),
                                  pars,
                                  initial.state$p0,
                                  initial.state$x0,
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
  
  return(mylist)
}

perform_prediction = function(parVec, self, private){
  
  # observation matrix
  obsMat = as.matrix(private$data[private$obs.names])
  
  # non-na observation matrix
  numeric_is_not_na_obsMat = t(apply(obsMat, 1, FUN=function(x) as.numeric(!is.na(x))))
  if(nrow(numeric_is_not_na_obsMat)==1) numeric_is_not_na_obsMat = t(numeric_is_not_na_obsMat)
  
  # number of non-na observations
  number_of_available_obs = apply(numeric_is_not_na_obsMat, 1, sum)
  
  # predict using c++ function
  mylist <- execute_ekf_prediction(private$Rcppfunction_f,
                                   private$Rcppfunction_g,
                                   private$Rcppfunction_dfdx,
                                   private$Rcppfunction_h,
                                   private$Rcppfunction_dhdx,
                                   private$Rcppfunction_hvar,
                                   as.matrix(private$data[private$obs.names]),
                                   as.matrix(private$data[private$input.names]),
                                   parVec,
                                   private$pred.initial.state$p0,
                                   private$pred.initial.state$x0,
                                   private$ode.timestep.size,
                                   private$ode.timesteps,
                                   numeric_is_not_na_obsMat,
                                   number_of_available_obs,
                                   private$number.of.states,
                                   private$number.of.observations,
                                   private$last.pred.index,
                                   private$n.ahead,
                                   private$ode.solver)
  
  
  return(mylist)
  
  
}
