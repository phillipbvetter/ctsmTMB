create_rcpp_function_strings = function(self, private){
  
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
  
  private$Rcppfunction_f = code
  
  ##################################################
  # drift jacobian
  ##################################################
  
  # calculate all the terms and substitute variables
  dfdx = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$state.names)){
      # term = hat2pow(Deriv::Deriv(private$diff.terms[[i]]$dt, x=private$state.names[j], cache.exp=F))
      term = hat2pow(ctsmTMB.Deriv(f=private$diff.terms[[i]]$dt, x=private$state.names[j]))
      new.term = do.call(substitute, list(term, subsList))
      dfdx = c(dfdx, sprintf("dfdx(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  
  code = sprintf("SEXP dfdx(Eigen::VectorXd stateVec, Eigen::VectorXd parVec, Eigen::VectorXd inputVec) {
                 Eigen::MatrixXd dfdx(%s,%s);
                 %s
                 return Rcpp::wrap(dfdx);
                 }",private$number.of.states, private$number.of.states, paste(dfdx,collapse=""))
  
  private$Rcppfunction_dfdx = code
  
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
  
  private$Rcppfunction_g <- code
  
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
  
  private$Rcppfunction_h <- code
  
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
  
  private$Rcppfunction_dhdx <- code
  
  
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
  
  private$Rcppfunction_hvar <- code
  
  
}
