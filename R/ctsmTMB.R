#' @title Methods for the 'ctsmTMB' R6 class
#' 
#' @description The following public methods are used to construct a stochastic state space model 
#' system, consisting of a set of stochastic differential equations (SDEs), and one or more algebraic observation 
#' equations (AOEs). The AOEs are used to infer information about the value of the (latent) states governed by the SDEs, and
#' thus must be functions of at least one state.
#' 
#' @returns The function returns an object of class \code{R6} and \code{ctsmTMB}, 
#' which can be used to define a stochastic state space system.
#' 
#' @docType package
#' @name ctsmTMB
#' @export
ctsmTMB = R6::R6Class(
  
  # Class name
  classname = "ctsmTMB",
  
  
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  # Public Methods
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  
  public = list(
    
    ########################################################################
    # INITIALIZE FIELDS
    ########################################################################
    #' @description 
    #' Initialize private fields
    initialize = function() {
      # modelname and path
      private$modelname = "ctsmTMB_model"
      private$cppfile.directory = paste(getwd(),"/","ctsmTMB_cppfiles", sep="")
      private$cppfile.path = paste(private$cppfile.directory,"/",private$modelname,sep="")
      private$cppfile.path.with.method = NULL
      private$modelname.with.method = NULL
      
      # model equations
      private$sys.eqs = NULL
      private$obs.eqs = NULL
      private$obs.var = NULL
      private$alg.eqs = NULL
      private$inputs = list(t=list(name="t",input=quote(t)))
      private$parameters = NULL
      private$initial.state = NULL
      private$pred.initial.state = list(mean=NULL,cov=NULL)
      private$tmb.initial.state.for.parameters = NULL
      private$iobs = NULL
      
      # after algebraics
      private$sys.eqs.trans = NULL
      private$obs.eqs.trans = NULL
      private$obs.var.trans = NULL
      
      # options
      private$method = "ekf"
      private$use.hessian = FALSE
      private$state.dep.diff = FALSE
      private$lamperti = list(transform="identity",states=NULL)
      private$compile = FALSE
      private$loss = list(loss=0L,c=3)
      private$tukey.pars = rep(0,4)
      private$silent = FALSE
      private$map = NULL
      private$control.nlminb = list()
      private$ode.solver = NULL
      private$unconstrained.optim = NULL
      
      # hidden
      private$lock.model = FALSE
      private$fixed.pars = NULL
      private$free.pars = NULL
      
      # names
      private$state.names = NULL
      private$obs.names = NULL
      private$obsvar.names = NULL
      private$input.names = "t"
      private$parameter.names = NULL
      
      # lengths
      private$number.of.states = 0
      private$number.of.observations = 0
      private$number.of.diffusions = 0
      private$number.of.pars = 0
      private$number.of.inputs = 0
      
      # differentials
      private$diff.processes = NULL
      private$diff.terms = NULL
      private$diff.terms.obs = NULL
      
      # data, nll, opt
      private$data = NULL
      private$nll = NULL
      private$opt = NULL
      private$fit = NULL
      
      # predict
      private$n.ahead = 0
      private$last.pred.index = 0
      
      # Rcpp Prediction
      private$Rcppfunction_f = NULL
      private$Rcppfunction_g = NULL
      private$Rcppfunction_dfdx = NULL
      private$Rcppfunction_h = NULL
      private$Rcppfunction_dhdx = NULL
      private$Rcppfunction_hvar = NULL
      
      # unscented transform
      private$ukf_alpha = NULL
      private$ukf_beta = NULL
      private$ukf_kappa = NULL
    },
    
    ########################################################################
    # ADD SYSTEMS
    ########################################################################
    #' @description 
    #' Define and add multiple stochastic differential equation governing the 
    #' process of individual state variables on the form 
    #' 
    #' \code{d<state> ~ f(t,<states>,<inputs>) * dt + g1(t,<states>,<inputs>) * dw1 
    #' + g2(t,<states>,<inputs>) * dw2 + ... + gN(t,<states>,<inputs>) * dwN}
    #'                                            
    #' where \code{f} is the drift, and \code{g1, g2, ..., gN} are diffusions, with 
    #' differential brownian motions dw1, dw2, ..., dwN.
    #' 
    #' @examples 
    #' # Specify Ornstein-Uhlenbeck Process
    #' add_systems(dx ~ theta * (mu - x + u) * dt + sigma * dw)
    #'              
    #' @param form formula specifying the stochastic differential equation to be 
    #' added to the system.
    #' @param ... additional formulas similar to \code{form} for specifying 
    #' multiple equations at once.
    #' 
    add_systems = function(form,...) {
      
      if(private$lock.model){
        stop("The model is locked after applying algebraics")
      }
      
      # store each provided formula
      lapply(c(form,...), function(form) {
        
        # Check if the system equation is valid
        result = check_system_eqs(form, self, private)
        
        # Check if name is not used for something else
        check_if_name_is_overwritable(result$name, "state", self, private)
        
        # Update equations and names
        private$sys.eqs[[result$name]] = result
        private$state.names = names(private$sys.eqs)
        
      })
      return(invisible(self))
    },
    
    ########################################################################
    # ADD OBSERVATIONS
    ########################################################################
    #' @description
    #' Define and add a relationship between an observed variable and system states. 
    #' The observation equation takes the form
    #' 
    #' \code{<observation> ~ h(t,<states>,<inputs>) + e)}
    #' 
    #' where \code{h} is the observation function, and \code{e} is normally 
    #' distributed noise with zero mean and variance to be specified. The 
    #' observation variable should be present in the data provided when calling
    #' \code{estimate(.data)} for parameter estimation.
    #'  
    #' @examples
    #' #Specify observation directly as a latent state
    #' add_observations(y ~ x)
    #' 
    #' Specify observation as the sum of exponentials of two latent states
    #' add_observations(y ~ exp(x1) + exp(x2))
    #' @param form formula class specifying the obsevation equation to be added to the system.
    #' @param ... additional formulas identical to \code{form} to specify multiple observation equations at a time.
    #' @param obsnames character vector specifying the name of the observation. When the observation left-hand side
    #' consists of more than just a single variable name (when its class is 'call' instead of 'name') it will be 
    #' given a name on the form obs__# where # is a number, unless obsnames is provided.
    add_observations = function(form,...,obsnames=NULL) {
      
      if(private$lock.model){
        stop("The model is locked after applying algebraics")
      }
      
      # Check obsnames
      if(!is.null(obsnames)){
        if(length(c(form,...))!=length(obsnames)){
          stop("You must supply as many observation names as there are observation equations")
        }
        if(!is.character(obsnames)){
          stop("The observation names in obsnames must be characters")
        }
      }
      
      # attach observation names to the each formula in list
      formlist = lapply(c(form,...), function(form) list(form=form))
      for(i in seq_along(formlist)){formlist[[i]]$name = obsnames[i]}
      
      # store each provided formula
      lapply(formlist, function(forms) {
        # Check if the equation is valid
        result = check_observation_eqs(forms, self, private)
        
        # Check if name is not used for something else
        check_if_name_is_overwritable(result$name, "obs", self, private)
        
        # Update equations and names
        private$obs.eqs[[result$name]] = result
        private$obs.names = names(private$obs.eqs)
        
        # Create space in the observation variances for the observation name
        private$obs.var[[result$name]] = list()
      })
      
      return(invisible(self))
    },
    
    ########################################################################
    # ADD OBSERVATION VARIANCES
    ########################################################################
    #' @description Specify the variance of an observation equation.
    #' 
    #' A defined observation variable \code{y} in e.g. \code{add_observations(y ~ 
    #' h(t,<states>,<inputs>)} is pertubed by Gaussian noise with zero mean and 
    #' variance 
    #' to-be specified using \code{add_observation_variances(y ~ p(t,<states>,<inputs>)}. 
    #' We can for instance declare \code{add_observation_variances(y ~ sigma_x^2} 
    #' where \code{sigma_x} is a fixed effect parameter to be declared through 
    #' \code{add_parameters}.
    #' 
    #' @param form formula class specifying the obsevation equation to be added 
    #' to the system.
    #' @param ... additional formulas identical to \code{form} to specify multiple 
    #' observation equations at a time.
    #' 
    add_observation_variances = function(form,...) {
      
      if(private$lock.model){
        stop("The model is locked after applying algebraics")
      }
      
      # store each provided formula
      lapply(c(form,...), function(form) {
        # Check if the equation is valid
        result = check_observation_variance_eqs(form, self, private)
        
        # Check if name is not used for something else
        check_if_name_is_overwritable(result$name, "obsvar", self, private)
        
        # Update equations and names
        private$obs.var[[result$name]] = result
        private$obsvar.names = names(private$obs.var)
        
      })
      
      return(invisible(self))
    },
    
    ########################################################################
    # ADD INPUTS
    ########################################################################
    #' @description Declare variables as data inputs
    #' 
    #' Declare whether a variable contained in system, observation or observation 
    #' variance equations is an input variable. If e.g. the system equation contains 
    #' an input variable \code{u} then it is declared using \code{add_inputs(u)}. 
    #' The input \code{u} must be contained in the data.frame \code{.data} provided 
    #' when calling the \code{estimate} or \code{predict} methods.
    #' 
    #' @param ... variable names that specifies the name of input variables in the defined system.
    #' 
    add_inputs =  function(...) {
      
      # if(private$lock.model){
      #   stop("The model is locked after applying algebraics")
      # }
      
      args = as.list(match.call()[-1])
      
      # lapply over all parsed inputs
      lapply(args, function(args) {
        
        # Check if the equation is valid
        result = check_inputs(args, self, private)
        
        # Check if name is not used for something else
        check_if_name_is_overwritable(result$name, "input", self, private)
        
        # Update equations and names
        private$inputs[[result$name]] = result
        private$input.names = names(private$inputs)
      })
      #
      return(invisible(self))
    },
    
    ########################################################################
    # ADD PARAMETERS
    ########################################################################
    #' @description Declare which variables that are (fixed effects) parameters in
    #' the specified model, and specify the initial optimizer value, as well as
    #' lower / upper bounds during optimization. There are two ways to declare parameters:
    #' 
    #' 1. You can declare parameters using formulas i.e. \code{add_parameters( 
    #' theta = c(1,0,10), mu = c(0,-10,10) )}. The first value is the initial 
    #' value for the optimizer, the second value is the lower optimization 
    #' bound and the third value is the upper optimization bound. 
    #' 
    #' 2. You can provide a 3-column matrix where rows corresponds to different 
    #' parameters, and the parameter names are provided as rownames of the matrix. 
    #' The columns values corresponds to the description in the vector format above.
    #'
    #' @param ... a named vector or matrix as described above.
    add_parameters = function(...) {
      
      if(nargs()==0L){stop("No arguments received")}
      
      arglist = list(...)
      argnames = names(arglist)
      
      
      # run over each parameter argument either a vector or a matrix
      lapply(seq_along(arglist), function(i) {
        
        # grab vector or matrix from list of parsed arguments
        par.entry = arglist[[i]]
        par.name = argnames[i]
        
        ##### VECTOR INPUTS #####
        if(is.vector(par.entry)){
          
          # check basics
          par.entry = check_parameter_vector(par.entry, par.name, self, private)
          
          # check name
          check_if_name_is_overwritable(par.name, "pars", self, private)
          
          # set expected names names
          expected.names = c("initial","lower","upper")
          
          # store in parameter list ordered
          private$parameters[[par.name]] = as.list(par.entry[expected.names])
          
          # set or remove a fixed parameter (NA-bounds)
          private$fixed.pars[[par.name]] = NULL
          if (all(is.na(par.entry[c("lower","upper")]))){
            private$fixed.pars[[par.name]] = private$parameters[[par.name]]
            private$fixed.pars[[par.name]][["factor"]] = factor(NA)
          }
          
          # update parameter names
          private$parameter.names = names(private$parameters)
          
          ##### MATRIX/DATA.FRAME INPUTS #####  
        } else if (is.matrix(par.entry) | is.data.frame(par.entry)){
          
          
          par.entry = check_parameter_matrix(par.entry, self, private)
          
          parnames = rownames(par.entry)
          
          # lapply over all matrix rows
          lapply( 1:nrow(par.entry), function(i) {
            parname = parnames[i]
            
            # basic validity checks
            check_parameter_vector(par.entry[i,], parname, self, private)
            
            # check name
            check_if_name_is_overwritable(parname, "pars", self, private)
            
            # store in parameter list
            private$parameters[[parname]] = list(initial = par.entry[i,"initial"], 
                                                 lower = par.entry[i,"lower"], 
                                                 upper = par.entry[i,"upper"])
            
            # set or remove a fixed parameter (NA-bounds)
            private$fixed.pars[[parname]] = NULL
            if(all(is.na(par.entry[i,c("lower","upper")]))){
              private$fixed.pars[[parname]] = private$parameters[[parname]]
              private$fixed.pars[[parname]][["factor"]] = factor(NA)
            }
            
            # update parameter names
            private$parameter.names = names(private$parameters)
            return(invisible(self))
          })
          
          ##### ELSE STOP #####  
        } else {
          
          stop("You can only supply parameter vectors or matrices/data.frames")
          
        }
        
      })
      
      # return
      return(invisible(self))
    },
    
    ########################################################################
    # ADD ALGEBRAICS
    ########################################################################
    #' @description Add algebraic relations.
    #' 
    #' Algebraic relations is a convenient way to transform parameters in your equations.
    #' In the Ornstein-Uhlenbeck process the rate parameter \code{theta} is always positive, so
    #' estimation in the log-domain is a good idea. Instead of writing \code{exp(theta)} directly
    #' in the system equation one can transform into the log domain using the algebraic relation
    #' \code{add_algebraics(theta ~ exp(logtheta))}. All instances of \code{theta} is replaced
    #' by \code{exp(logtheta)} when compiling the C++ function. Note that you must provide values
    #' for \code{logtheta} now instead of \code{theta} when declaring parameters through 
    #' \code{add_parameters}
    #' 
    #' @param form formula specifying the stochastic differential equation(s) to be added to the system.
    #' @param ... additional formulas similar to \code{form} for specifying multiple equations at once.
    add_algebraics = function(form,...) {
      
      # You are not allowed to change the model after you add algebraics
      private$lock.model = TRUE
      
      lapply(c(form,...), function(form) {
        
        # Check if the equation is valid
        result = check_algebraics(form, self, private)
        
        # Update equations
        private$alg.eqs[[result$name]] = result
        
        # Remove potentially redefined parameters
        remove_parameter(result$name, self, private)
      })
      
      return(invisible(self))
    },
    
    ########################################################################
    # SET INITIAL STATE
    ########################################################################
    #' @description Declare the initial state values i.e. mean and covariance for the system states.
    #' 
    #' @param initial.state a named list of two entries 'x0' and 'p0' containing the initial state and covariance of the state
    #' @param estimate boolean value which indicates whether or not the initial conditions
    #' shall be estimated as fixed effects parameters. The provided mean and covariance are then
    #' used as initial guesses
    #' 
    set_initial_state = function(initial.state, estimate=FALSE) {
      if (is.null(private$sys.eqs)) {
        stop("Please specify system equations first")
      }
      
      if (!is.list(initial.state)) {
        stop("Please provide a list!")
      }
      
      if (length(initial.state) != 2) {
        stop("Please provide a list of length two")
      }
      
      # unpack list items
      x0 = initial.state[[1]]
      p0 = initial.state[[2]]
      
      if (!is.numeric(x0)) {
        stop("The mean vector is not a numeric")
      }
      
      if (any(is.na(x0))) {
        stop("The mean vector contains NAs.")
      }
      
      if (length(x0)!=length(private$sys.eqs)) {
        stop("The initial state vector should have length ",length(private$sys.eqs))
      }
      
      if (!all(dim(p0)==c(length(private$sys.eqs),length(private$sys.eqs)))) {
        stop("The covariance matrix should be square with dimension ", length(private$sys.eqs))
      }
      
      # convert scalar to matrix
      if(!is.matrix(p0) & is.numeric(p0) & length(p0)==1){
        p0 = p0 * diag(1)
      }
      
      if (!is.numeric(p0)) {
        stop("The covariance matrix is not a numeric")
      }
      
      if (any(is.na(p0))) {
        stop("The covariance matrix contains NAs")
      }
      
      if (any(eigen(p0)$values < 0)){
        stop("The covariance matrix is not positive semi-definite")
      }
      
      if (!isSymmetric.matrix(p0)){
        stop("The covariance matrix is symmetric")
      }
      
      private$initial.state = list(x0=x0, p0=as.matrix(p0))
      
      # if(estimate){
      # then we should add the states as fixed effects parameters?
      # }
      
      return(invisible(self))
    },
    
    ########################################################################
    # SET LAMPERTI TRANSFORMATION
    ########################################################################
    #' @description Set a Lamperti Transformation
    #'
    #' If the provided system equations have state dependent diffusion in a few available ways
    #' then it is advantageous to perform a transformation to remove the state dependence. This 
    #' comes at the cost of a more complicated drift function. The following types of state-dependence
    #' is currently supported
    #'
    #' 1. 'identity' - when the diffusion is state-independent (default)
    #' 2. 'log' - when the diffusion is proportional to to x * dw
    #' 3. 'logit' - when the diffusion is proportional to x * (1-x) * dw
    #' 4. 'sqrt-logit' - when the diffusion is proportional to sqrt(x * (1-x)) * dw
    #' 
    #' @param transforms character vector - one of either "identity, "log", "logit", "sqrt-logit"
    #' @param states a vector of the state names for which the specified transformations should be applied to. 
    set_lamperti = function(transforms, states=NULL) {
      
      # remove repeated entries
      states = unique(states)
      
      # Check if transformation is a string
      if (!(is.character(transforms))) {
        stop("Error: You must pass a (vector of) string(s)")
      }
      
      # select all states if states=NULL
      if(is.null(states)){
        states = private$state.names
      }
      
      # check if parsed state names are strings
      if (!(is.character(states))) {
        stop("Error: You must pass a vector of state names")
      }
      
      #  do state names exist?
      bool = !(states %in% names(private$sys.eqs))
      if (any(bool)) {
        stop("The following state names don't exist: \n\t ",paste(states[bool],collapse=", "))
      }
      
      # The length of states and transformations must be equal
      if( length(transforms) != length(states)){
        
        # Recycle if 1 transformation is supplied
        if(length(transforms)==1){
          
          warning("You provided fewer transforms than states - recycling transformation")
          transforms = rep(transforms,length(states))
          
        } else {
          
          # else throw an error
          stop("Error: Mismatching number of transformations and states. You must pass either 
               1) one transformation for all states, 
               2) as many transformations as states - one for each of them.")
        }
        
      }
      
      if(length(transforms) < length(states)){
        
        warning("You provided fewer transforms than states - recycling transformation")
        transforms = rep(transforms,length(states))
        
      }
      
      # check if requested transform is valid
      available_transforms = c("identity","log","logit","sqrt-logit")
      if (!all(transforms %in% available_transforms)) {
        
        stop("You requested a transform that is not available. Choose among these:
             1. 'identity'
             2. 'log'
             3. 'logit'
             4. 'sqrt-logit'")
        
      }
      
      # Store the transformation
      private$lamperti = list(transforms=transforms, states=states)
      
      # return
      return(invisible(self))
    },
    
    ########################################################################
    # SET MODEL NAME
    ########################################################################
    #' @description Set modelname used to create the C++ file for TMB
    #'
    #' When calling \code{TMB::MakeADFun} the (negative log) likelihood function 
    #' is created in the directory specified by the \code{set_cppfile_directory} 
    #' method with name \code{<modelname>.cpp}
    #' 
    #' @param name string defining the model name.
    set_modelname = function(name) {
      
      # was a string passed?
      if (!is.character(name)) {
        stop("The modelname must be a string")
      }
      
      # set options
      private$modelname = name
      private$cppfile.path = paste(private$cppfile.directory,"/",name,sep="")
      
      # return
      return(invisible(self))
    },
    
    ########################################################################
    # SET CPP FILE DIRECTORY
    ########################################################################
    #' @description Set the path directory where the constructed C++ file is created.
    #' You should specify the entire path, unless you want to construct a subfolder
    #' in the current working directory - then you can call e.g. 
    #' \code{set_cppfile_directory("folder_in_current_wd")}.
    #'
    #' @param directory string specifying a local directory
    set_cppfile_directory = function(directory) {
      
      # Check if string
      if (!is.character(directory)) {
        stop("You must pass a string")
      }
      
      # create directory if it does not exist
      if (!dir.exists(directory)) {
        message("The specified directory does not exist - creating it")
        dir.create(directory)
      }
      private$cppfile.directory = directory
      
      # update private$cppfile.path by calling set_modelname
      self$set_modelname(private$modelname)
      
      # return
      return(invisible(self))
    },
    
    ########################################################################
    # SET MAXIMUM A POSTERIORI 
    ########################################################################
    #' @description Enable maximum a posterior (MAP) estimation.
    #'
    #' Adds a maximum a posterior contribution to the (negative log) likelihood 
    #' function by evaluating the fixed effects parameters in a multivariate Gaussian 
    #' with \code{mean} and \code{covariance} as provided.
    #' 
    #' @param mean mean vector of the Gaussian prior parameter distribution
    #' @param cov covariance matrix of the Gaussian prior parameter distribution
    set_map = function(mean,cov) {
      
      # Test the inputs
      if (!is.numeric(mean)) {
        stop("The MAP mean vector is not numeric")
      }
      
      if (length(mean)!=length(private$parameters)) {
        stop("The MAP parameter vector should have length ",length(private$parameters))
      }
      
      if (!is.matrix(cov)) {
        stop("The MAP covariance matrix is not of class matrix")
      }
      
      if (!all(dim(cov)==rep(length(private$parameters),2))) {
        stop("The MAP covariance matrix should be square with dimension ", length(private$parameters))
      }
      
      # Store the mean and covariance
      private$map = list(mean=mean,cov=cov)
      
      # Return
      return(invisible(self))
    },
    
    ########################################################################
    # GET SYSTEMS
    ########################################################################
    #' @description Retrieve system equations.
    get_systems = function() {
      
      # extract system formulas
      syseqs = lapply(private$sys.eqs,function(x) x$form)
      
      # return
      return(syseqs)
    },
    
    ########################################################################
    # GET OBSERVATIONS
    ########################################################################
    #' @description Retrieve observation equations.
    get_observations = function() {
      
      # extract observation formulas
      obseqs = lapply(private$obs.eqs,function(x) x$form)
      
      # return
      return(obseqs)
    },
    
    ########################################################################
    # GET OBSERVATION VARIANCES
    ########################################################################
    #' @description Retrieve observation variances
    get_observation_variances = function() {
      
      # extract observation variance formulas
      obsvar = lapply(private$obs.var,function(x) x$form)
      
      # return
      return(obsvar)
    },
    
    ########################################################################
    # GET ALGEBRAICS
    ########################################################################
    #' @description Retrieve algebraic relations
    get_algebraics = function() {
      
      # extract algebraic relation formulas
      algs = lapply(private$alg.eqs,function(x) x$form)
      
      # return
      return(algs)
    },
    
    ########################################################################
    # GET ALGEBRAICS
    ########################################################################
    #' @description Retrieve initially set state and covariance
    get_initial_state = function() {
      
      
      # extract algebraic relation formulas
      initial.state = private$initial.state
      
      # return
      return(initial.state)
    },
    
    ########################################################################
    # GET PARAMETER MATRIX
    ########################################################################
    #' @description Get initial (and estimated) parameters.
    #' @param type one of "all", free" or "fixed" parameters.
    #' @param value one of "all", initial", "estimate", "lower" or "upper"
    get_parameters = function(type="all", value="all") {
      
      if(is.null(private$parameters)){
        return(NULL)
      }
      
      # create return matrix
      .df = data.frame(matrix(NA,nrow=length(private$parameters),ncol=5))
      names(.df) = c("type","estimate", "initial","lower","upper")
      rownames(.df) = private$parameter.names
      
      # put parameters into it
      .df[,c("initial","lower","upper")] = t(sapply(private$parameters,unlist))
      .df[["type"]] = "free"
      .df[["estimate"]] = NA
      .df[names(private$fixed.pars),"type"] = "fixed"
      .df[names(private$fixed.pars),"estimate"] = sapply(private$fixed.pars,function(x) x$initial)
      
      # if the fit exists then assign the free (estimate) parameters the estimated values from fit$par.fixed
      if(!is.null(private$fit$par.fixed)){
        .df[names(private$free.pars),"estimate"] = private$fit$par.fixed[names(private$free.pars)]
      }
      # if(is.null(private$fit)){
      #   # remove estimate if not fit has been generated
      #   # .df = .df[,-2]
      # } else {
      #   # if the fit exists then assign the free (estimate) parameters the estimated values
      #   # from fit$par.fixed
      #   .df[names(private$free.pars),"estimate"] = private$fit$par.fixed[names(private$free.pars)]
      # }
      
      # Filter rows by free or fixed parameter types
      .df = switch(type,
                   free = {
                     .df[.df[["type"]] == "free",]
                   },
                   fixed = {
                     .df[.df[["type"]] == "fixed",]
                   },
                   all = {
                     .df
                   })
      
      
      # Filter columns by value 
      .df = switch(value,
                   initial = {
                     .df[,"initial",drop=T]
                   },
                   lower = {
                     .df[,"lower",drop=T]
                   },
                   upper = {
                     .df[,"upper",drop=T]
                   },
                   estimate = {
                     .df[,"estimate",drop=T]
                   },
                   all = {
                     .df
                   })
      
      # return
      return(.df)
    },
    
    ########################################################################
    # CONSTRUCT NEG. LIKELIHOOD FUNCTION HANDLERS FROM TMB
    ########################################################################
    #' @description Construct and extract function handlers for the negative
    #' log likelihood function.
    #'
    #' The handlers from \code{TMB}'s \code{MakeADFun} are constructed and returned. 
    #' This enables the user to e.g. choose their own optimization algorithm, or just
    #' have more control of the optimization workflow.
    #' 
    #' @param data a data.frame containing time-vector 't', observations and inputs. 
    #' The observations can take \code{NA}-values.  
    #' @param ode.timestep the time-step used in the filtering schemes. The
    #' time-step has two different uses depending on the chosen method.
    #' 
    #' 1. Kalman Filters: The time-step is used when numerically solving the 
    #' moment differential equations.
    #' 2. Laplace Approximation: The time-step is used in the Euler-Maruyama
    #' simulation scheme for simulating a sample path of the stochastic differential 
    #' equation, which serves to link together the latent (random effects) states.
    #' 
    #' The defined step-size is used to calculate the number of steps between 
    #' observation time-points as 
    #' defined by the provided \code{data}. If the calculated number of steps is larger than N.01 where N 
    #' is an integer, then the time-step is reduced such that exactly N+1 steps is taken between observations  
    #' The step-size is used in the two following ways depending on the
    #' chosen method:
    #' 1. Kalman filters: The time-step is used as the step-size in the
    #' numerical Forward-Euler scheme to compute the prior state mean and
    #' covariance estimate as the final time solution to the first and second
    #' order moment differential equations.
    #' 2. TMB method: The time-step is used as the step-size in the Euler-Maruyama
    #' scheme for simulating a sample path of the stochastic differential equation,
    #' which serves to link together the latent (random effects) states.
    #' @param ode.solver Sets the ODE solver used in the Kalman Filter methods for solving the moment 
    #' differential equations. The default "euler" is the Forward Euler method, alternatively the classical
    #' 4th order Runge Kutta method is available via "rk4".
    #' @param compile boolean value. The default (\code{FALSE}) is to not compile the C++ objective
    #' function but assume it is already compiled and corresponds to the specified model object. It is
    #' the user's responsibility to ensure correspondence between the specified model and the precompiled
    #' C++ object. If a precompiled C++ object is not found in the specified directory i.e. 
    #' in \code{<cppfile_directory>/<modelname>/(dll/so)} then the compile flag is set to \code{TRUE}.
    #' If the user makes changes to system equations, observation equations, observation variances, 
    #' algebraic relations or lamperi transformations then the C++ object should be recompiled.
    #' @param method character vector - one of either "ekf", "ukf" or "tmb". Sets the estimation 
    #' method. The package has three available methods implemented:
    #' 1. The natural TMB-style formulation where latent states are considered random effects
    #' and are integrated out using the Laplace approximation. This method only yields the gradient
    #' of the (negative log) likelihood function with respect to the fixed effects for optimization.
    #' The method is slower although probably has some precision advantages, and allows for non-Gaussian
    #' observation noise (not yet implemented). One-step / K-step residuals are not yet available in
    #' the package.
    #' 2. (Continous-Discrete) Extended Kalman Filter where the system dynamics are linearized
    #' to handle potential non-linearities. This is computationally the fastest method.
    #' 3. (Continous-Discrete) Unscented Kalman Filter. This is a higher order non-linear Kalman Filter
    #' which improves the mean and covariance estimates when the system display high nonlinearity, and
    #' circumvents the necessity to compute the jacobian of the drift and observation functions.
    #' 
    #' All package features are currently available for the kalman filters, while TMB is limited to
    #' parameter estimation. In particular, it is straight-forward to obtain k-step-ahead predictions
    #' with these methods (use the \code{predict} S3 method), and stochastic simulation is also available 
    #' in the cases where long prediction horizons are sought, where the normality assumption will be 
    #' inaccurate.
    #' @param unscented_hyperpars the three hyper-parameters \code{alpha}, \code{beta} and \code{kappa} defining
    #' the unscented transformation.
    #' @param loss character vector. Sets the loss function type (only implemented for the kalman filter
    #' methods). The loss function is per default quadratic in the one-step residauls as is natural 
    #' when the Gaussian (negative log) likelihood is evaluated, but if the tails of the 
    #' distribution is considered too small i.e. outliers are weighted too much, then one 
    #' can choose loss functions that accounts for this. The three available types available:
    #' 
    #' 1. Quadratic loss (\code{quadratic}).
    #' 2. Quadratic-Linear (\code{huber})
    #' 3. Quadratic-Constant (\code{tukey})
    #' 
    #' The cutoff for the Huber and Tukey loss functions are determined from a provided cutoff 
    #' parameter \code{loss_c}. The implementations of these losses are approximations (pseudo-huber and sigmoid 
    #' approxmation respectively) for smooth derivatives.
    #' @param loss_c cutoff value for huber and tukey loss functions. Defaults to \code{c=3}
    #' @param silent logical value whether or not to suppress printed messages such as 'Checking Data',
    #' 'Building Model', etc. Default behaviour (FALSE) is to print the messages.
    construct_nll = function(data,
                             method = "ekf",
                             ode.solver = "rk4",
                             ode.timestep = diff(data$t),
                             loss = "quadratic",
                             loss_c = 3,
                             unscented_hyperpars = list(alpha=1, beta=0, kappa=3-private$number.of.states),
                             compile=FALSE,
                             silent=FALSE){
      
      # set flags
      private$set_compile(compile)
      private$set_method(method)
      private$set_ode_solver(ode.solver)
      private$set_timestep(ode.timestep)
      private$set_simulation_timestep(ode.timestep)
      private$set_loss(loss,loss_c)
      
      # build model
      if(!silent) message("Building model...")
      build_model(self, private)
      
      # check and set data
      if(!silent) message("Checking data...")
      check_and_set_data(data, unscented_hyperpars, self, private)
      
      # construct neg. log-likelihood
      if(!silent) message("Constructing objective function...")
      comptime <- system.time(
        construct_makeADFun(self, private)
      )
      
      comptime = format(round(as.numeric(comptime["elapsed"])*1e4)/1e4,digits=5,scientific=F)
      if(!private$silent) message("...took: ", comptime, " seconds.")
      
      # return
      if(!silent) message("Succesfully returned function handlers")
      return(private$nll)
    },
    
    ########################################################################
    # ESTIMATE FUNCTION
    ########################################################################
    #' @description Estimate the fixed effects parameters in the specified model.
    #' 
    #' @param data data.frame containing time-vector 't', observations and inputs. The observations
    #' can take \code{NA}-values.  
    #' @param use.hessian boolean value. The default (\code{TRUE}) causes the optimization algorithm
    #' \code{stats::nlminb} to use the fixed effects hessian of the (negative log) likelihood when
    #' performing the optimization. This feature is only available for the kalman filter methods 
    #' without any random effects.
    #' @param ode.timestep numeric value. Sets the time step-size in numerical filtering schemes. 
    #' The defined step-size is used to calculate the number of steps between observation time-points as 
    #' defined by the provided \code{data}. If the calculated number of steps is larger than N.01 where N 
    #' is an integer, then the time-step is reduced such that exactly N+1 steps is taken between observations  
    #' The step-size is used in the two following ways depending on the
    #' chosen method:
    #' 1. Kalman filters: The time-step is used as the step-size in the
    #' numerical Forward-Euler scheme to compute the prior state mean and
    #' covariance estimate as the final time solution to the first and second
    #' order moment differential equations.
    #' 2. TMB method: The time-step is used as the step-size in the Euler-Maruyama
    #' scheme for simulating a sample path of the stochastic differential equation,
    #' which serves to link together the latent (random effects) states.
    #' @param ode.solver Sets the ODE solver used in the Kalman Filter methods for solving the moment 
    #' differential equations. The default "euler" is the Forward Euler method, alternatively the classical
    #' 4th order Runge Kutta method is available via "rk4".
    #' @param compile boolean value. The default (\code{FALSE}) is to not compile the C++ objective
    #' function but assume it is already compiled and corresponds to the specified model object. It is
    #' the user's responsibility to ensure correspondence between the specified model and the precompiled
    #' C++ object. If a precompiled C++ object is not found in the specified directory i.e. 
    #' in \code{<cppfile_directory>/<modelname>/(dll/so)} then the compile flag is set to \code{TRUE}.
    #' If the user makes changes to system equations, observation equations, observation variances, 
    #' algebraic relations or lamperi transformations then the C++ object should be recompiled.
    #' @param method character vector - one of either "ekf", "ukf" or "tmb". Sets the estimation 
    #' method. The package has three available methods implemented:
    #' 1. The natural TMB-style formulation where latent states are considered random effects
    #' and are integrated out using the Laplace approximation. This method only yields the gradient
    #' of the (negative log) likelihood function with respect to the fixed effects for optimization.
    #' The method is slower although probably has some precision advantages, and allows for non-Gaussian
    #' observation noise (not yet implemented). One-step / K-step residuals are not yet available in
    #' the package.
    #' 2. (Continous-Discrete) Extended Kalman Filter where the system dynamics are linearized
    #' to handle potential non-linearities. This is computationally the fastest method.
    #' 3. (Continous-Discrete) Unscented Kalman Filter. This is a higher order non-linear Kalman Filter
    #' which improves the mean and covariance estimates when the system display high nonlinearity, and
    #' circumvents the necessity to compute the jacobian of the drift and observation functions.
    #' 
    #' All package features are currently available for the kalman filters, while TMB is limited to
    #' parameter estimation. In particular, it is straight-forward to obtain k-step-ahead predictions
    #' with these methods (use the \code{predict} S3 method), and stochastic simulation is also available 
    #' in the cases where long prediction horizons are sought, where the normality assumption will be 
    #' inaccurate.
    #' @param unscented_hyperpars the three hyper-parameters \code{alpha}, \code{beta} and \code{kappa} defining
    #' the unscented transformation.
    #' @param unconstrained.optim boolean value. When TRUE then the optimization is carried out unconstrained i.e.
    #' without any of the parameter bounds specified during \code{add_parameters}.
    #' @param loss character vector. Sets the loss function type (only implemented for the kalman filter
    #' methods). The loss function is per default quadratic in the one-step residauls as is natural 
    #' when the Gaussian (negative log) likelihood is evaluated, but if the tails of the 
    #' distribution is considered too small i.e. outliers are weighted too much, then one 
    #' can choose loss functions that accounts for this. The three available types available:
    #' 
    #' 1. Quadratic loss (\code{quadratic}).
    #' 2. Quadratic-Linear (\code{huber})
    #' 3. Quadratic-Constant (\code{tukey})
    #' 
    #' The cutoff for the Huber and Tukey loss functions are determined from a provided cutoff 
    #' parameter \code{loss_c}. The implementations of these losses are approximations (pseudo-huber and sigmoid 
    #' approxmation respectively) for smooth derivatives.
    #' @param laplace.residuals boolean - whether or not to calculate one-step ahead residuls
    #' using the method of \link[TMB]{oneStepPredict}.
    #' @param loss_c cutoff value for huber and tukey loss functions. Defaults to \code{c=3}
    #' @param control list of control parameters parsed to \code{nlminb} as its \code{control} argument. 
    #' See \code{?stats::nlminb} for more information
    #' @param silent logical value whether or not to suppress printed messages such as 'Checking Data',
    #' 'Building Model', etc. Default behaviour (FALSE) is to print the messages.
    estimate = function(data, 
                        method = "ekf",
                        ode.solver = "rk4",
                        ode.timestep = diff(data$t),
                        loss = "quadratic",
                        loss_c = 3,
                        unscented_hyperpars = list(alpha=1, beta=0, kappa=3-private$number.of.states),
                        control = list(trace=1,iter.max=1e5,eval.max=1e5),
                        use.hessian = FALSE,
                        laplace.residuals = FALSE,
                        unconstrained.optim = FALSE,
                        compile = FALSE,
                        silent = FALSE){
      
      # settings and flags
      private$set_compile(compile)
      private$set_method(method)
      private$set_ode_solver(ode.solver)
      private$set_timestep(ode.timestep)
      private$set_simulation_timestep(ode.timestep)
      private$use_hessian(use.hessian)
      private$set_loss(loss, loss_c)
      private$set_control(control)
      private$set_unconstrained_optim(unconstrained.optim)
      private$set_silence(silent)
      
      # build model
      if(!private$silent) message("Building model...")
      build_model(self ,private)
      
      # check and set data
      if(!private$silent) message("Checking data...")
      check_and_set_data(data, unscented_hyperpars, self, private)
      
      # construct neg. log-likelihood function
      if(!private$silent) message("Constructing objective function and derivative tables...")
      comptime <- system.time(
      construct_makeADFun(self, private)
      )
      comptime = format(round(as.numeric(comptime["elapsed"])*1e4)/1e4,digits=5,scientific=F)
      if(!private$silent) message("...took: ", comptime, " seconds.")
      
      # estimate
      if(!private$silent) message("Minimizing the negative log-likelihood...")
      optimize_negative_loglikelihood(self, private)
      
      # exit if optimization failed
      if(is.null(private$opt)){
        return(invisible(NULL))
      }
      
      # create return fit
      if(!private$silent) message("Returning results...")
      create_return_fit(self, private, laplace.residuals)
      if(!private$silent) message("Finished!")
      
      # return cloned fit
      cloned.self = self$clone(deep=TRUE)
      cloned.private = cloned.self$.__enclos_env__$private
      return(invisible(cloned.private$fit))
    },
    ########################################################################
    # PREDICT FUNCTION
    ########################################################################
    #' @description Perform prediction/filtration to obtain state mean and covariance estimates. The predictions are
    #' obtained by solving the moment equations \code{n.ahead} steps forward in time when using the current step posterior 
    #' state estimate as the initial condition. 
    #' 
    #' @return A data.frame that contains for each time step the posterior state estimate at that time.step (\code{k = 0}), and the
    #' prior state predictions (\code{k = 1,...,n.ahead}). If \code{return.covariance = TRUE} then the state covariance/correlation 
    #' matrix is returned, otherwise only the marginal variances are returned.
    #' 
    #' @param data data.frame containing time-vector 't', observations and inputs. The observations
    #' can take \code{NA}-values.
    #' @param pars fixed parameter vector parsed to the objective function for prediction/filtration. The default
    #' parameter values used are the initial parameters provided through \code{add_parameters}, unless the \code{estimate}
    #' function has been run, then the default values will be those at the found optimum.
    #' @param k.ahead integer specifying the desired number of time-steps (as determined by the provided
    #' data time-vector) for which predictions are made (integrating the moment ODEs forward in time without 
    #' data updates).
    #' @param return.k.ahead numeric vector of integers specifying which n.ahead predictions to that
    #' should be returned.
    #' @param return.covariance booelan value to indicate whether the covariance (instead of the correlation) 
    #' should be returned.
    #' @param initial.state a named list of two entries 'x0' and 'p0' containing the initial state and covariance of the state
    #' @param ode.timestep numeric value. Sets the time step-size in numerical filtering schemes. 
    #' The defined step-size is used to calculate the number of steps between observation time-points as 
    #' defined by the provided \code{data}. If the calculated number of steps is larger than N.01 where N 
    #' is an integer, then the time-step is reduced such that exactly N+1 steps is taken between observations  
    #' The step-size is used in the two following ways depending on the
    #' chosen method:
    #' 1. Kalman filters: The time-step is used as the step-size in the
    #' numerical Forward-Euler scheme to compute the prior state mean and
    #' covariance estimate as the final time solution to the first and second
    #' order moment differential equations.
    #' 2. TMB method: The time-step is used as the step-size in the Euler-Maruyama
    #' scheme for simulating a sample path of the stochastic differential equation,
    #' which serves to link together the latent (random effects) states.
    #' @param ode.solver Sets the ODE solver used in the Kalman Filter methods for solving the moment 
    #' differential equations. The default "euler" is the Forward Euler method, alternatively the classical
    #' 4th order Runge Kutta method is available via "rk4".
    #' @param method
    #' 1. The natural TMB-style formulation where latent states are considered random effects
    #' and are integrated out using the Laplace approximation. This method only yields the gradient
    #' of the (negative log) likelihood function with respect to the fixed effects for optimization.
    #' The method is slower although probably has some precision advantages, and allows for non-Gaussian
    #' observation noise (not yet implemented). One-step / K-step residuals are not yet available in
    #' the package.
    #' 2. (Continous-Discrete) Extended Kalman Filter where the system dynamics are linearized
    #' to handle potential non-linearities. This is computationally the fastest method.
    #' 3. (Continous-Discrete) Unscented Kalman Filter. This is a higher order non-linear Kalman Filter
    #' which improves the mean and covariance estimates when the system display high nonlinearity, and
    #' circumvents the necessity to compute the jacobian of the drift and observation functions.
    #' 
    #' All package features are currently available for the kalman filters, while TMB is limited to
    #' parameter estimation. In particular, it is straight-forward to obtain k-step-ahead predictions
    #' with these methods (use the \code{predict} S3 method), and stochastic simulation is also available 
    #' in the cases where long prediction horizons are sought, where the normality assumption will be 
    #' inaccurate
    #' @param unscented_hyperpars the three hyper-parameters \code{alpha}, \code{beta} and \code{kappa} defining
    #' the unscented transformation.
    #' @param silent logical value whether or not to suppress printed messages such as 'Checking Data',
    #' 'Building Model', etc. Default behaviour (FALSE) is to print the messages.
    #' @param n.sims number of simulations
    #' @param simulation.timestep timestep used in the euler-maruyama scheme
    #' 
    simulate = function(data,
                        method = "ekf",
                        ode.timestep = diff(data$t),
                        ode.solver = "rk4",
                        pars = NULL,
                        initial.state = self$get_initial_state(),
                        n.sims = 100,
                        simulation.timestep = diff(data$t),
                        k.ahead = 1,
                        return.k.ahead = NULL,
                        unscented_hyperpars = list(alpha=1, beta=0, kappa=3-private$number.of.states),
                        silent = FALSE){
      
      
      
      if(method!="ekf"){ stop("The simulate function is currently only implemented for method = 'ekf'.") }
      
      ###### SET FLAGS #######
      private$set_method(method)
      private$set_ode_solver(ode.solver)
      private$set_timestep(ode.timestep)
      private$set_simulation_timestep(simulation.timestep)
      private$set_pred_initial_state(initial.state)
      private$set_silence(silent)
      
      
      ###### BUILD MODEL #######
      if(!private$silent) message("Building model...")
      build_model(self, private, prediction=TRUE)
      
      ###### CHECK AND SET DATA, PARS, ETC.  #######
      if(!private$silent) message("Checking data...")
      check_and_set_data(data, unscented_hyperpars, self, private)
      private$set_n_ahead_and_last_pred_index(k.ahead)
      pars = set_parameters(pars, silent, self, private)
      
      ###### PERFORM PREDICTION #######
      if(!private$silent) message("Compiling C++ functions...")
      create_rcpp_statespace_functions(self, private)
      
      if(!private$silent) message("Simulating...")
      predict.list = perform_rcpp_ekf_simulation(self, 
                                                 private, 
                                                 pars, 
                                                 private$pred.initial.state, 
                                                 n.sims)
      
      # construct return data.frame
      if(!private$silent) message("Constructing return data.frame...")
      list.out = construct_simulate_rcpp_dataframe(pars,
                                                   predict.list,
                                                   data,
                                                   return.k.ahead,
                                                   n.sims,
                                                   self,
                                                   private)
      
      # return
      if(!private$silent) message("Finished.")
      return(invisible(list.out))
    },
    ########################################################################
    # PREDICT FUNCTION
    ########################################################################
    #' @description Perform prediction/filtration to obtain state mean and covariance estimates. The predictions are
    #' obtained by solving the moment equations \code{n.ahead} steps forward in time when using the current step posterior 
    #' state estimate as the initial condition. 
    #' 
    #' @return A data.frame that contains for each time step the posterior state estimate at that time.step (\code{k = 0}), and the
    #' prior state predictions (\code{k = 1,...,n.ahead}). If \code{return.covariance = TRUE} then the state covariance/correlation 
    #' matrix is returned, otherwise only the marginal variances are returned.
    #' 
    #' @param data data.frame containing time-vector 't', observations and inputs. The observations
    #' can take \code{NA}-values.
    #' @param pars fixed parameter vector parsed to the objective function for prediction/filtration. The default
    #' parameter values used are the initial parameters provided through \code{add_parameters}, unless the \code{estimate}
    #' function has been run, then the default values will be those at the found optimum.
    #' @param k.ahead integer specifying the desired number of time-steps (as determined by the provided
    #' data time-vector) for which predictions are made (integrating the moment ODEs forward in time without 
    #' data updates).
    #' @param return.k.ahead numeric vector of integers specifying which n.ahead predictions to that
    #' should be returned.
    #' @param return.covariance booelan value to indicate whether the covariance (instead of the correlation) 
    #' should be returned.
    #' @param initial.state a named list of two entries 'x0' and 'p0' containing the initial state and covariance of the state
    #' @param ode.timestep numeric value. Sets the time step-size in numerical filtering schemes. 
    #' The defined step-size is used to calculate the number of steps between observation time-points as 
    #' defined by the provided \code{data}. If the calculated number of steps is larger than N.01 where N 
    #' is an integer, then the time-step is reduced such that exactly N+1 steps is taken between observations  
    #' The step-size is used in the two following ways depending on the
    #' chosen method:
    #' 1. Kalman filters: The time-step is used as the step-size in the
    #' numerical Forward-Euler scheme to compute the prior state mean and
    #' covariance estimate as the final time solution to the first and second
    #' order moment differential equations.
    #' 2. TMB method: The time-step is used as the step-size in the Euler-Maruyama
    #' scheme for simulating a sample path of the stochastic differential equation,
    #' which serves to link together the latent (random effects) states.
    #' @param ode.solver Sets the ODE solver used in the Kalman Filter methods for solving the moment 
    #' differential equations. The default "euler" is the Forward Euler method, alternatively the classical
    #' 4th order Runge Kutta method is available via "rk4".
    #' @param method
    #' 1. The natural TMB-style formulation where latent states are considered random effects
    #' and are integrated out using the Laplace approximation. This method only yields the gradient
    #' of the (negative log) likelihood function with respect to the fixed effects for optimization.
    #' The method is slower although probably has some precision advantages, and allows for non-Gaussian
    #' observation noise (not yet implemented). One-step / K-step residuals are not yet available in
    #' the package.
    #' 2. (Continous-Discrete) Extended Kalman Filter where the system dynamics are linearized
    #' to handle potential non-linearities. This is computationally the fastest method.
    #' 3. (Continous-Discrete) Unscented Kalman Filter. This is a higher order non-linear Kalman Filter
    #' which improves the mean and covariance estimates when the system display high nonlinearity, and
    #' circumvents the necessity to compute the jacobian of the drift and observation functions.
    #' 
    #' All package features are currently available for the kalman filters, while TMB is limited to
    #' parameter estimation. In particular, it is straight-forward to obtain k-step-ahead predictions
    #' with these methods (use the \code{predict} S3 method), and stochastic simulation is also available 
    #' in the cases where long prediction horizons are sought, where the normality assumption will be 
    #' inaccurate.
    #' @param unscented_hyperpars the three hyper-parameters \code{alpha}, \code{beta} and \code{kappa} defining
    #' the unscented transformation.
    #' @param silent logical value whether or not to suppress printed messages such as 'Checking Data',
    #' 'Building Model', etc. Default behaviour (FALSE) is to print the messages.
    #' 
    predict = function(data,
                       method = "ekf",
                       ode.timestep = diff(data$t),
                       ode.solver = "rk4",
                       pars = NULL,
                       initial.state = self$get_initial_state(),
                       k.ahead = 1,
                       return.k.ahead = 0:k.ahead,
                       return.covariance = TRUE,
                       unscented_hyperpars = list(alpha=1, beta=0, kappa=3-private$number.of.states),
                       silent = FALSE){
      
      
      
      if(method!="ekf"){ stop("The predict function is currently only implemented for method = 'ekf'.") }
      
      ###### SET FLAGS #######
      private$set_method(method)
      private$set_ode_solver(ode.solver)
      private$set_timestep(ode.timestep)
      private$set_simulation_timestep(ode.timestep)
      private$set_pred_initial_state(initial.state)
      private$set_silence(silent)
      
      ###### BUILD MODEL #######
      if(!private$silent) message("Building model...")
      build_model(self, private, prediction=TRUE)
      
      ###### CHECK AND SET DATA, PARS, ETC.  #######
      if(!private$silent) message("Checking data...")
      check_and_set_data(data, unscented_hyperpars, self, private)
      private$set_n_ahead_and_last_pred_index(k.ahead)
      pars = set_parameters(pars, silent, self, private)
      
      
      ##### COMPILE C++ FUNCTIONS #######
      if(!private$silent) message("Compiling C++ functions...")
      create_rcpp_statespace_functions(self, private)
      
      ##### PERFORM PREDICTION #######
      if(!private$silent) message("Predicting...")
      predict.out = perform_prediction(pars, self, private)
      predict.list = perform_rcpp_ekf_prediction(self, private, pars)
      
      ##### CREATE RETURN DATA.FRAME #######
      if(!private$silent) message("Constructing return data.frame...")
      df.out = construct_predict_rcpp_dataframe(pars,
                                                predict.list,
                                                data,
                                                return.covariance,
                                                return.k.ahead,
                                                self,
                                                private)
      
      ##### RETURN #######
      if(!private$silent) message("Finished.")
      return(invisible(df.out))
    },
    ########################################################################
    # PRINT
    ########################################################################
    #' @description Function to print the model object
    print = function() {
      
      n = length(private$sys.eqs)
      m = length(private$obs.eqs)
      p = length(private$inputs)-1
      ng = max(length(unique(unlist(lapply(private$sys.eqs, function(x) x$diff))))-1,0)
      q = length(private$alg.eqs)
      par = length(private$parameters)
      fixedpars = length(private$fixed.pars)
      freepars = length(private$free.pars)
      
      # If the model is empty
      cat("Stochastic State Space Model:")
      basic.data = c(private$modelname,n,ng,m,p,par)
      row.names = c("Name", "States","Diffusions",
                    "Observations","Inputs",
                    "Parameters")
      mat=data.frame(basic.data,row.names=row.names,fix.empty.names=F)
      print(mat,quote=FALSE)
      
      # STATE EQUATIONS
      if (n>0) {
        cat("\nSystem Equations:\n\n")
        lapply(private$sys.eqs,function(x) cat("\t",deparse1(x$form),"\n"))
      }
      
      # OBS EQUATIONS
      if (m>0) {
        cat("\nObservation Equations:\n\n")
        for (i in 1:length(private$obs.eqs)) {
          bool = private$obs.names[i] %in% private$obsvar.names
          bool2 = is.null(private$obs.var[[i]])
          if (bool & bool2) {
            cat("\t",paste(names(private$obs.eqs)[i],": ",sep=""),deparse1(private$obs.eqs[[i]]$form),"+ e", "\t","e ~ N(0,",paste0(deparse1(private$obs.var[[i]]$rhs),")"),"\n")
          } else {
            cat("\t",paste(names(private$obs.eqs)[i],": ",sep=""),deparse1(private$obs.eqs[[i]]$form),"+ e", "\t","e ~ N(0,?)","\n")
          }
        }
      }
      
      # ALGEBRAICS
      if (q>0) {
        cat("\nAlgebraic Relations:\n\n")
        for(i in 1:length(private$alg.eqs)){
          cat("\t",deparse1(private$alg.eqs[[i]]$form),"\n")
        }
      }
      
      # INPUTS
      if (p>1) {
        cat("\nInputs:\n")
        cat("\t", paste(private$input.names[!private$input.names %in% "t"],collapse=", "))
      }
      
      # PARAMETERS
      if (freepars>0) {
        cat("\n\nFree Parameters:\n")
        cat("\t", paste(names(private$free.pars),collapse=", "))
      }
      
      # FIXED PARAMETERS
      if (fixedpars>0) {
        cat("\n\nFixed Parameters:\n")
        cat("\t", paste(names(private$fixed.pars),collapse=", "))
      }
      
      # return
      return(invisible(self))
    },
    
    ########################################################################
    # SUMMARY
    ########################################################################
    #' @description Summary function for fit
    #' @param correlation boolean value. The default (\code{FALSE}) is to not provide the fixed effects parameter
    #' correlation matrix. 
    summary = function(correlation=FALSE) {
      
      # check if model was estimated
      if (is.null(private$fit)) {
        message("Please estimate your model to get a fit summary.")
        return(invisible(NULL))
      }
      
      # return summary fit
      sumfit = summary(private$fit,correlation=correlation)
      
      return(invisible(sumfit$parameters))
    },
    
    ########################################################################
    # PLOT USING GGPLOT
    ########################################################################
    #' @description Function to print the model object
    #' @param plot.obs a vector to indicate which observations should be plotted for. If multiple
    #' are chosen a list of plots for each observation is returned.
    #' @param pacf logical to indicate whether or not the partial autocorrelations should be returned.
    #' The default is FALSE in which case a histogram is returned instead.
    #' @param extended logical. if TRUE additional information is printed
    #' @param ggtheme ggplot2 theme to use for creating the ggplot.
    plot = function(plot.obs=1, 
                    pacf=FALSE,
                    extended=FALSE,
                    ggtheme=getggplot2theme()){
      
      # check if the estimate method has been run i.e. and private$fit has been created
      if (is.null(private$fit)) {
        message("Please estimate your model in order to plot residuals.")
        return(invisible(NULL))
      }
      
      plotlist = plot(private$fit, 
                      plot.obs=plot.obs,
                      pacf=pacf,
                      extended=extended, 
                      ggtheme=ggtheme)
      
      # return
      return(invisible(plotlist))
    }
  ),
  
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  # Private Methods
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  
  private = list(
    
    # Model
    modelname = character(0),
    cppfile.directory = NULL,
    cppfile.path = character(0),
    cppfile.path.with.method = NULL,
    modelname.with.method = NULL,
    
    # Equations
    sys.eqs = NULL,
    obs.eqs = NULL,
    obs.var = NULL,
    alg.eqs = NULL,
    inputs = NULL,
    parameters = NULL,
    initial.state = NULL,
    pred.initial.state = NULL,
    tmb.initial.state.for.parameters = NULL,
    iobs = NULL,
    
    # after algebraics
    sys.eqs.trans = NULL,
    obs.eqs.trans = NULL,
    obs.var.trans = NULL,
    
    # names
    state.names = NULL,
    obs.names = NULL,
    obsvar.names = NULL,
    input.names = NULL,
    parameter.names = NULL,
    
    # options
    method = NULL,
    use.hessian = NULL,
    state.dep.diff = NULL,
    lamperti = NULL,
    compile = NULL,
    loss = NULL,
    tukey.pars = NULL,
    silent = NULL,
    map = NULL,
    control.nlminb = NULL,
    ode.timestep = NULL,
    ode.timestep.size = NULL,
    ode.timesteps = NULL,
    ode.timesteps.cumsum = NULL,
    simulation.timestep = NULL,
    simulation.timesteps = NULL,
    simulation.timestep.size = NULL,
    ode.solver = NULL,
    unconstrained.optim = NULL,
    
    # hidden
    lock.model = FALSE,
    fixed.pars = NULL,
    free.pars = NULL,
    
    # lengths
    number.of.states = NULL,
    number.of.observations = NULL,
    number.of.diffusions = NULL,
    number.of.pars = NULL,
    number.of.inputs = NULL,
    
    # differentials
    diff.processes = NULL,
    diff.terms = NULL,
    diff.terms.obs = NULL,
    
    # data, nll, opt
    data = NULL,
    nll = NULL,
    opt = NULL,
    sdr = NULL,
    fit = NULL,
    
    # predict
    n.ahead = NULL,
    last.pred.index = NULL,
    
    # Rcpp Prediction
    Rcppfunction_f = NULL,
    Rcppfunction_g = NULL,
    Rcppfunction_dfdx = NULL,
    Rcppfunction_h = NULL,
    Rcppfunction_dhdx = NULL,
    Rcppfunction_hvar = NULL,
    
    # unscented transform
    ukf_alpha = NULL,
    ukf_beta = NULL,
    ukf_kappa = NULL,
    
    ########################################################################
    # ADD TRANSFORMED SYSTEM EQS
    ########################################################################
    add_trans_systems = function(formlist) {
      
      # result = check_system_eqs(form, self, private)
      # private$sys.eqs.trans[[result$name]] = result
      
      form = formlist$form
      rhs = form[[3]]
      lhs = form[[2]]
      name = formlist$name
      
      # extract all variables
      bool = unique(all.vars(rhs)) %in% private$diff.processes
      variables = unique(all.vars(rhs))[!bool]
      
      # create transformed system equation
      private$sys.eqs.trans[[name]] = list(
        name = name,
        form = form,
        rhs = rhs,
        allvars = variables,
        diff = private$sys.eqs[[name]]$diff
      )
      
      # return
      return(invisible(self))
    },
    
    ########################################################################
    # ADD TRANSFORMED OBS EQS
    ########################################################################
    add_trans_observations = function(formlist) {
      
      # result = check_observation_eqs(forms, self, private)
      # private$obs.eqs.trans[[result$name]] = result
      
      form = formlist$form
      name = formlist$name
      
      private$obs.eqs.trans[[name]] = list(
        name = name,
        form = form,
        rhs = form[[3]],
        lhs = form[[2]],
        allvars = all.vars(form[[3]])
      )
      
      # return
      return(invisible(self))
    },
    
    ########################################################################
    # ADD TRANSFORMED OBS VAR EQS
    ########################################################################
    # lamperti transform functions
    add_trans_observation_variances = function(formlist) {
      
      # result = check_observation_variance_eqs(form, self, private)
      # private$obs.var.trans[[result$name]] = result
      
      form = formlist$form
      name = formlist$name
      
      private$obs.var.trans[[name]] = list(
        name = name,
        form = form,
        rhs = form[[3]],
        lhs = form[[2]],
        allvars = all.vars(form[[3]])
      )
      
      # return
      return(invisible(self))
    },
    
    ########################################################################
    # SET COMPILE
    ########################################################################
    set_compile = function(bool) {
      
      # check logical
      if (!is.logical(bool)) {
        stop("You must pass a logical value")
      }
      
      # set flag
      private$compile = bool
      
      # return
      return(invisible(self))
    },
    
    ########################################################################
    # SET SILENT
    ########################################################################
    set_silence = function(bool) {
      
      # check logical
      if (!is.logical(bool)) {
        stop("You must pass a logical value")
      }
      
      # set flag
      private$silent = bool
      
      # return
      return(invisible(self))
    },
    
    ########################################################################
    # SET METHOD
    ########################################################################
    # set method
    set_method = function(method) {
      
      # check string
      if (!(is.character(method))) {
        stop("You must pass a string")
      }
      
      # check if method is available
      available_methods = c("ekf","ukf", "ekf_rtmb","laplace")
      if (!(method %in% available_methods)) {
        stop("That method is not available. Please choose one of:
             1. 'ekf' - Extended Kalman Filter in C++ (Requires Compilation, but faster than 'ekf_rtmb')
             2. 'ekf_rtmb' - Extended Kalman Filter with RTMB (No Compilation)
             3. 'ukf' - Unscented Kalman Filter with C++ (Requires Compilation)
             4. 'laplace' - Laplace Approximation using Random Effects Formulation with RTMB (No Compilation)"
        )
      }
      
      # set flag
      private$method = method
      
      # set file with method flag
      private$modelname.with.method = paste0(private$modelname,sprintf("_%s",private$method))
      private$cppfile.path.with.method = paste0(private$cppfile.path,sprintf("_%s",private$method))
      
      # return
      return(invisible(self))
    },
    ########################################################################
    # SET UNCONSTRAINED OPTIMIZATION
    ########################################################################
    # set predict
    set_unconstrained_optim = function(bool) {
      
      # check string
      if (!(is.logical(bool))) {
        stop("You must pass a logical")
      }
      
      # set flag
      private$unconstrained.optim = bool
      
      # return
      return(invisible(self))
    },
    ########################################################################
    # SET ODE TIME-STEP
    ########################################################################
    set_timestep = function(dt) {
      
      # must be numeric
      if (!is.numeric(dt)) {
        stop("The timestep should be a numeric value.")
      }
      
      private$ode.timestep = dt
    },
    ########################################################################
    # SET SIMULATION TIME-STEP
    ########################################################################
    set_simulation_timestep = function(dt) {
      
      # must be numeric
      if (!is.numeric(dt)) {
        stop("The timestep should be a numeric value.")
      }
      
      private$simulation.timestep = dt
    },
    
    ########################################################################
    # UTILITY FUNCTION: FOR SETTING PREDICTION AHEAD AND LAST PRED INDEX
    ########################################################################
    # SET k step ahead and last pred index for obj$predict
    set_n_ahead_and_last_pred_index = function(n.ahead) {
      
      # check if n.ahead is positive with length 1
      if (!(is.numeric(n.ahead)) | !(length(n.ahead==1)) | !(n.ahead >= 1)) {
        stop("n.ahead must be a non-negative numeric integer")
      }
      
      # Find last prediction index to avoid exciting boundary
      last.pred.index = nrow(private$data) - n.ahead
      if(last.pred.index < 1){
        # message("The provided k.ahead is too large, setting it to the maximum value nrow(data)-1.")
        n.ahead = nrow(private$data) - 1
        last.pred.index = 1
      }
      
      # set values
      private$n.ahead = n.ahead
      private$last.pred.index = nrow(private$data) - n.ahead
      
      # return values
      return(invisible(self))
    },
    
    ########################################################################
    # SET USE OF HESSIAN
    ########################################################################
    # USE HESSIAN FUNCTION
    use_hessian = function(bool) {
      
      # check logical
      if (!is.logical(bool)) {
        stop("The entry must be logical")
      }
      
      # set flag
      private$use.hessian = bool
      
      # return
      return(invisible(self))
    },
    
    ########################################################################
    # SET ODE SOLVER 
    ########################################################################
    set_ode_solver = function(ode.solver){
      
      if(length(ode.solver) != 1){
        stop("Please select only one ODE solver algorithm")
      }
      
      # is numeric
      bool = ode.solver %in% c("euler","rk4")
      if(!bool){
        stop("Please select one of the following ODE solvers:
             1. 'euler' - Forward Euler
             2. 'rk4' - Runge-Kutta 4th Order")
      }
      
      # set flag
      switch(ode.solver,
             euler = { private$ode.solver = 1},
             rk4 = { private$ode.solver = 2 }
      )
      
      # return
      return(invisible(self))
    },
    
    ########################################################################
    # SET INITIAL PREDICTION STATE / COVARIANCE
    ########################################################################
    set_pred_initial_state = function(initial.state) {
      
      if(is.null(initial.state)){
        stop("Please provide an initial state for the mean and covariance")
      }
      
      if (is.null(private$sys.eqs)) {
        stop("Please specify system equations first")
      }
      
      if (!is.list(initial.state)) {
        stop("Please provide a list of length two!")
      }
      
      if (length(initial.state) != 2) {
        stop("Please provide a list of length two")
      }
      
      # unpack list items
      x0 = initial.state[[1]]
      p0 = initial.state[[2]]
      
      if (!is.numeric(x0)) {
        stop("The mean vector is not a numeric")
      }
      
      if (any(is.na(x0))) {
        stop("The mean vector contains NAs.")
      }
      
      if (length(x0)!=length(private$sys.eqs)) {
        stop("The initial state vector should have length ",length(private$sys.eqs))
      }
      
      if (!all(dim(p0)==c(length(private$sys.eqs),length(private$sys.eqs)))) {
        stop("The covariance matrix should be square with dimension ", length(private$sys.eqs))
      }
      
      # convert scalar to matrix
      if(!is.matrix(p0) & is.numeric(p0) & length(p0)==1){
        p0 = p0 * diag(1)
      }
      
      if (!is.numeric(p0)) {
        stop("The covariance matrix is not a numeric")
      }
      
      if (any(is.na(p0))) {
        stop("The covariance matrix contains NAs")
      }
      
      if (any(eigen(p0)$values < 0)){
        stop("The covariance matrix is not positive semi-definite")
      }
      
      if (!isSymmetric.matrix(p0)){
        stop("The covariance matrix is symmetric")
      }
      
      # private$initial.state = list(x0=x0, p0=as.matrix(p0))
      
      # set field
      private$pred.initial.state = list(x0=x0, p0=as.matrix(p0))
      
      # return
      return(invisible(self))
    },
    
    ########################################################################
    # SET LOSS FUNCTION
    ########################################################################
    # SET LOSS FUNCTION
    set_loss = function(loss, loss_c=3) {
      
      # check for string
      if (!(is.character(loss))) {
        stop("You must pass a string")
      }
      
      # check if method is available
      available_losses = c("quadratic","huber","tukey")
      if (!(loss %in% available_losses)) {
        stop("That method is not available. Please choose one of the following:
             1. 'quadratic' - default quadratic loss
             2. 'huber' - quadratic-linear pseudo-huber loss
             3. 'tukey' - quadratic-constant tukey loss")
      }
      
      # IF QUADRATIC:
      loss.flag = 0L
      
      # IF HUBER:
      if(loss=="huber"){
        loss.flag = 2L
      } 
      
      # IF TUKEY: compute tukey parameters by nonlinear optimization
      # we simplify choose tukey parameters that minimize the distance to a
      # sigmoid function
      if (loss=="tukey") {
        
        loss.flag = 1L
        
        # compute tukey approx coefficients
        rtukey = seq(0,100,by=1e-2)
        ctukey = loss_c
        
        # define tukey function
        funtukey = function(r){
          ifelse(r^2 <= ctukey^2,
                 ctukey^2/6 * (1-(1-(r/ctukey)^2)^3),
                 ctukey^2/6
          )
        }
        
        # define least squares loss function between tukey and sigmoid function
        tukeyloss = function(.pars){
          res = sum((funtukey(rtukey) - .pars[4]*(ctsmTMB:::invlogit2(rtukey,a=.pars[1],b=.pars[2])+.pars[3]))^2)
        }
        
        # minimize least squares between tukey and sigmoid
        tukeyopt = stats::nlminb(start=rep(1,4),objective=tukeyloss)
        
        # set the found parameters
        private$tukey.pars = tukeyopt$par
      } 
      
      # set flag
      private$loss = list(loss=loss.flag, c=loss_c)
      
      # return
      return(invisible(self))
    },
    
    ########################################################################
    # SET NLMIN CONTROL OPTIONS
    ########################################################################
    set_control = function(control) {
      
      # is the control a list?
      if (!(is.list(control))) {
        stop("The control argument must be a list. See ?stats::nlminb for control options")
      }
      
      # set flag
      private$control.nlminb = control
      
      # return
      return(invisible(self))
    },
    
    ########################################################################
    # SET UNSCENTED TRANSFORMATION HYPERPARAMAETERS
    ########################################################################
    set_ukf_hyperpars = function(parlist) {
      
      # is the control a list?
      if (!(is.list(parlist))) {
        stop("Please provide a named list containing 'alpha', 'beta' and 'kappa'")
      }
      
      # set parameters
      private$ukf_alpha = parlist[["alpha"]]
      private$ukf_beta = parlist[["beta"]]
      private$ukf_kappa = parlist[["kappa"]]
      
      # return
      return(invisible(self))
    }
  )
)

