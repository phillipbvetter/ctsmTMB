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
#' @examples
#' library(ctsmTMB)
#' model <- ctsmTMB$new()
#' 
#' # adding a single system equations
#' model$addSystem(dx ~ theta * (mu+u-x) * dt + sigma_x*dw)
#' 
#' # adding an observation equation and setting variance
#' model$addObs(y ~ x)
#' model$setVariance(y ~ sigma_y^2)
#' 
#' # add model input
#' model$addInput(u)
#' 
#' # add parameters
#' model$setParameter(
#'   theta   = c(initial = 1, lower=1e-5, upper=50),
#'   mu      = c(initial=1.5, lower=0, upper=5),
#'   sigma_x = c(initial=1, lower=1e-10, upper=30),
#'   sigma_y = 1e-2
#' )
#' 
#' # set the model initial state
#' model$setInitialState(list(1,1e-1))
#' 
#' # extract the likelihood handlers
#' nll <- model$likelihood(data=Ornstein)
#' 
#' # calculate likelihood, gradient and hessian w.r.t parameters
#' nll$fn(nll$par)
#' nll$gr(nll$par)
#' nll$he(nll$par)
#' 
#' # estimate the parameters using an extended kalman filter
#' fit <- model$estimate(data=Ornstein)
#' 
#' # perform moment predictions
#' pred <- model$predict(data=Ornstein)
#' 
#' # perform stochatic simulations
#' sim <- model$simulate(data=Ornstein, n.sims=10)
#' 
#' @name ctsmTMB
#' @export
ctsmTMB = R6::R6Class(
  
  # Class name
  classname = "ctsmTMB",
  
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  # Public Methods
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
      # modelname, directory and path (directory+name)
      private$modelname = "ctsmTMB_model"
      private$cppfile.directory = NULL
      private$cppfile.path = NULL
      private$cppfile.path.with.method = NULL
      private$modelname.with.method = NULL
      
      # estimation, prediction or simulation?
      private$procedure = NULL
      
      # model equations
      private$sys.eqs = NULL
      private$obs.eqs = NULL
      private$obs.var = NULL
      private$alg.eqs = NULL
      private$inputs = list(t=list(name="t",input=quote(t)))
      private$parameters = NULL
      private$initial.state = NULL
      private$pred.initial.state = list(mean=NULL,cov=NULL)
      private$tmb.initial.state = NULL
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
      private$tukey.pars = rep(0,6)
      private$silent = FALSE
      private$map = NULL
      private$control.nlminb = list()
      private$ode.solver = NULL
      private$unconstrained.optim = NULL
      private$estimate.initial = FALSE
      private$initial.variance.scaling = 1
      
      # rebuild
      private$rebuild.model = TRUE
      private$rebuild.ad = TRUE
      private$rebuild.data = TRUE
      private$old.data = list()
      
      # hidden
      private$fixed.pars = list()
      private$free.pars = list()
      private$pars = NULL
      
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
      private$diff.terms.drift = NULL
      
      # data, nll, opt
      private$data = NULL
      private$nll = NULL
      private$opt = NULL
      private$fit = NULL
      private$prediction = NULL
      private$simulation = NULL
      private$filt = NULL
      private$smooth = NULL
      
      # predict
      private$n.ahead = 0
      private$last.pred.index = 0
      
      # rtmb
      private$rtmb.function.strings = NULL
      private$rtmb.function.strings.indexed = NULL
      private$rekf.function.strings = NULL
      private$rtmb.nll.strings = NULL
      private$rcpp.function.strings = NULL
      
      # rcpp functions
      private$rcpp_function_ptr = NULL
      
      # unscented transform
      private$ukf_hyperpars = list()
    },
    
    ########################################################################
    # GET OBJECT PRIVATE FIELDS
    ########################################################################
    #' @description 
    #' Extract the private fields of a ctsmTMB model object. Primarily used for
    #' debugging.
    .private = function(){
      return(invisible(private))
    },
    
    ########################################################################
    # ADD SYSTEMS
    ########################################################################
    #' @description 
    #' Define stochastic differential equation(s) on the form
    #' 
    #' \code{d<state> ~ f(t,<states>, <inputs>) * dt + g(t, <states>, <inputs>) * dw}
    #' 
    #' 
    #' @param form a formula specifying a stochastic differential equation
    #' @param ... additional formulas similar to \code{form} for specifying 
    #' multiple equations at once.
    addSystem = function(form,...) {
      
      # adding a state triggers a model rebuild
      private$rebuild.model <- TRUE
      
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
    #' Define algebraic observation equations on the form
    #' 
    #' \code{<observation> ~ h(t, <states>, <inputs>) + e)}
    #' 
    #' where \code{h} is the observation function, and \code{e} is normally 
    #' distributed noise with zero mean. 
    #' 
    #' This function only specifies the observation name, and its mean through \code{h}.
    #' 
    #' @param form a formula specifying an observation equation
    #' @param ... additional formulas similar to \code{form} for specifying 
    #' multiple equations at once.
    #' @param obsnames character vector specifying the name of the observation. 
    #' This is used when the left-hand side of `form` consists of more than just 
    #' a single variable (of class 'call').
    addObs = function(form,...,obsnames=NULL) {
      
      # Check obsnames
      if(!is.null(obsnames)){
        if(length(c(form,...))!=length(obsnames)){
          stop("You must supply as many observation names as there are observation equations.")
        }
        if(!is.character(obsnames)){
          stop("The observation names in 'obsnames' must be strings.")
        }
      }
      
      # trigger a rebuild
      private$rebuild.model <- TRUE
      
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
        # if it doesnt already exist
        holder <- private$obs.var[[result$name]]
        if(length(holder)==0){
          private$obs.var[[result$name]] = list()
        }
      })
      
      return(invisible(self))
    },
    
    ########################################################################
    # ADD OBSERVATION VARIANCES
    ########################################################################
    #' @description Specify the variance of an observation equation.
    #' 
    #' A defined observation variable \code{y} in e.g. \code{addObs(y ~ 
    #' h(t,<states>,<inputs>)} is perturbed by Gaussian noise with zero mean and 
    #' variance 
    #' to-be specified using \code{setVariance(y ~ p(t,<states>,<inputs>)}. 
    #' We can for instance declare \code{setVariance(y ~ sigma_x^2} 
    #' where \code{sigma_x} is a fixed effect parameter to be declared through 
    #' \code{setParameter}.
    #' 
    #' @param form formula class specifying the observation equation to be added 
    #' to the system.
    #' @param ... additional formulas identical to \code{form} to specify multiple 
    #' observation equations at a time.
    #' 
    setVariance = function(form,...) {
      
      # trigger a rebuild
      private$rebuild.model <- TRUE
      
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
    #' an input variable \code{u} then it is declared using \code{addInput(u)}. 
    #' The input \code{u} must be contained in the data.frame \code{.data} provided 
    #' when calling the \code{estimate} or \code{predict} methods.
    #' 
    #' @param ... variable names that specifies the name of input variables in the defined system.
    #' 
    addInput =  function(...) {
      
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
    #' the specified model, and specify the initial optimizer guess, as well as
    #' lower / upper bounds during optimization. There are two ways to declare parameters:
    #' 
    #' 1. You can declare parameters using formulas i.e. \code{setParameter( 
    #' theta = c(1,0,10), mu = c(0,-10,10) )}. The first value is the initial 
    #' value for the optimizer, the second value is the lower optimization 
    #' bound and the third value is the upper optimization bound. 
    #' 
    #' 2. You can provide a 3-column matrix where rows corresponds to different 
    #' parameters, and the parameter names are provided as rownames of the matrix. 
    #' The columns values corresponds to the description in the vector format above.
    #' 
    #' @param ... a named vector or matrix as described above.
    setParameter = function(...) {
      
      if(nargs() == 0L) stop("setParameter requires at least one parameter vector or matrix")
      
      arglist = list(...)
      argnames = names(arglist)
      
      
      # run over each parameter argument either a vector or a matrix
      for(i in 1:length(arglist)){
        
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
          
          # begin if-statement to seperate fix and free parameters
          if (all(is.na(par.entry[c("lower","upper")]))){
            
            # if the entry is new, then recompile the ad graph
            if(is.null(private$fixed.pars[[par.name]])){
              private$rebuild.ad <- TRUE
            }
            
            # set the parameter values
            private$fixed.pars[[par.name]] = private$parameters[[par.name]]
            private$fixed.pars[[par.name]][["factor"]] = factor(NA)
            
            # remove the parameter from the free list (in case it was there previously)
            private$free.pars[[par.name]] = NULL
          } else {
            
            # if the entry is new, then recompile the ad graph
            if(is.null(private$free.pars[[par.name]])){
              private$rebuild.ad <- TRUE
            }
            
            # set the parameter
            private$free.pars[[par.name]] = private$parameters[[par.name]]
            
            # remove the parameter from the fixed list (in case it was there previously)
            private$fixed.pars[[par.name]] = NULL
          }
          
          # print(private$fixed.pars)
          
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
            private$free.pars[[parname]] = NULL
            if(all(is.na(par.entry[i,c("lower","upper")]))){
              private$fixed.pars[[parname]] = private$parameters[[parname]]
              private$fixed.pars[[parname]][["factor"]] = factor(NA)
            } else {
              private$free.pars[[parname]] = private$parameters[[parname]]
            }
            
            # update parameter names
            private$parameter.names = names(private$parameters)
            return(invisible(self))
          })
          
          ##### ELSE STOP #####  
        } else {
          
          stop("setParameter only expects vectors or matrices")
          
        }
        
      }
      
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
    #' \code{setAlgebraics(theta ~ exp(logtheta))}. All instances of \code{theta} is replaced
    #' by \code{exp(logtheta)} when compiling the C++ function. Note that you must provide values
    #' for \code{logtheta} now instead of \code{theta} when declaring parameters through 
    #' \code{setParameter}
    #' 
    #' @param form algebraic formula
    #' @param ... additional formulas
    setAlgebraics = function(form,...) {
      
      # trigger a rebuild
      private$rebuild.model <- TRUE
      
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
    #' @param initial.state a named list of two entries 'x0' and 'p0' containing the initial state and covariance of the state.
    #' 
    setInitialState = function(initial.state) {
      
      if (is.null(private$sys.eqs)) {
        stop("Please specify system equations first")
      }
      
      if (!is.list(initial.state) || length(initial.state)!=2) {
        stop("Please provide a list of length 2!")
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
      
      names(initial.state) <- c("x0","p0")
      # Store old initial state
      bool <- identical(initial.state, private$initial.state)
      if(!bool){
        private$rebuild.ad <- TRUE
      }
      
      # Store initial state
      private$initial.state = list(x0=x0, p0=as.matrix(p0))
      
      
      return(invisible(self))
    },
    
    ########################################################################
    # SET INITIAL SCALING OF THE COVARIANCE MATRIX
    ########################################################################
    #' @description 
    #' A scalar value that is multiplied onto the estimated
    #' initial state covariance matrix. The scaling is only applied when the
    #' initial state/cov is estimated, not when it is set by the user.
    #' @param scaling a numeric scalar value.
    
    setInitialVarianceScaling = function(scaling){
      
      if(!is.numeric(scaling)){
        stop("The scaling must be a scalar numerical value.")
      }
      
      private$initial.variance.scaling = scaling
      return(invisible(NULL))
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
    setLamperti = function(transforms, states=NULL) {
      
      # trigger a rebuild
      private$rebuild.model <- TRUE
      
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
    #' is created in the directory specified by the \code{setCppfilesDirectory} 
    #' method with name \code{<modelname>.cpp}
    #' 
    #' @param name string defining the model name.
    setModelname = function(name) {
      
      # was a string passed?
      if (!is.character(name)) {
        stop("The modelname must be a string")
      }
      
      # set name field
      private$modelname = name
      
      # update path-field
      # private$cppfile.path <- file.path(private$cppfile.directory, private$modelname)
      
      # return
      return(invisible(self))
    },
    
    ########################################################################
    # SET MAXIMUM A POSTERIORI 
    ########################################################################
    #' @description Enable maximum a posterior (MAP) estimation.
    #'
    #' Adds a maximum a posterior contribution to the (negative log) likelihood 
    #' function by  evaluating the fixed effects parameters in a multivariate Gaussian 
    #' with \code{mean} and \code{covariance} as provided.
    #' 
    #' @param mean mean vector of the Gaussian prior parameter distribution
    #' @param cov covariance matrix of the Gaussian prior parameter distribution
    setMAP = function(mean,cov) {
      
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
    getSystems = function() {
      
      # extract system formulas
      syseqs = lapply(private$sys.eqs,function(x) x$form)
      
      # return
      return(syseqs)
    },
    
    ########################################################################
    # GET OBSERVATIONS
    ########################################################################
    #' @description Retrieve observation equations.
    getObservations = function() {
      
      # extract observation formulas
      obseqs = lapply(private$obs.eqs,function(x) x$form)
      
      # return
      return(obseqs)
    },
    
    ########################################################################
    # GET OBSERVATION VARIANCES
    ########################################################################
    #' @description Retrieve observation variances
    getVariances = function() {
      
      # extract observation variance formulas
      obsvar = lapply(private$obs.var,function(x) x$form)
      
      # return
      return(obsvar)
    },
    
    ########################################################################
    # GET ALGEBRAICS
    ########################################################################
    #' @description Retrieve algebraic relations
    getAlgebraics = function() {
      
      # extract algebraic relation formulas
      algs = lapply(private$alg.eqs,function(x) x$form)
      
      # return
      return(algs)
    },
    
    ########################################################################
    # GET ALGEBRAICS
    ########################################################################
    #' @description Retrieve initially set state and covariance
    getInitialState = function() {
      
      
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
    getParameters = function(type="all", value="all") {
      
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
    # GET ESTIMATION
    ########################################################################
    #' @description Retrieve initially set state and covariance
    getEstimate = function() {
      
      
      # extract algebraic relation formulas
      if(is.null(private$fit)){
        message("There are no estimation results to be exctracted - run 'estimate'.")
        return(invisible(NULL))
      }
      
      # return
      self.clone <- self$clone()
      fit <- private$fit
      fit$private = self.clone$.__enclos_env__$private
      return(invisible(fit))
    },
    
    
    ########################################################################
    # GET NEG LOG LIKE (MakeADFUN)
    ########################################################################
    #' @description Retrieve initially set state and covariance
    getLikelihood = function() {
      
      
      # extract algebraic relation formulas
      if(is.null(private$nll)){
        message("There is no likelihood function to be exctrated - run 'estimate' or 'likelihood'.")
        return(invisible(NULL))
      }
      
      # return
      return(invisible(private$nll))
    },
    
    ########################################################################
    # GET PREDICTION
    ########################################################################
    #' @description Retrieve initially set state and covariance
    getPrediction = function() {
      
      
      # extract algebraic relation formulas
      if(is.null(private$prediction)){
        message("There are no prediction results to be extracted - run 'predict'.")
        return(invisible(NULL))
      }
      
      # return
      return(invisible(private$prediction))
    },
    
    ########################################################################
    # GET SIMULATION
    ########################################################################
    #' @description Retrieve initially set state and covariance
    getSimulation = function() {
      
      
      # extract algebraic relation formulas
      if(is.null(private$simulation)){
        message("There are no simulation results to be extracted - run 'simulate'.")
        return(invisible(NULL))
      }
      
      # return
      return(invisible(private$simulation))
    },
    ########################################################################
    # FILTERING
    ########################################################################
    #' @description Perform state filtering (or smoothing for the 'laplace' method)
    #' 
    #' @param data data.frame containing time-vector 't', observations and inputs. The observations
    #' can take \code{NA}-values.
    #' @param pars fixed parameter vector parsed to the objective function for prediction/filtration. The default
    #' parameter values used are the initial parameters provided through \code{setParameter}, unless the \code{estimate}
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
    #' @param method character vector specifying the filtering method used for state/likelihood calculations. 
    #' Must be one of either "lkf", "ekf", "laplace".
    #' @param estimate.initial.state boolean value. When TRUE the initial state and covariance matrices are
    #' estimated as the stationary solution of the linearized mean and covariance differential equations. When the
    #' system contains time-varying inputs, the first element of these is used.
    #' @param loss character vector. Sets the loss function type (only implemented for the kalman filter
    #' methods). The loss function is per default quadratic in the one-step residuals as is natural 
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
    #' approximation respectively) for smooth derivatives.
    #' @param laplace.residuals boolean - whether or not to calculate one-step ahead residuals
    #' using the method of \link[TMB]{oneStepPredict}.
    #' @param loss_c cutoff value for huber and tukey loss functions. Defaults to \code{c=3}
    #' @param silent logical value whether or not to suppress printed messages such as 'Checking Data',
    #' 'Building Model', etc. Default behaviour (FALSE) is to print the messages.
    #' @param use.cpp a boolean to indicate whether to use C++ to perform calculations
    #' @param ... additional arguments
    filter = function(data,
                      pars = NULL,
                      method = "ekf",
                      ode.solver = "euler",
                      ode.timestep = diff(data$t),
                      loss = "quadratic",
                      loss_c = NULL,
                      laplace.residuals = FALSE,
                      estimate.initial.state = FALSE,
                      use.cpp = FALSE,
                      silent = FALSE,
                      ...){
      
      # set flags
      args = list(
        method = method,
        ode.solver = ode.solver,
        ode.timestep = ode.timestep,
        laplace.residuals = laplace.residuals,
        estimate.initial.state = estimate.initial.state,
        silent = silent
      )
      set_flags("filter", args, self, private)
      
      # build model
      build_model(self, private)
      
      # check and set data
      check_and_set_data(data, self, private)
      private$set_loss(loss, loss_c)
      
      # set parameters
      set_parameters(pars, self, private)
      
      # compile cpp functions if using cpp
      if(use.cpp){
        compile_rcpp_functions(self, private)
      }
      
      # filter
      perform_filtering(self, private, use.cpp)
      
      # create return fit
      create_filter_results(self, private, laplace.residuals)
      
      # return
      if(!private$silent) message("Finished!")
      return(invisible(private$filt))
    },
    
    ########################################################################
    # SMOOTHING
    ########################################################################
    #' @description Perform state filtering (or smoothing for the 'laplace' method)
    #' 
    #' @param data data.frame containing time-vector 't', observations and inputs. The observations
    #' can take \code{NA}-values.  
    #' @param pars fixed parameter vector parsed to the objective function for prediction/filtration. The default
    #' parameter values used are the initial parameters provided through \code{setParameter}, unless the \code{estimate}
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
    #' @param method character vector specifying the filtering method used for state/likelihood calculations. 
    #' Must be one of either "lkf", "ekf", "laplace".
    #' @param estimate.initial.state boolean value. When TRUE the initial state and covariance matrices are
    #' estimated as the stationary solution of the linearized mean and covariance differential equations. When the
    #' system contains time-varying inputs, the first element of these is used.
    #' @param loss character vector. Sets the loss function type (only implemented for the kalman filter
    #' methods). The loss function is per default quadratic in the one-step residuals as is natural 
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
    #' approximation respectively) for smooth derivatives.
    #' @param laplace.residuals boolean - whether or not to calculate one-step ahead residuals
    #' using the method of \link[TMB]{oneStepPredict}.
    #' @param loss_c cutoff value for huber and tukey loss functions. Defaults to \code{c=3}
    #' @param silent logical value whether or not to suppress printed messages such as 'Checking Data',
    #' 'Building Model', etc. Default behaviour (FALSE) is to print the messages.
    #' @param ... additional arguments
    smoother = function(data,
                        pars = NULL,
                        method = "ekf",
                        ode.solver = "euler",
                        ode.timestep = diff(data$t),
                        loss = "quadratic",
                        loss_c = NULL,
                        laplace.residuals = FALSE,
                        estimate.initial.state = FALSE,
                        silent = FALSE,
                        ...){
      
      # set flags
      args = list(
        method = method,
        ode.solver = ode.solver,
        ode.timestep = ode.timestep,
        laplace.residuals = laplace.residuals,
        estimate.initial.state = estimate.initial.state,
        silent = silent
      )
      set_flags("smoother", args, self, private)
      
      # build model
      build_model(self, private)
      
      # check and set data
      check_and_set_data(data, self, private)
      private$set_loss(loss, loss_c)
      
      # set parameters
      set_parameters(pars, self, private)
      
      # construct nll AD function if the method is laplace
      if(any(private$method==c("laplace","laplace2"))){
        construct_makeADFun(self, private)
      }
      
      # smooth
      perform_smoothing(self, private)
      
      # create return fit
      create_smooth_results(self, private, laplace.residuals)
      
      # return
      if(!private$silent) message("Finished!")
      return(invisible(private$smooth))
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
    #' @param method character vector specifying the filtering method used for state/likelihood calculations. 
    #' Must be one of either "lkf", "ekf", "laplace".
    #' @param unconstrained.optim boolean value. When TRUE then the optimization is carried out unconstrained i.e.
    #' without any of the parameter bounds specified during \code{setParameter}.
    #' @param estimate.initial.state boolean value. When TRUE the initial state and covariance matrices are
    #' estimated as the stationary solution of the linearized mean and covariance differential equations. When the
    #' system contains time-varying inputs, the first element of these is used.
    #' @param loss character vector. Sets the loss function type (only implemented for the kalman filter
    #' methods). The loss function is per default quadratic in the one-step residuals as is natural 
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
    #' approximation respectively) for smooth derivatives.
    #' @param laplace.residuals boolean - whether or not to calculate one-step ahead residuals
    #' using the method of \link[TMB]{oneStepPredict}.
    #' @param loss_c cutoff value for huber and tukey loss functions. Defaults to \code{c=3}
    #' @param control list of control parameters parsed to \code{nlminb} as its \code{control} argument. 
    #' See \code{?stats::nlminb} for more information
    #' @param silent logical value whether or not to suppress printed messages such as 'Checking Data',
    #' 'Building Model', etc. Default behaviour (FALSE) is to print the messages.
    #' @param trace 0 or 1, passed to \code{control} to determine whether to print optimization information at each step.
    #' @param ... additional arguments
    estimate = function(data, 
                        method = "ekf",
                        ode.solver = "euler",
                        ode.timestep = diff(data$t),
                        loss = "quadratic",
                        loss_c = NULL,
                        trace = 1,
                        control = list(trace=trace,iter.max=1e5,eval.max=1e5),
                        use.hessian = FALSE,
                        laplace.residuals = FALSE,
                        unconstrained.optim = FALSE,
                        estimate.initial.state = FALSE,
                        silent = FALSE,
                        ...){
      
      # set flags
      args = list(
        method = method,
        ode.solver = ode.solver,
        ode.timestep = ode.timestep,
        control = control,
        use.hessian = use.hessian,
        laplace.residuals = laplace.residuals,
        unconstrained.optim = unconstrained.optim,
        estimate.initial.state = estimate.initial.state,
        silent = silent
        
      )
      set_flags("estimation", args, self, private)
      
      
      # build model
      build_model(self, private)
      
      # check and set data
      check_and_set_data(data, self, private)
      private$set_loss(loss, loss_c)
      
      # construct nll AD function
      construct_makeADFun(self, private)
      
      # estimate
      perform_estimation(self, private)
      
      # exit if optimization failed
      if(is.null(private$opt)){
        return(invisible(NULL))
      }
      
      # create return fit
      create_fit(self, private, laplace.residuals)
      
      # return
      if(!private$silent) message("Finished!")
      return(invisible(private$fit))
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
    #' @param method character vector specifying the filtering method used for state/likelihood calculations. 
    #' Must be one of either "lkf", "ekf", "laplace".
    #' @param loss character vector. Sets the loss function type (only implemented for the kalman filter
    #' methods). The loss function is per default quadratic in the one-step residuals as is natural 
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
    #' approximation respectively) for smooth derivatives.
    #' @param loss_c cutoff value for huber and tukey loss functions. Defaults to \code{c=3}
    #' @param estimate.initial.state boolean value. When TRUE the initial state and covariance matrices are
    #' estimated as the stationary solution of the linearized mean and covariance differential equations. When the
    #' system contains time-varying inputs, the first element of these is used.
    #' @param silent logical value whether or not to suppress printed messages such as 'Checking Data',
    #' 'Building Model', etc. Default behaviour (FALSE) is to print the messages.
    #' @param ... additional arguments
    likelihood = function(data,
                          method = "ekf",
                          ode.solver = "euler",
                          ode.timestep = diff(data$t),
                          loss = "quadratic",
                          loss_c = NULL,
                          estimate.initial.state = FALSE,
                          silent=FALSE,
                          ...){
      
      # set flags
      args = list(
        method = method,
        ode.solver = ode.solver,
        ode.timestep = ode.timestep,
        estimate.initial.state = estimate.initial.state,
        silent = silent
      )
      set_flags("construction", args, self, private)
      
      # build model
      build_model(self, private)
      
      # check and set data
      check_and_set_data(data, self, private)
      private$set_loss(loss, loss_c)
      
      # construct nll AD function
      construct_makeADFun(self, private)
      
      # return
      if(!silent) message("Succesfully returned function handlers")
      return(invisible(private$nll))
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
    #' parameter values used are the initial parameters provided through \code{setParameter}, unless the \code{estimate}
    #' function has been run, then the default values will be those at the found optimum.
    #' @param k.ahead integer specifying the desired number of time-steps (as determined by the provided
    #' data time-vector) for which predictions are made (integrating the moment ODEs forward in time without 
    #' data updates).
    #' @param return.k.ahead numeric vector of integers specifying which n.ahead predictions to that
    #' should be returned.
    #' @param return.covariance boolean value to indicate whether the covariance (instead of the correlation) 
    #' should be returned.
    #' @param estimate.initial.state bool - stationary estimation of initial mean and covariance
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
    #' @param method The prediction method
    #' @param silent logical value whether or not to suppress printed messages such as 'Checking Data',
    #' 'Building Model', etc. Default behaviour (FALSE) is to print the messages.
    #' @param use.cpp a boolean to indicate whether to use C++ to perform calculations
    #' @param ... additional arguments
    predict = function(data,
                       pars = NULL,
                       method = "ekf",
                       ode.solver = "euler",
                       ode.timestep = diff(data$t),
                       k.ahead = nrow(data)-1,
                       return.k.ahead = 0:k.ahead,
                       return.covariance = TRUE,
                       initial.state = self$getInitialState(),
                       estimate.initial.state = private$estimate.initial,
                       use.cpp = FALSE,
                       silent = FALSE,
                       ...){
      
      if(method!="ekf"){ stop("The predict function is currently only implemented for method = 'ekf'.") }
      
      # set flags
      args = list(
        method = method,
        ode.solver = ode.solver,
        ode.timestep = ode.timestep,
        initial.state = initial.state,
        estimate.initial.state = estimate.initial.state,
        silent = silent
      )
      set_flags("prediction", args, self, private)
      
      # build model
      build_model(self, private)
      
      # set data
      check_and_set_data(data, self, private)
      
      # set parameters
      set_parameters(pars, self, private)
      set_k_ahead(k.ahead, self, private)
      
      # predict
      if(use.cpp){
        compile_rcpp_functions(self, private)
        ekf_rcpp_prediction(self, private)
      } else {
        ekf_r_prediction(self, private)
      }
      
      # return
      create_return_prediction(return.covariance, return.k.ahead, use.cpp, self, private)
      
      # return
      if(!private$silent) message("Finished!")
      return(invisible(private$prediction))
    },
    
    ########################################################################
    # SIMULATE FUNCTION
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
    #' parameter values used are the initial parameters provided through \code{setParameter}, unless the \code{estimate}
    #' function has been run, then the default values will be those at the found optimum.
    #' @param k.ahead integer specifying the desired number of time-steps (as determined by the provided
    #' data time-vector) for which predictions are made (integrating the moment ODEs forward in time without 
    #' data updates).
    #' @param return.k.ahead numeric vector of integers specifying which n.ahead predictions to that
    #' should be returned.
    #' @param return.covariance boolean value to indicate whether the covariance (instead of the correlation) 
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
    #' @param estimate.initial.state bool - stationary estimation of initial mean and covariance
    #' @param method
    #' 1. The natural TMB-style formulation where latent states are considered random effects
    #' and are integrated out using the Laplace approximation. This method only yields the gradient
    #' of the (negative log) likelihood function with respect to the fixed effects for optimization.
    #' The method is slower although probably has some precision advantages, and allows for non-Gaussian
    #' observation noise (not yet implemented). One-step / K-step residuals are not yet available in
    #' the package.
    #' 2. (Continuous-Discrete) Extended Kalman Filter where the system dynamics are linearized
    #' to handle potential non-linearities. This is computationally the fastest method.
    #' 3. (Continuous-Discrete) Unscented Kalman Filter. This is a higher order non-linear Kalman Filter
    #' which improves the mean and covariance estimates when the system display high nonlinearity, and
    #' circumvents the necessity to compute the Jacobian of the drift and observation functions.
    #' 
    #' All package features are currently available for the kalman filters, while TMB is limited to
    #' parameter estimation. In particular, it is straight-forward to obtain k-step-ahead predictions
    #' with these methods (use the \code{predict} S3 method), and stochastic simulation is also available 
    #' in the cases where long prediction horizons are sought, where the normality assumption will be 
    #' inaccurate
    #' @param silent logical value whether or not to suppress printed messages such as 'Checking Data',
    #' 'Building Model', etc. Default behaviour (FALSE) is to print the messages.
    #' @param n.sims number of simulations
    #' @param simulation.timestep timestep used in the euler-maruyama scheme
    #' @param use.cpp a boolean to indicate whether to use C++ to perform calculations
    #' @param cpp.seed an integer seed value to control RNG normal draws on the C++ side.
    #' @param ... additional arguments
    simulate = function(data,
                        pars = NULL,
                        use.cpp = FALSE,
                        cpp.seed = NULL,
                        method = "ekf",
                        ode.solver = "rk4",
                        ode.timestep = diff(data$t),
                        simulation.timestep = diff(data$t),
                        k.ahead = nrow(data)-1,
                        return.k.ahead = 0:k.ahead,
                        n.sims = 100,
                        initial.state = self$getInitialState(),
                        estimate.initial.state = private$estimate.initial,
                        silent = FALSE,
                        ...){
      
      if(method!="ekf"){ stop("The simulate function is currently only implemented for method = 'ekf'.") }
      
      # set flags
      args = list(
        method = method,
        ode.solver = ode.solver,
        ode.timestep = ode.timestep,
        simulation.timestep = simulation.timestep,
        initial.state = initial.state,
        estimate.initial.state = estimate.initial.state,
        silent = silent,
        cpp.seed = cpp.seed
      )
      set_flags("simulation", args, self, private)
      
      # build model
      build_model(self, private)
      
      # check data
      check_and_set_data(data, self, private)
      
      # set parameters
      set_k_ahead(k.ahead, self, private)
      set_parameters(pars, self, private)
      
      # simulate
      if(use.cpp){
        compile_rcpp_functions(self, private)
        rcpp_simulation(self, private, n.sims)
      } else {
        ekf_r_simulation(self, private, n.sims)
      }
      
      # construct return data.frame
      create_return_simulation(return.k.ahead, n.sims, self, private)
      
      # return
      if(!private$silent) message("Finished.")
      return(invisible(private$simulation))
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
      cat("ctsmTMB model object:")
      # basic.data = c(private$modelname,n,ng,m,p,par)
      # row.names = c("Name", "States","Diffusions",
      #               "Observations","Inputs",
      #               "Parameters")
      basic.data = c(n,ng,m,p,par)
      row.names = c("States","Diffusions",
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
          bool2 = !is.null(private$obs.var[[i]])
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
    }
  ),
  
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  # Private Methods
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  
  private = list(
    
    # model stats
    modelname = character(0),
    cppfile.directory = NULL,
    cppfile.path = character(0),
    cppfile.path.with.method = NULL,
    modelname.with.method = NULL,
    
    # estimation, prediction or simulation?
    procedure = NULL,
    
    # model equations
    sys.eqs = NULL,
    obs.eqs = NULL,
    obs.var = NULL,
    alg.eqs = NULL,
    inputs = NULL,
    parameters = NULL,
    initial.state = NULL,
    pred.initial.state = NULL,
    tmb.initial.state = NULL,
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
    estimate.initial = NULL,
    initial.variance.scaling = NULL,
    
    # rebuild
    rebuild.model = FALSE,
    rebuild.ad = FALSE,
    rebuild.data = FALSE,
    old.data = list(),
    
    # hidden
    fixed.pars = NULL,
    free.pars = NULL,
    pars = NULL,
    
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
    diff.terms.drift = NULL,
    
    # data, nll, opt
    data = NULL,
    nll = NULL,
    opt = NULL,
    sdr = NULL,
    fit = NULL,
    prediction = NULL,
    simulation = NULL,
    filt = NULL,
    smooth = NULL,
    
    # predict
    n.ahead = NULL,
    last.pred.index = NULL,
    
    # function strings
    rtmb.function.strings = NULL,
    rtmb.function.strings.indexed = NULL,
    rekf.function.strings = NULL,
    rtmb.nll.strings = NULL,
    rcpp.function.strings = NULL,
    
    # rcpp
    rcpp_function_ptr = NULL,
    
    # unscented transform
    ukf_hyperpars = NULL,
    
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
        diff.dt = ctsmTMB.Deriv(f=rhs, x="dt"),
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
    set_procedure = function(str) {
      
      # check logical
      if (!is.character(str)) {
        stop("The procedure must be a string - estimation / prediction / simulation")
      }
      
      # set flag
      switch(str,
             filter = {private$procedure <- "filter"},
             smoother = {private$procedure <- "smoother"},
             estimation = {private$procedure <- "estimation"},
             prediction = {private$procedure <- "prediction"},
             simulation = {private$procedure <- "simulation"},
             construction = {private$procedure <- "construction"}
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
      available_methods = c("lkf", "ekf", "ukf" , "laplace", "laplace2")
      if (!(method %in% available_methods)) {
        stop("That method is not available. Please choose one of:
             1. 'lkf' - Linear Kalman Filter
             2. 'ekf' - Extended Kalman Filter
             3. 'ukf' - Unscented Kalman Filter
             4. 'laplace' - Laplace Approximation using Random Effects Formulation (X),
             5. 'laplace2' - Laplace Approximation using Random Effects Formulation (XdB)"
        )
      }
      
      # set flag
      private$method = method
      
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
    #       set_ode_solver = function(ode.solver){
    #         
    #         if(any(private$method == c("lkf","laplace"))){
    #           return(invisible(self))
    #         }
    #         
    #         if(private$method=="ekf"){
    #           available.ode.solvers <- c("euler", 
    #                                      "rk4", 
    #                                      "lsoda", 
    #                                      "lsode", 
    #                                      "lsodes", 
    #                                      "lsodar", 
    #                                      "vode", 
    #                                      "daspk",
    #                                      "ode23", 
    #                                      "ode45", 
    #                                      "radau", 
    #                                      "bdf", 
    #                                      "bdf_d", 
    #                                      "adams", 
    #                                      "impAdams", 
    #                                      "impAdams_d")
    #         } else {
    #           available.ode.solvers <- c("euler","rk4")
    #         }
    #         
    #         # is numeric
    #         bool = ode.solver %in% available.ode.solvers
    #         if(!bool){
    #           stop("You must choose one of the following ode solvers:\n",
    #                paste(available.ode.solvers,collase=" "))
    #         }
    #         
    #         # If using an RTMBode solver check if RTMBode is available
    #         RTMBode.solvers <- c("lsoda", 
    #                              "lsode", 
    #                              "lsodes", 
    #                              "lsodar", 
    #                              "vode", 
    #                              "daspk",
    #                              "ode23", 
    #                              "ode45", 
    #                              "radau", 
    #                              "bdf", 
    #                              "bdf_d", 
    #                              "adams", 
    #                              "impAdams", 
    #                              "impAdams_d")
    #         bool = ode.solver %in% RTMBode.solvers
    #         if(bool){
    #           check.for.package <- requireNamespace("RTMBode", quietly=TRUE)
    #           if(!check.for.package){
    #             stop("The RTMBode package is not installed. Please install the package with:
    #       install.packages('RTMBode', repos = c('https://kaskr.r-universe.dev', 'https://cloud.r-project.org'))
    # or visit https://github.com/kaskr/RTMB for more information."
    #             )
    #           }
    #         }
    #         
    #         
    #         # set flag
    #         private$ode.solver <- switch(ode.solver,
    #                                      euler = 1,
    #                                      rk4 = 2,
    #                                      # otherwise
    #                                      ode.solver
    #         )
    #         
    #         # return
    #         return(invisible(self))
    #       },
    
    ########################################################################
    # SET ODE SOLVER 
    ########################################################################
    set_ode_solver = function(ode.solver){
      
      # these meethods dont use ode solvers
      if(any(private$method == c("lkf","laplace","laplace2"))){
        return(invisible(self))
      }
      
      # check input
      available.ode.solvers <- c("euler","rk4")
      bool = ode.solver %in% available.ode.solvers
      if(!bool){
        stop("You must choose one of the following ode solvers:\n\t",
             paste(available.ode.solvers, collapse = ", "))
      }
      
      # set solver
      private$ode.solver <- switch(ode.solver,
                                   euler = 1,
                                   rk4 = 2)
      
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
      
      # set field
      private$pred.initial.state = list(x0=x0, p0=as.matrix(p0))
      
      # return
      return(invisible(self))
    },
    
    ########################################################################
    # SET LOSS FUNCTION
    ########################################################################
    # SET LOSS FUNCTION
    set_loss = function(loss, loss_c) {
      
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
      
      if(is.null(loss_c)){
        loss_c <- qchisq(0.95, df=private$number.of.observations)
      }
      
      if(loss_c <= 0){
        stop("The loss threshold must be positive")
      }
      
      # set flag
      private$loss = list(loss=loss, c=loss_c)
      
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
    # SET INITIAL STATE ESTIMATION
    ########################################################################
    set_initial_state_estimation = function(bool){
      
      if(!is.logical(bool)){
        stop("The initial state estimation must be TRUE or FALSE.")
      }
      
      private$estimate.initial = bool
      return(invisible(NULL))
    },
    
    ########################################################################
    # SET UNSCENTED TRANSFORMATION HYPERPARAMAETERS
    ########################################################################
    set_ukf_hyperpars = function(parlist) {
      
      # is the control a list?
      if (!(is.list(parlist))) {
        stop("Please provide a named list containing 'alpha', 'beta' and 'kappa'")
      }
      
      # check if entries are numerics
      if (!is.numeric(parlist[["alpha"]])){
        stop("'alpha' must be a numeric")
      }
      if (!is.numeric(parlist[["beta"]])){
        stop("'beta' must be a numeric")
      }
      if (!is.numeric(parlist[["kappa"]])){
        stop("'kappa' must be a numeric")
      }
      
      # set parameters
      private$ukf_hyperpars = do.call(c,unname(parlist))
      
      # return
      return(invisible(self))
    },
    
    ########################################################################
    # SET CPP SEED FOR SIMULATIONS
    ########################################################################
    set_cpp_seed = function(seed) {
      
      # return if null (unset)
      if(is.null(seed)){
        return(invisible(self))
      }
      
      # check for numeric
      if(!is.numeric(seed)){
        stop("The cpp.seed should be a scalar numeric value")
      }
      
      # take first entry if longer than 1
      if(!length(seed) < 1){
        seed <- seed[[1]]
      }
      
      # set the seed with Rcpp function
      ziggsetseed(seed)
      
      # return
      return(invisible(self))
    }
    
  )
)


