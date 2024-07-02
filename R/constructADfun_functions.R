# Functions for construction likelihood functions using TMB / RTMB

#######################################################
# MAIN CONSTRUCT MAKEADFUN FUNCTION THAT CALL OTHERS
#######################################################

construct_makeADFun = function(self, private){
  
  # TMB::openmp(max=TRUE, autopar=TRUE, DLL=private$modelname.with.method)
  
  if(private$method == "ekf"){
    construct_kalman_makeADFun(self, private)
  }
  if(private$method == "ukf"){
    construct_kalman_makeADFun(self, private)
  }
  if(private$method == "ekf_rtmb"){
    construct_rtmb_ekf_makeADFun(self, private)
  }
  
  if(private$method=="laplace"){
    construct_rtmb_laplace_makeADFun(self, private)
  }
  
  return(invisible(self))
}

#######################################################
# CONSTRUCT KALMAN MAKEADFUN
#######################################################
construct_kalman_makeADFun = function(self, private){
  
  ################################################
  # Data
  ################################################
  
  # add mandatory entries to data
  tmb.data = list(
    
    # methods and purpose
    estimation_method = switch(private$method, ekf = 1, ukf = 2),
    ode_solver = private$ode.solver,
    
    # initial
    stateVec = private$initial.state$x0,
    covMat = private$initial.state$p0,
    
    # time-steps
    ode_timestep_size = private$ode.timestep.size,
    ode_timesteps = private$ode.timesteps,
    ode_cumsum_timesteps = private$ode.timesteps.cumsum,
    
    # loss function
    loss_function = private$loss$loss,
    loss_threshold_value = private$loss$c,
    tukey_loss_parameters = private$tukey.pars,
    
    # system size
    number_of_state_eqs = private$number.of.states,
    number_of_obs_eqs = private$number.of.observations,
    number_of_diffusions = private$number.of.diffusions,
    
    # inputs
    inputMat = as.matrix(private$data[private$input.names]),
    
    # observations
    obsMat = as.matrix(private$data[private$obs.names])
  )
  
  # unscented parameters
  ukf_hyperpars_list = list()
  if(private$method=="ukf")
  {
    ukf_hyperpars_list = list(
      ukf_alpha = private$ukf_alpha,
      ukf_beta = private$ukf_beta,
      ukf_kappa = private$ukf_kappa
    )
  }
  
  # MAP Estimation?
  tmb.map.data = list(
    MAP_bool = 0L
  )
  if (!is.null(private$map)) {
    bool = self$get_parameters()[,"type"] == "free"
    tmb.map.data = list(
      MAP_bool = 1L,
      map_mean__ = private$map$mean[bool],
      map_cov__ = private$map$cov[bool,bool],
      map_ints__ = as.numeric(bool),
      sum_map_ints__ = sum(as.numeric(bool))
    )
  }
  
  # construct final data list
  data = c(tmb.data, private$iobs, tmb.map.data, ukf_hyperpars_list)
  
  ################################################
  # Parameters
  ################################################
  
  parameters = lapply(private$parameters, function(x) x[["initial"]]) # Initial parameter values
  
  ################################################
  # Construct Neg. Log-Likelihood
  ################################################
  
  nll = TMB::MakeADFun(data = data,
                       parameters = parameters,
                       map = lapply(private$fixed.pars, function(x) x$factor),
                       DLL = private$modelname.with.method,
                       silent = TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
}

#######################################################
# CONSTRUCT KALMAN MAKEADFUN WITH RTMB
#######################################################

construct_rtmb_ekf_makeADFun = function(self, private)
{
  
  ################################################
  # Data
  ################################################
  
  # methods and purpose
  ode_solver = private$ode.solver
  
  # initial
  stateVec = cbind(private$initial.state$x0)
  covMat = private$initial.state$p0
  
  # loss function
  loss_function = private$loss$loss
  loss_threshold_value = private$loss$c
  tukey_loss_parameters = private$tukey.pars
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  ode_cumsum_timesteps = private$ode.timesteps.cumsum
  
  # system size
  number_of_state_eqs = private$number.of.states
  number_of_obs_eqs = private$number.of.observations
  number_of_diffusions = private$number.of.diffusions
  
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  # unscented parameters
  # ukf_hyperpars_list = list()
  # if(private$method=="ukf")
  # {
  #   ukf_hyperpars_list = list(
  #     ukf_alpha = private$ukf_alpha,
  #     ukf_beta = private$ukf_beta,
  #     ukf_kappa = private$ukf_kappa
  #   )
  # }
  
  # MAP Estimation?
  MAP_bool = 0L
  if (!is.null(private$map)) {
    bool = self$get_parameters()[,"type"] == "free"
    MAP_bool = 1L
    map_mean__ = private$map$mean[bool]
    map_cov__ = private$map$cov[bool,bool]
    map_ints__ = as.numeric(bool)
    sum_map_ints__ = sum(as.numeric(bool))
  }
  
  ################################################
  # Initial Parameters
  ################################################
  
  parameters = lapply(private$parameters, function(x) x[["initial"]])
  
  ################################################
  # Define general functions
  ################################################
  
  # Log-Determinant Hack
  logdet <- RTMB::ADjoint(
    function(x) {
      dim(x) <- rep(sqrt(length(x)), 2)
      determinant(x, log=TRUE)$modulus
    },
    function(x, y, dy) {
      dim(x) <- rep(sqrt(length(x)), 2)
      t(RTMB::solve(x)) * dy
    },
    name = "logdet")
  
  # ODE Solver
  ode_integrator = function(covMat, stateVec, parVec, inputVec, dinputVec, dt, ode_solver){
    
    # Initials
    X0 = stateVec
    P0 = covMat
    
    # Explicit Forward Euler
    if(ode_solver==1){
      X1 = X0 + f__(stateVec, parVec, inputVec) * dt
      P1 = P0 + cov_ode_1step(covMat, stateVec, parVec, inputVec) * dt
    }
    
    # Classical 4th Order Runge-Kutta Method
    if(ode_solver==2){
      
      # 1. Approx Slope at Initial Point
      k1 = f__(stateVec, parVec, inputVec)
      c1 = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # 2. First Approx Slope at Midpoint
      # inputVec = inputVec + 0.5 * dinputVec
      stateVec = X0 + 0.5 * dt * k1
      covMat   = P0 + 0.5 * dt * c1
      k2       = f__(stateVec, parVec, inputVec)
      c2       = cov_ode_1step(covMat, stateVec, parVec, inputVec)   
      
      # 3. Second Approx Slope at Midpoint
      stateVec = X0 + 0.5 * dt * k2
      covMat   = P0 + 0.5 * dt * c2
      k3       = f__(stateVec, parVec, inputVec)
      c3       = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # 4. Approx Slope at End Point
      # inputVec = inputVec + 0.5 * dinputVec
      stateVec = X0 + dt * k3
      covMat   = P0 + dt * c3
      k4       = f__(stateVec, parVec, inputVec)
      c4       = cov_ode_1step(covMat, stateVec, parVec, inputVec)
      
      # ODE UPDATE
      X1 = X0 + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0 * dt
      P1 = P0 + (c1 + 2.0*c2 + 2.0*c3 + c4)/6.0 * dt
    }
    return(invisible(list(X1,P1)))
  }
  
  # Covariance ODE 1-Step
  cov_ode_1step = function(covMat, stateVec, parVec, inputVec){
    A = dfdx__(stateVec, parVec, inputVec)
    G = g__(stateVec, parVec, inputVec)
    AcovMat = A %*% covMat
    return(AcovMat + t(AcovMat) + G %*% t(G))
  }
  
  
  ################################################
  # Define general functions
  ################################################
  
  # Get function elements in list
  elements = get_rtmb_function_elements(self, private)
  
  sde.functions.txt = paste('###### DEFINE FUNCTIONS #######
# drift function
f__ = function(stateVec, parVec, inputVec){
ans = c(F_ELEMENTS)
return(ans)
}

# jacobian drift function
dfdx__ = function(stateVec, parVec, inputVec){
ans = RTMB::matrix(c(DFDX_ELEMENTS), nrow=NUMBER_OF_STATES, ncol=NUMBER_OF_STATES, byrow=T)
return(ans)
}

# diffusion function
g__ = function(stateVec, parVec, inputVec){
ans = RTMB::matrix(c(G_ELEMENTS), nrow=NUMBER_OF_STATES, ncol=NUMBER_OF_DIFFUSIONS, byrow=T)
return(ans)
}

# obs function
h__ = function(stateVec, parVec, inputVec){
ans = c(H_ELEMENTS)
return(ans)
}

# jacobian obs function
dhdx__ = function(stateVec, parVec, inputVec){
ans = RTMB::matrix(c(DHDX_ELEMENTS),nrow=NUMBER_OF_OBSERVATIONS, ncol=NUMBER_OF_STATES, byrow=T)
return(ans)
}

# variance
hvar__ = function(stateVec, parVec, inputVec){
ans = RTMB::matrix(c(HVAR_ELEMENTS),nrow=NUMBER_OF_OBSERVATIONS, ncol=NUMBER_OF_OBSERVATIONS, byrow=T)
return(ans)
}')
  
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="NUMBER_OF_STATES", 
                                               replacement=deparse(as.numeric(private$number.of.states)))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="NUMBER_OF_DIFFUSIONS", 
                                               replacement=deparse(as.numeric(private$number.of.diffusions)))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="NUMBER_OF_OBSERVATIONS", 
                                               replacement=deparse(as.numeric(private$number.of.observations)))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="F_ELEMENTS", 
                                               replacement=paste(elements$f,collapse=","))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="DFDX_ELEMENTS", 
                                               replacement=paste(elements$dfdx,collapse=","))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="G_ELEMENTS", 
                                               replacement=paste(elements$g,collapse=","))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="H_ELEMENTS", 
                                               replacement=paste(elements$h,collapse=","))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="DHDX_ELEMENTS", 
                                               replacement=paste(elements$dhdx,collapse=","))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="HVAR_ELEMENTS", 
                                               replacement=paste(elements$hvar,collapse=","))
  eval(parse(text=sde.functions.txt))
  
  ################################################
  # Define main likelihood function
  ################################################
  main.function.txt = paste('ekf.nll = function(p){
####### INITIALIZATION #######
# storage variables
xPrior <- pPrior <- xPost <- pPost <- Innovation <- InnovationCovariance <- vector("list",length=nrow(obsMat))
xPrior[[1]] <- xPost[[1]] <- stateVec
pPrior[[1]] <- pPost[[1]] <- covMat

# constants
halflog2pi = log(2*pi)/2
I0 <- diag(NUMBER_OF_STATES)
E0 <- V0 <- diag(NUMBER_OF_OBSERVATIONS)

# set negative log-likelihood
nll = 0

# fixed effects parameter vector
parVec = c(FIXED_PARAMETERS)

# ###### TIME LOOP START #######
for(i in 1:(nrow(obsMat)-1)){

# # Define inputs and use first order input interpolation
inputVec = inputMat[i,]
dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]

###### TIME UPDATE - ODE SOLVER #######
for(j in 1:ode_timesteps[i]){
sol = ode_integrator(covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size[i], ode_solver)
stateVec = sol[[1]]
covMat = sol[[2]]
}
xPrior[[i+1]] = stateVec
pPrior[[i+1]] = covMat

######## DATA UPDATE - KALMAN ALGORITHM ########
obsVec = obsMat[i+1,]
inputVec = inputMat[i+1,]
obsVec_bool = !is.na(obsVec)

######## DATA UPDATE - IF ANY DATA AVAILABLE ######## 
if(any(obsVec_bool)){
s = sum(obsVec_bool)
y = obsVec[obsVec_bool]
E = E0[obsVec_bool,obsVec_bool,drop=FALSE]
H = h__(stateVec, parVec, inputVec)
C = E %*% dhdx__(stateVec, parVec, inputVec)
e = y - E %*% H
V0 = hvar__(stateVec, parVec, inputVec)
V = E %*% V0 %*% t(E)
R = C %*% covMat %*% t(C) + V
Ri = RTMB::solve(R)
K = covMat %*% t(C) %*% Ri

# Update State
stateVec = stateVec + K %*% e
covMat = (I0 - K %*% C) %*% covMat %*% t(I0 - K %*% C) + K %*% V %*% t(K)
nll = nll + 0.5 * logdet(R) + 0.5 * t(e) %*% Ri %*% e + halflog2pi * s

# Store innovation and covariance
Innovation[[i+1]] = e
InnovationCovariance[[i+1]] = R
}
xPost[[i+1]] = stateVec
pPost[[i+1]] = covMat
}

# ###### MAXIMUM A POSTERIOR #######

# ###### REPORT #######
RTMB::REPORT(Innovation)
RTMB::REPORT(InnovationCovariance)
RTMB::REPORT(xPrior)
RTMB::REPORT(xPost)
RTMB::REPORT(pPrior)
RTMB::REPORT(pPost)

# ###### RETURN #######
return(nll)
}')
  
  main.function.txt = stringr::str_replace_all(main.function.txt, 
                                               pattern="NUMBER_OF_STATES", 
                                               replacement=deparse(as.numeric(private$number.of.states)))
  main.function.txt = stringr::str_replace_all(main.function.txt, 
                                               pattern="NUMBER_OF_DIFFUSIONS", 
                                               replacement=deparse(as.numeric(private$number.of.diffusions)))
  main.function.txt = stringr::str_replace_all(main.function.txt, 
                                               pattern="NUMBER_OF_OBSERVATIONS", 
                                               replacement=deparse(as.numeric(private$number.of.observations)))
  main.function.txt = stringr::str_replace_all(main.function.txt, 
                                               pattern="STATE_NAMES", 
                                               replacement=paste("p$",private$state.names,sep="",collapse=","))
  main.function.txt = stringr::str_replace_all(main.function.txt, 
                                               pattern="FIXED_PARAMETERS", 
                                               replacement=paste("p$",private$parameter.names,sep="",collapse=","))
  eval(parse(text=main.function.txt))
  
  ################################################
  # Construct Neg. Log-Likelihood
  ################################################
  # private$nll = ekf.nll
  # stop("halting")
  
  nll = RTMB::MakeADFun(func = ekf.nll,
                        parameters=parameters,
                        map = lapply(private$fixed.pars, function(x) x$factor),
                        silent=TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
  
}

#######################################################
# CONSTRUCT LAPLACE MAKEADFUN WITH RTMB
#######################################################

construct_rtmb_laplace_makeADFun = function(self, private)
{
  
  ################################################
  # Data
  ################################################
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  ode_cumsum_timesteps = private$ode.timesteps.cumsum
  
  # system size
  number_of_state_eqs = private$number.of.states
  number_of_obs_eqs = private$number.of.observations
  number_of_diffusions = private$number.of.diffusions
  
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  # indices in state parameter vectors corresponding to indices in observations / inputs
  # add 1 because too lazy to change private$iobs from 0-index to 1-indexed.
  iobs = lapply(private$iobs,function(x) x+1)
  
  ################################################
  # Parameters
  ################################################
  
  parameters = c(
    lapply(private$parameters, function(x) x[["initial"]]),
    private$tmb.initial.state.for.parameters
  )
  
  ################################################
  # Define function
  ################################################
  
  elements = get_rtmb_function_elements(self, private)
  
  sde.functions.txt = paste('###### DEFINE FUNCTIONS #######
# drift function
f__ = function(stateVec, parVec, inputVec){
ans = c(F_ELEMENTS)
return(ans)
}

# diffusion function
g__ = function(stateVec, parVec, inputVec){
ans = RTMB::matrix(c(G_ELEMENTS), nrow=NUMBER_OF_STATES, ncol=NUMBER_OF_DIFFUSIONS, byrow=T)
return(ans)
}

h__ = function(stateVec, parVec, inputVec){
ans = c(H_ELEMENTS)
return(ans)
}

hvar__ = function(stateVec, parVec, inputVec){
ans = c(HVAR_ELEMENTS)
return(ans)
}')
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="NUMBER_OF_STATES", 
                                               replacement=deparse(as.numeric(private$number.of.states)))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="NUMBER_OF_DIFFUSIONS", 
                                               replacement=deparse(as.numeric(private$number.of.diffusions)))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="NUMBER_OF_OBSERVATIONS", 
                                               replacement=deparse(as.numeric(private$number.of.observations)))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="F_ELEMENTS", 
                                               replacement=paste(elements$f,collapse=","))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="G_ELEMENTS", 
                                               replacement=paste(elements$g,collapse=","))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="H_ELEMENTS", 
                                               replacement=paste(elements$h,collapse=","))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="HVAR_ELEMENTS", 
                                               replacement=paste(elements$hvar,collapse=","))
  eval(parse(text=sde.functions.txt))
  
  main.function.txt = paste('laplace.nll = function(p){
# set negative log-likelihood
nll = 0

# small identity matrix
small_identity = diag(1e-8, nrow=NUMBER_OF_STATES, ncol=NUMBER_OF_STATES)

# extract state random effects
stateMat = cbind(STATE_NAMES)

# fixed effects parameter vector
parVec = cbind(FIXED_PARAMETERS)

###### TIME LOOP START #######
for(i in 1:(nrow(obsMat)-1)) {

# Define inputs and use first order input interpolation
inputVec = inputMat[i,]
dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]

###### BETWEEN TIME POINTS LOOP START #######
for(j in 1:ode_timesteps[i]){

# grab current and next state
x_now = stateMat[ode_cumsum_timesteps[i]+j,]
x_next = stateMat[ode_cumsum_timesteps[i]+j+1,]

# compute drift (vector) and diffusion (matrix)
f = f__(x_now, parVec, inputVec)
g = g__(x_now, parVec, inputVec)
inputVec = inputVec + dinputVec

# assume multivariate gauss distribution according to euler-step
# and calculate the likelihood
z = x_next - (x_now + f * ode_timestep_size[i])
v = (g %*% t(g) + small_identity) * ode_timestep_size[i]
nll = nll - RTMB::dmvnorm(z, Sigma=v, log=TRUE)
}
###### BETWEEN TIME POINTS LOOP END #######
}
###### TIME LOOP END #######

obsMat = RTMB::OBS(obsMat)
###### DATA UPDATE START #######
for(i in 1:NUMBER_OF_OBSERVATIONS){
for(j in 1:length(iobs[[i]])){

# Get index where observation is available
k = iobs[[i]][j]

# Get corresponding input, state and observation
inputVec = inputMat[k,]
stateVec = stateMat[ode_cumsum_timesteps[k]+1,]
obsScalar = obsMat[k,i]

# Observation equation and varianace
h_x = h__(stateVec, parVec, inputVec)[i]
hvar_x = hvar__(stateVec, parVec, inputVec)[i]

# likelihood contribution
nll = nll - RTMB::dnorm(obsScalar, mean=h_x, sd=sqrt(hvar_x), log=TRUE)
}}
###### DATA UPDATE END #######

# return
return(invisible(nll))
}')
  
  main.function.txt = stringr::str_replace_all(main.function.txt, 
                                               pattern="NUMBER_OF_STATES", 
                                               replacement=deparse(as.numeric(private$number.of.states)))
  main.function.txt = stringr::str_replace_all(main.function.txt, 
                                               pattern="NUMBER_OF_DIFFUSIONS", 
                                               replacement=deparse(as.numeric(private$number.of.diffusions)))
  main.function.txt = stringr::str_replace_all(main.function.txt, 
                                               pattern="NUMBER_OF_OBSERVATIONS", 
                                               replacement=deparse(as.numeric(private$number.of.observations)))
  main.function.txt = stringr::str_replace_all(main.function.txt, 
                                               pattern="STATE_NAMES", 
                                               replacement=paste("p$",private$state.names,sep="",collapse=","))
  main.function.txt = stringr::str_replace_all(main.function.txt, 
                                               pattern="FIXED_PARAMETERS", 
                                               replacement=paste("p$",private$parameter.names,sep="",collapse=","))
  eval(parse(text=main.function.txt))
  
  ################################################
  # Construct Neg. Log-Likelihood
  ################################################
  
  nll = RTMB::MakeADFun(func = laplace.nll, 
                        parameters=parameters, 
                        random=private$state.names,
                        map = lapply(private$fixed.pars, function(x) x$factor),
                        silent=TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
  
}
