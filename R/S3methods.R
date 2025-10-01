#######################################################
# Print - S3 Method
#######################################################

#' Basic print of ctsmTMB objects
#' @param x  an object of class 'ctsmTMB'
#' @param ... additional arguments (not in use)
#' @examples
#' library(ctsmTMB)
#' model <- ctsmTMB$new()
#' 
#' # print empty model
#' print(model)
#' 
#' # add elements to model and see new print
#' model$addSystem(dx ~ theta * (mu+u-x) * dt + sigma_x*dw)
#' model$addObs(y ~ x)
#' model$setVariance(y ~ sigma_y^2)
#' model$addInput(u)
#' model$setParameter(
#'   theta   = c(initial = 1, lower=1e-5, upper=50),
#'   mu      = c(initial=1.5, lower=0, upper=5),
#'   sigma_x = c(initial=1, lower=1e-10, upper=30),
#'   sigma_y = 1e-2
#' )
#' print(model)
#' @returns Print of ctsmTMB model object
#' @export
print.ctsmTMB = function(x,...) {
  obj <- x
  output <- obj$print()
  
  return(invisible(output))
}

#' Basic print of objects ctsmTMB fit objects
#' @param x a ctsmTMB fit object
#' @param ... additional arguments
#' @examples
#' library(ctsmTMB)
#' model <- ctsmTMB$new()
#' 
#' # create model
#' model$addSystem(dx ~ theta * (mu+u-x) * dt + sigma_x*dw)
#' model$addObs(y ~ x)
#' model$setVariance(y ~ sigma_y^2)
#' model$addInput(u)
#' model$setParameter(
#'   theta   = c(initial = 1, lower=1e-5, upper=50),
#'   mu      = c(initial=1.5, lower=0, upper=5),
#'   sigma_x = c(initial=1, lower=1e-10, upper=30),
#'   sigma_y = 1e-2
#' )
#' model$setInitialState(list(1,1e-1))
#' 
#' # fit model to data
#' fit <- model$estimate(Ornstein)
#' 
#' # print fit
#' print(fit)
#' @returns Print of ctsmTMB fit object
#' @export
print.ctsmTMB.fit = function(x,...) {
  fit <- x
  if(is.null(fit$sd.fixed)) {
    fit$sd.fixed <- rep(NA,length(fit$par.fixed))
    fit$tvalue <- rep(NA,length(fit$par.fixed))
    fit$Pr.tvalue <- rep(NA,length(fit$par.fixed))
  }
  mat = cbind(fit$par.fixed, fit$sd.fixed, fit$tvalue, fit$Pr.tvalue)
  colnames(mat) = c("Estimate","Std. Error","t value","Pr(>|t|)")
  cat("Coefficent Matrix \n")
  stats::printCoefmat(mat)
  
  return(invisible(mat))
}

#' Basic summary of ctsmTMB fit object
#' @param object a ctsmTMB fit object
#' @param correlation boolean indicating whether or not to display the
#' parameter correlation structure
#' @param ... additional arguments
#' @examples
#' library(ctsmTMB)
#' model <- ctsmTMB$new()
#' 
#' # create model
#' model$addSystem(dx ~ theta * (mu+u-x) * dt + sigma_x*dw)
#' model$addObs(y ~ x)
#' model$setVariance(y ~ sigma_y^2)
#' model$addInput(u)
#' model$setParameter(
#'   theta   = c(initial = 1, lower=1e-5, upper=50),
#'   mu      = c(initial=1.5, lower=0, upper=5),
#'   sigma_x = c(initial=1, lower=1e-10, upper=30),
#'   sigma_y = 1e-2
#' )
#' model$setInitialState(list(1,1e-1))
#' 
#' # fit model to data
#' fit <- model$estimate(Ornstein)
#' 
#' # print model summary
#' summary(fit, correlation=TRUE)
#' @returns a summary of the estimated ctsmTMB model fit
#' @export
summary.ctsmTMB.fit = function(object, 
                               correlation = FALSE,
                               ...) {
  
  fit <- object #method consistency (argument must be called 'object')
  
  if (!is.logical(correlation)) {
    stop("correlation must be logical")
  }
  
  mat = cbind(fit$par.fixed,fit$sd.fixed,fit$tvalue,fit$Pr.tvalue)
  colnames(mat) = c("Estimate","Std. Error","t value","Pr(>|t|)")
  cat("Coefficent Matrix \n")
  stats::printCoefmat(mat)
  if (correlation){
    cat("\nCorrelation Matrix\n")
    S = diag(1/fit$sd.fixed)
    cor = S %*% fit$cov.fixed %*% S
    rownames(cor) = names(fit$par.fixed)
    colnames(cor) = names(fit$par.fixed)
    cor <- format(round(cor, 2), nsmall = 2, digits = 6)
    cor[!lower.tri(cor,diag=T)] <- ""
    print(cor, drop = FALSE, quote = FALSE)
    # cor[!lower.tri(cor)] <- ""
    # print(cor[-1, -(dim(cor)[1]), drop = FALSE], quote = FALSE)
  }
  
  return(invisible(list(parameters=mat)))
}

#######################################################
# Plot - S3 Method
#######################################################

#' Plot of k-step predictions from a ctsmTMB prediction object
#' @param x a ctsmTMB.pred object
#' @param y not used
#' @param k.ahead an integer indicating which k-ahead predictions to plot
#' @param state.name a string indicating which states to plot
#' @param type one of 'states' or 'observations', to plot
#' @param against name of an observations to plot predictions against
#' @param ... additional arguments
#' @examples
#' library(ctsmTMB)
#' model <- ctsmTMB$new()
#' 
#' # create model
#' model$addSystem(dx ~ theta * (mu+u-x) * dt + sigma_x*dw)
#' model$addObs(y ~ x)
#' model$setVariance(y ~ sigma_y^2)
#' model$addInput(u)
#' model$setParameter(
#'   theta   = c(initial = 1, lower=1e-5, upper=50),
#'   mu      = c(initial=1.5, lower=0, upper=5),
#'   sigma_x = c(initial=1, lower=1e-10, upper=30),
#'   sigma_y = 1e-2
#' )
#' model$setInitialState(list(1,1e-1))
#' 
#' # fit model to data
#' fit <- model$estimate(Ornstein)
#' 
#' # perform moment predictions
#' pred <- model$predict(Ornstein)
#' 
#' # plot the k.ahead=10 predictions
#' plot(pred, against="y.data")
#' 
#' 
#' # plot filtered states
#' plot(fit, type="states", against="y")
#' @returns A plot of predicted states
#' @export
plot.ctsmTMB.pred = function(x, 
                             y,
                             k.ahead = unique(x[["states"]][["k.ahead"]]),
                             state.name = NULL,
                             type="states",
                             against=NULL,
                             ...) {
  
  # method consistency
  states <- x[["states"]]
  obs <- x[["observations"]]
  
  # check n.ahead
  if(!any(k.ahead %in% unique(states$k.ahead))){
    stop("n.ahead not found in the prediction data frame")
  }
  
  # set state name to plot
  if(is.null(state.name)){
    state.name = colnames(states)[6]
  }
  
  # filter and plot
  bool = states$k.ahead %in% k.ahead
  x = states[bool, "t.j"]
  y = states[bool, state.name]
  plot(x=x, y=y, type="l",...)
  
  if(!is.null(against)){
    z = obs[bool, against]
    graphics::points(x=x, y=z,col="red",pch=16)
  }
  
  
  return(invisible(NULL))
}

#' This function creates residual plots for an estimated ctsmTMB object
#' @param x A R6 ctsmTMB fit object
#' @param print.plot a single integer determining which element out of all
#' states/observations (depending on the argument to \code{type}).
#' @param type a character vector either 'residuals' or 'states' determining what
#' to plot.
#' @param state.type a character vector either 'prior', 'posterior' or 'smoothed'
#' determining what kind of states to plot.
#' @param against.obs name of an observation to plot state predictions against.
#' @param ggtheme ggplot2 theme to use for creating the ggplot.
#' @param ylims limits on the y-axis for residual time-series plot
#' @param residual.burnin integer N to remove the first N residuals
#' @param residual.vs.obs.and.inputs the residual plots also include a new window
#' with time-series plots of residuals, associated observations and inputs
#' @param ... additional arguments
#' @examples
#' library(ctsmTMB)
#' model <- ctsmTMB$new()
#' 
#' # create model
#' model$addSystem(dx ~ theta * (mu+u-x) * dt + sigma_x*dw)
#' model$addObs(y ~ x)
#' model$setVariance(y ~ sigma_y^2)
#' model$addInput(u)
#' model$setParameter(
#'   theta   = c(initial = 1, lower=1e-5, upper=50),
#'   mu      = c(initial=1.5, lower=0, upper=5),
#'   sigma_x = c(initial=1, lower=1e-10, upper=30),
#'   sigma_y = 1e-2
#' )
#' model$setInitialState(list(1,1e-1))
#' 
#' # fit model to data
#' fit <- model$estimate(Ornstein)
#' 
#' # plot residuals
#' \dontrun{plot(fit)}
#' 
#' # plot filtered states
#' \dontrun{plot(fit, type="states")}
#' @returns a (list of) ggplot residual plot(s)
#' @export
plot.ctsmTMB.fit = function(x,
                            print.plot=1,
                            type="residuals",
                            state.type="prior",
                            against.obs=NULL,
                            ggtheme=getggplot2theme(),
                            ylims = c(NA,NA),
                            residual.burnin=0L,
                            residual.vs.obs.and.inputs = FALSE,
                            ...) {
  
  fit <- x
  private <- fit$private
  
  if (!(inherits(ggtheme,"theme") & inherits(ggtheme,"gg"))) {
    stop("The provided theme is not a ggtheme")
  }
  if(!is.numeric(print.plot)){
    stop("print.plot must be a numeric value")
  }
  if(!any(state.type %in% c("prior","posterior","smoothed"))){
    stop("The state.type must be one of 'prior', 'posterior' or 'smoothed'.")
  }
  
  plots <- plots2 <- list()
  
  # residual plots
  ########################################
  if(type=="residuals"){
    
    if(is.null(fit$residuals)){
      if(private$method=="laplace"){
        stop("No residuals to plot. Did you calculate residuals with argument 'laplace.residuals=TRUE' when calling 'estimate'?")
      }
      stop("Error: no residuals were found in the fit object.")
    }
    
    mycolor <- "steelblue"
    t <- fit$residuals$residuals[["t"]]
    if(residual.burnin) t <- tail(t,-residual.burnin)
    
    for (i in 1:private$number.of.observations) {
      
      e = fit$residuals$normalized[[private$obs.names[i]]]
      e0 <- fit$residuals$residuals[[private$obs.names[i]]]
      if(residual.burnin){
        e <- tail(e,-residual.burnin)
        e0 <- tail(e0,-residual.burnin)
      }
      id = !is.na(e)
      e = e[id]
      t = t[id]
      nam = private$obs.names[i]
      
      # time vs residuals
      plot.res =
        # ggplot2::ggplot(data=data.frame(t,e)) +
        # ggplot2::geom_point(ggplot2::aes(x=t,y=e),color=mycolor) +
        ggplot2::ggplot(data=data.frame(t,e0)) +
        ggplot2::geom_line(ggplot2::aes(x=t,y=e0),color=mycolor) +
        ggtheme +
        coord_cartesian(ylim=ylims) + 
        ggplot2::labs(
          title = paste("Time Series of Residuals: "),
          y = "",
          x = "Time"
        )
      # quantile plots
      plot.qq =
        ggplot2::ggplot(data=data.frame(e)) +
        ggplot2::stat_qq(ggplot2::aes(sample=e),color=mycolor) +
        ggplot2::stat_qq_line(ggplot2::aes(sample=e),lty="dashed") +
        ggtheme +
        ggplot2::labs(
          title = paste("Normal Q-Q Plot: "),
          y = "Sample Quantiles",
          x = "Theoretical Quantiles"
        )
      # histogram
      plot.hist =
        # ggplot2::ggplot(data=data.frame(e)) +
        # ggplot2::geom_histogram(ggplot2::aes(x=e,y=ggplot2::after_stat(density)),bins=100,color="black",fill=mycolor) +
        ggplot2::ggplot(data=data.frame(e0)) +
        ggplot2::geom_histogram(ggplot2::aes(x=e0,y=ggplot2::after_stat(density)),bins=100,color="black",fill=mycolor) +
        ggtheme +
        ggplot2::labs(
          title = paste("Histogram: "),
          y = "",
          x = ""
        )
      # acf
      myacf = stats::acf(e,na.action=na.pass,plot=FALSE)
      plot.acf =
        ggplot2::ggplot(data=data.frame(lag=myacf$lag[-1],acf=myacf$acf[-1])) +
        ggplot2::geom_errorbar(ggplot2::aes(x=lag,ymax=acf,ymin=0),width=0,color=mycolor) +
        ggplot2::geom_hline(yintercept=0) +
        ggplot2::geom_hline(yintercept=c(-2/sqrt(myacf$n.used),2/sqrt(myacf$n.used)),color="blue",lty="dashed") +
        ggtheme +
        ggplot2::coord_cartesian(xlim=c(0,ceiling(max(myacf$lag)/10)*10)) +
        ggplot2::labs(
          title = paste("Auto-Correlation: "),
          y = "",
          x = "Lag"
        )
      # pacf
      mypacf = stats::pacf(e,na.action=na.pass,plot=FALSE)
      plot.pacf =
        ggplot2::ggplot(data=data.frame(lag=mypacf$lag[-1], pacf=mypacf$acf[-1])) +
        ggplot2::geom_errorbar(ggplot2::aes(x=lag,ymax=pacf,ymin=0),width=0,color=mycolor) +
        ggplot2::geom_hline(yintercept=0) +
        ggplot2::geom_hline(yintercept=c(-2/sqrt(mypacf$n.used),2/sqrt(mypacf$n.used)),color="blue",lty="dashed") +
        ggtheme +
        ggplot2::coord_cartesian(xlim=c(0,ceiling(max(mypacf$lag)/10)*10)) +
        ggplot2::labs(
          title = paste("Partial Auto-Correlation: "),
          y = "",
          x = "Lag"
        )
      # cpgram
      # if (requireNamespace("ggfortify", quietly = TRUE)) {
      plot.cpgram =
        ggfortify_ggcpgram(e,colour=mycolor) +
        ggtheme +
        ggplot2::labs(
          title = paste("Cumulative Periodogram: "),
          y = "",
          x = "Lag"
        )
      # } else {
      #   plot.cpgram = ggplot2::ggplot() + 
      #     ggplot2::geom_text(
      #       aes(x=0.5,
      #           y=0.5,
      #           label="Install 'ggfortify' \n for \n Cumulative Periodogram")
      #     ) +
      #     theme_void()
      # }
      
      plots[[i]] = patchwork::wrap_plots(plot.res,
                                         plot.hist,
                                         plot.qq,
                                         plot.cpgram,
                                         plot.acf,
                                         plot.pacf ,
                                         nrow=3) +
        patchwork::plot_annotation(title=paste("Residual Analysis for ", nam),
                                   subtitle=paste("modelname:", private$modelname),
                                   theme = ggtheme + 
                                     ggplot2::theme(text=ggplot2::element_text(size=12, face="bold"))
        )
      
      if(residual.vs.obs.and.inputs){
        
        y.obs <- private$data[[nam]]
        y <- private$data[[private$obs.names[i]]]
        if(residual.burnin){
          y.obs <- tail(y.obs, -residual.burnin)
          y <- tail(y, -residual.burnin)
        }
        
        # time vs observations
        plot.obs =
          ggplot2::ggplot(data=data.frame(t,y.obs)) +
          # ggplot2::geom_point(ggplot2::aes(x=t,y=y,color="Predictions"),linewidth=0.5) +
          ggplot2::geom_line(ggplot2::aes(x=t,y=y.obs,color="Observations"),linewidth=1.5) +
          ggplot2::geom_line(ggplot2::aes(x=t,y=y,color="Predictions"),linewidth=0.5) +
          ggtheme +
          theme(
            # legend.position=c(0.9, 0.9),
            # legend.direction = "vertical",
            # legend.box.background = element_rect(fill="grey"),
            # legend.box ="vertical",
            legend.margin = ggplot2::margin(0,0,0,0),
            legend.box.margin = ggplot2::margin(0,0,0,0),
            legend.spacing.y= unit(0,"mm")
          ) +
          # theme(legend.position="inside") +
          scale_color_manual(values=c(mycolor,"tomato")) +
          ggplot2::labs(
            title = "Observation",
            y = nam,
            x = "",
            color=""
          )
        
        # time vs input plots
        # input plots are constant across obs plot index "i"
        if(i==1){
          input.plots <- list()
          input.names <- private$input.names
          k <- private$number.of.inputs-1
          for(j in 1:k){
            y.input <- private$data[[input.names[j+1]]]
            if(residual.burnin){
              y.input <- tail(y.input, -residual.burnin)
            }
            input.plots[[j]] <- 
              ggplot2::ggplot(data=data.frame(t,y.input)) +
              ggplot2::geom_line(ggplot2::aes(x=t,y=y.input),color=mycolor) +
              ggtheme +
              ggplot2::labs(
                y = input.names[j+1],
                x = ""
              )
          }
          # only the bottom plot x-axis should be titled "Time"
          input.plots[[k]]$labels$x = "Time"
          input.plots[[1]]$labels$title = "Inputs"
          plot.inputs <- patchwork::wrap_plots(input.plots,ncol=1)
        }
        
        plot.res2 <- plot.res
        plot.res2$labels$x <- ""
        plot.res2$labels$title <- "Residuals"
        
        # save current residual vs input plot
        plots2[[i]] <- patchwork::wrap_plots(plot.res2,
                                             plot.obs,
                                             plot.inputs,
                                             heights=c(1,1,k)
        ) +
          patchwork::plot_annotation(title=paste("Residual vs Inputs for", nam),
                                     subtitle=paste("modelname:", private$modelname),
                                     theme = ggtheme + 
                                       ggplot2::theme(text=ggplot2::element_text(size=12, face="bold")))
      }
    }
    
    
  }
  
  # state plots
  ########################################
  if(type=="states"){
    
    mycolors <- c("steelblue","tomato")
    t <- fit$private$data$t
    
    for (i in 1:private$number.of.states) {
      nam <- private$state.names[i]
      y.mean <- fit$states$mean[[state.type]][[nam]]
      y.sd <- fit$states$sd[[state.type]][[nam]]
      y.lab <- sprintf("%s (%s)", capitalize_first(state.type), nam)
      
      tempdata <- data.frame(t=t, y=y.mean, sd=y.sd)
      plots[[i]] <- ggplot2::ggplot() +
        ggplot2::geom_ribbon(data=tempdata, ggplot2::aes(x=t, ymin=y-2*sd, ymax=y+2*sd), fill="grey", alpha=0.6) +
        ggplot2::geom_line(data=tempdata, ggplot2::aes(x=t, y=y, color=nam)) +
        { 
          if(i==print.plot && !is.null(against.obs))
          {
            y.obs <- fit$private$data[[against.obs]]
            tempobsdata <- data.frame(t=t, y=y.obs)
            ggplot2::geom_point(data=tempobsdata, ggplot2::aes(x=t, y=y, color=against.obs), size=0.5)
          }
        } +
        ggplot2::labs(
          title = y.lab,
          color="",
          x = "Time",
          y = "",
        ) +
        {
          if(i==print.plot && !is.null(against.obs))
          {
            ggplot2::scale_color_manual(values=mycolors, labels=c(against.obs, nam))
          } else {
            ggplot2::scale_color_manual(values=mycolors, labels=nam)
          }
        } +
        getggplot2theme()
      
    }
  }
  
  # print the first plot to the console
  grDevices::dev.new()
  print(plots[[print.plot]])
  
  # print second plot if requested
  if(length(plots2) > 0){
    grDevices::dev.new()
    print(plots2[[print.plot]])
  }
  
  # return plot list
  return(invisible(c(plots,plots2)))
}

#' Performs full multi-dimensional profile likelihood calculations
#' @param fitted a ctmsTMB fit object
#' @param ... various arguments (not in use)
#' @param parlist a named-list of parameters to profile over. The user can either
#' supply grid-values in the list or leave it empty. If the any one list is empty
#' then grid-values will be calculated using the estimated parameter mean value
#' and standard deviation.
#' @param grid.size a vector of \code{length(parlist)} indicating the number
#' of grid-points along each parameter direction. This is only used if the
#' \code{parlist} is empty.
#' @param grid.qnt a vector of \code{length(parlist)} determining the width of
#' the grid points from the mean value in multiples of the standard deviation.
#' @param hessian a boolean indicating whether to use the hessian or not during
#' the profile optimization.
#' @param silent boolean whether or not to mute current iteration number
#' the \code{control} argument.
#' @param control a list of optimization output controls (see \link[stats]{nlminb})
#' @examples
#' library(ctsmTMB)
#' model <- ctsmTMB$new()
#' 
#' # create model
#' model$addSystem(dx ~ theta * (mu+u-x) * dt + sigma_x*dw)
#' model$addObs(y ~ x)
#' model$setVariance(y ~ sigma_y^2)
#' model$addInput(u)
#' model$setParameter(
#'   theta   = c(initial = 1, lower=1e-5, upper=50),
#'   mu      = c(initial=1.5, lower=0, upper=5),
#'   sigma_x = c(initial=1, lower=1e-10, upper=30),
#'   sigma_y = 1e-2
#' )
#' model$setInitialState(list(1,1e-1))
#' 
#' # fit model to data
#' fit <- model$estimate(Ornstein)
#' 
#' # calculate profile likelihood
#' out <- profile(fit,parlist=list(theta=NULL))
#' @note The implementation was modified from that of
#' https://github.com/kaskr/adcomp/blob/master/TMB/R/tmbprofile.R
#' @export
profile.ctsmTMB.fit = function(fitted,
                               parlist,
                               grid.size = rep(10, length(parlist)),
                               grid.qnt = rep(3, length(parlist)),
                               hessian=FALSE,
                               silent=FALSE,
                               control=list(trace=0,iter.max=1e3,eval.max=1e3),
                               ...
){
  
  fit <- fitted
  
  if(missing(fit)){
    stop("Please supply a fit from a ctsmTMB model.")
  }
  
  bool <- !(names(parlist) %in% names(fit$par.fixed))
  if(any(bool)){
    stop("The following parameter name(s) do not exist in the model:\n", paste(names(parlist)[bool],collapse=", "))
  }
  
  # 1. Get likelihood function
  nll <- fit$private$nll
  
  # Create parameter grid
  parnames <- names(parlist)
  id <- names(fit$par.fixed) %in% parnames
  if(all(lapply(parlist,length) == 0)){
    for(i in seq_along(parnames)){
      p = fit$par.fixed[parnames[i]]
      s = fit$sd.fixed[parnames[i]]
      if(is.na(s)) stop("No std. error available for ",parnames[i], ". Please supply your own grid.")
      parlist[[parnames[i]]] = seq(p-grid.qnt[i]*s, p+grid.qnt[i]*s, length.out=grid.size[i])
    }
  }
  len <- length(parlist)
  X <- do.call(expand.grid, rev(parlist))
  names(X) <- rev(parnames)
  X <- X[,len:1,drop=FALSE]
  X <- cbind(X, nll=0)
  n <- nrow(X)
  prof.nll <- numeric(n)
  prof.pars <- vector("list",length=n)
  opt.list <- vector("list",length=n)
  
  # 2. Create parameter filter matrix C that extracts only the
  # non-profiled entries which needs to be optimized over.
  # C maps from length(par.fixed) to length(par.fixed) - length(parnames)
  # by multiplication.
  C <- diag(length(fit$par.fixed))[,!id,drop=F]
  
  # Create optimization function that takes initial guess x0
  # and finds the maximum profile likelihood estimate in the
  # reduced parameter space
  f.optim <- function(x0){
    
    # objective function
    f <- function(x){
      y <- par + as.numeric(C %*% x)
      nll$fn(y)
    }
    
    # gradient
    gr <- function(x){
      y <- par + as.numeric(C %*% x)
      as.numeric(nll$gr(y) %*% C)
    }
    
    # hessian
    he <- function(x){
      y <- par + as.numeric(C %*% x)
      # t(C) %*% nll$he(y)[!id,!id] %*% C
      nll$he(y)[!id,!id]
    }
    
    # optimize
    if(hessian) {
      opt <- stats::nlminb(x0, f, gr, he, control=control)
    } else {
      opt <- stats::nlminb(x0, f, gr, control=control)
    }
    return(opt)
  }
  
  # Robustify f.optim to handle NA cases
  f <- function(x0){
    y <- try_withWarningRecovery(f.optim(x0), silent=TRUE)
    if(inherits(y,"try-error")) y <- NA
    return(y)
  }
  
  names.in.parlist <- names(parlist)
  
  par <- fit$par.fixed
  x0 <- fit$par.fixed[!id]
  par[!id] <- 0
  # Run for-loop over all pairs
  for(i in 1:nrow(X)){
    # par[id] <- unlist(X[i, -(len+1)])
    par[names.in.parlist] = unlist(X[i,names.in.parlist])
    par.temp <- par
    #
    if(!silent) cat("Iteration:", i, "/", n, "\n")
    opt <- f(x0)
    #
    opt.list[[i]] <- opt
    # store results and set new initial parameter guess
    if(any(is.na(opt))){
      prof.nll[i] <- NA
      prof.pars[[i]] <- NA
      x0 <- fit$par.fixed[!id]
    } else {
      # prof.nll[i] <- opt$objective
      X[i,len+1] <- opt$objective
      par.temp[!id] <- opt$par
      prof.pars[[i]] <- par.temp
      
      # This next start guess is optimum found 
      # in the previous iteration
      # This may be bad when parameter grid jumps??
      x0 <- opt$par
      if(opt$convergence==1) {
        x0 <- fit$par.fixed[!id]
      }
    }
  }
  
  # return -----------------
  # normalize and convert nll into likelihood
  X[,len+1] <- exp(fit$nll - X[,len+1])
  names(X)[ncol(X)] <- "likelihood"
  
  # create return list
  returnlist = list(
    profile.grid.and.likelihood = X,
    parameter.values = parlist,
    parameter.pairs = prof.pars,
    optimiser.results = opt.list,
    full.likelihood.optimum = fit$par.fixed
  )
  # set class so that S3 methods work
  class(returnlist) <- "ctsmTMB.profile"
  
  return(returnlist)
}

#' Plot a profile likelihood ctsmTMB object
#' @param x a profile.ctsmTMB object
#' @param y not in use
#' @param include.opt boolean which indicates whether or not to include the
#' total likelihood optimizer in the plot.
#' @param ... additional arguments
#' @examples
#' library(ctsmTMB)
#' model <- ctsmTMB$new()
#' 
#' # create model
#' model$addSystem(dx ~ theta * (mu+u-x) * dt + sigma_x*dw)
#' model$addObs(y ~ x)
#' model$setVariance(y ~ sigma_y^2)
#' model$addInput(u)
#' model$setParameter(
#'   theta   = c(initial = 1, lower=1e-5, upper=50),
#'   mu      = c(initial=1.5, lower=0, upper=5),
#'   sigma_x = c(initial=1, lower=1e-10, upper=30),
#'   sigma_y = 1e-2
#' )
#' model$setInitialState(list(1,1e-1))
#' 
#' # fit model to data
#' fit <- model$estimate(Ornstein)
#' 
#' # calculate profile likelihood
#' # out <- profile(fit,parlist=list(theta=NULL))
#' out <- profile(fit,parlist=list(theta=NULL, mu=NULL))
#' 
#' # plot profile
#' plot(out)
#' @export
plot.ctsmTMB.profile = function(x, y, include.opt=TRUE,...){
  
  list <- x
  l <- length(list$parameter.values)
  df <- list$profile.grid.and.likelihood
  par.names <- head(names(df), l)
  opt <- list$full.likelihood.optimum[par.names]
  
  if(l==1L){
    
    threshold <- exp(-qchisq(0.95, df=1)/2)
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_hline(ggplot2::aes(yintercept=threshold, color="Likelihood CI"), lty="dashed") +
      ggplot2::geom_line(ggplot2::aes(x=df[[par.names[1]]],y=df$likelihood)) +
      {if(include.opt) ggplot2::geom_vline(ggplot2::aes(xintercept=opt,color="Full Likelihood Optimum"),lty="dashed")} +
      ggplot2::labs(color="",
                    title="Profile Likelihood",
                    y = "",
                    x = par.names[1]) +
      getggplot2theme()
    
  } else if (l==2L){
    
    quants <- c(0.75, 0.90, 0.95, 0.99)
    chisquared_contours <- exp(-qchisq(quants, df=2)/2)
    labelfun <- function(x){
      y <- (1-x) * 100
      y <- paste0(y,"%")
    }
    
    if (requireNamespace("geomtextpath", quietly = TRUE)) {
      # Use geomtextpath if available
      p <- ggplot2::ggplot() +
        geomtextpath::geom_textcontour(
          data = df, 
          ggplot2::aes(
            x = .data[[par.names[1]]], 
            y = .data[[par.names[2]]], 
            z = .data$likelihood, 
            label = labelfun(after_stat(level)),
            color = after_stat(level)
          ),
          breaks = chisquared_contours,
          linewidth = 0.6
        )
    } else {
      # Fallback without labels if geomtextpath is not available
      message("Please install 'geomtextpath' to get contour labels.")
      p <- ggplot2::ggplot() +
        ggplot2::geom_contour(
          data = df,
          ggplot2::aes(
            x = .data[[par.names[1]]],
            y = .data[[par.names[2]]],
            z = .data$likelihood,
            color = after_stat(level)
          ),
          breaks = chisquared_contours,
          linewidth = 0.6
        )
    }
    
    # add labels
    p <- p +
      ggplot2::labs(
        color="",
        title="Profile Likelihood",
        x = par.names[1],
        y = par.names[2]
      ) +
      scale_color_continuous(guide="none") +
      getggplot2theme()
    
    if (include.opt) {
      p <- p + ggplot2::geom_point(
        data = data.frame(x = opt[1], y = opt[2]), 
        ggplot2::aes(x = x, y = y), 
        size = 1
      )
    }
    
  } else {
    
    stop("Profile likelihood plotting for more than two parameters is not supported.")
    
  }
  
  return(p)
}

