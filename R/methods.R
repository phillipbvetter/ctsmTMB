#######################################################
# Print - S3 Method
#######################################################


#' Basic print of ctsmTMB objects
#' @param x a R6 ctsmTMB model object
#' @param ... additional arguments
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
#' @returns Print of ctsmTMB fit object
#' @export
print.ctsmTMB.fit = function(x,...) {
  fit <- x #method consistency (argument must be called x)
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


#######################################################
# Summary - S3 Method
#######################################################

#' Basic summary of objects of class 'ctsmTMB'
#' @param object a R6 ctsmTMB model object
#' @param correlation a boolean to indicate whether to return parameter correlations
#' @param ... additional arguments
#' @returns a summary of the model object
#' @export
summary.ctsmTMB = function(object,
                           correlation = FALSE,
                           ...) {
  obj = object$summary(correlation)
  
  return(invisible(obj))
}

#' Basic summary of ctsmTMB fit object
#' @param object a ctsmTMB fit object
#' @param correlation boolean indicating whether or not to display the
#' parameter correlation structure
#' @param ... additional arguments
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
# GGPLOT2 FUNCTIONS FOR USE IN PLOTS
#######################################################

getggplot2theme = function() {
  mytheme =
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text = ggplot2::element_text("Avenir Next Condensed",size=12),
      legend.text = ggplot2::element_text(size=12),
      axis.text = ggplot2::element_text(size=12),
      strip.text = ggplot2::element_text(face="bold",size=12),
      # panel.grid.major = element_blank(),
      # panel.grid.minor = element_blank(),
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.position = "top",
      plot.title = ggplot2::element_text(hjust=0.5)
    )
  return(mytheme)
}

# getggplot2colors = function(n) {
#   hues = seq(15, 375, length = n + 1)
#   ggcolors = grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
#   return(ggcolors)
# }

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
#' @param x A R6 ctsmTMB object
#' @param plot.obs a vector of integers to indicate which observations should be
#' plotted. When multiple are requested a list of plots, one for each observation
#' is returned instead.
#' @param ggtheme ggplot2 theme to use for creating the ggplot.
#' @param ... additional arguments
#' @returns a (list of) ggplot residual plot(s)
#' @export
plot.ctsmTMB = function(x,
                        plot.obs=1,
                        ggtheme=getggplot2theme(),
                        ...) {
  
  object <- x
  object$plot(plot.obs=plot.obs, ggtheme=ggtheme)
  
  return(invisible(NULL))
}

#' This function creates residual plots for an estimated ctsmTMB object
#' @param x A R6 ctsmTMB fit object
#' @param plot.obs a vector of integers to indicate which observations should be
#' plotted. When multiple are requested a list of plots, one for each observation
#' is returned instead.
#' @param ggtheme ggplot2 theme to use for creating the ggplot.
#' @param ... additional arguments
#' @returns a (list of) ggplot residual plot(s)
#' @export
plot.ctsmTMB.fit = function(x,
                            plot.obs=1,
                            ggtheme=getggplot2theme(),
                            ...) {
  
  fit <- x
  private <- fit$private
  
  if (!(inherits(ggtheme,"theme") & inherits(ggtheme,"gg"))) {
    stop("ggtheme must be a ggplot2 theme")
  }
  if(!is.numeric(plot.obs)){
    stop("plot.obs must be a numeric (or integer) value")
  }
  
  if(is.null(fit$residuals)){
    if(private$method=="laplace"){
      stop("No residuals to plot. Did you calculate residuals with argument 'laplace.residuals=TRUE' when calling 'estimate'?")
    }
    stop("Error: no residuals were found in the fit object.")
  }
  
  # retrieve user default parameter settings
  
  # use ggplot to plot
  # mycolor = getggplot2colors(2)[2]
  mycolor <- "steelblue"
  # if (use.ggplot) {
  plots = list()
  t = fit$residuals$normalized[["t"]]
  for (i in 1:private$number.of.observations) {
    e = fit$residuals$normalized[[private$obs.names[i]]]
    id = !is.na(e)
    e = e[id]
    t = t[id]
    nam = private$obs.names[i]
    
    # time vs residuals
    plot.res =
      ggplot2::ggplot(data=data.frame(t,e)) +
      ggplot2::geom_point(ggplot2::aes(x=t,y=e),color=mycolor) +
      ggtheme +
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
      ggplot2::ggplot(data=data.frame(e)) +
      ggplot2::geom_histogram(ggplot2::aes(x=e,y=ggplot2::after_stat(density)),bins=100,color="black",fill=mycolor) +
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
    plot.cpgram =
      ggfortify::ggcpgram(e,colour=mycolor) +
      ggtheme +
      ggplot2::labs(
        title = paste("Cumulative Periodogram: "),
        y = "",
        x = "Lag"
      )
    
    plots[[i]] = patchwork::wrap_plots(plot.res,
                                       plot.hist,
                                       plot.qq,
                                       plot.cpgram,
                                       plot.acf,
                                       plot.pacf ,
                                       nrow=3) +
      patchwork::plot_annotation(title=paste("Residual Analysis for ", nam),
                                 theme= ggtheme + ggplot2::theme(text=ggplot2::element_text("Avenir Next Condensed", size=18, face="bold"))
      )
    
  }
  
  # return plot list, and print one of the plots
  print(plots[[plot.obs]])
  return(invisible(plots))
}

#' #' Performs full multi-dimensional profile likelihood calculations
#' @param fit a ctmsTMB fit object
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
#' @param trace the optimization output flag (see \link[stats]{nlminb}) given to
#' the \code{control} argument.
#' @param control a list of optimization output controls (see \link[stats]{nlminb})
#' @note The implemetation was modified from that of
#' https://github.com/kaskr/adcomp/blob/master/TMB/R/tmbprofile.R
#' @export
profile.ctsmTMB.fit = function(fit,
                               parlist,
                               grid.size = rep(10, length(parlist)),
                               grid.qnt = rep(3, length(parlist)),
                               hessian=FALSE,
                               trace=0,
                               control=list(trace=trace,iter.max=1e3,eval.max=1e3)
){
  
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
      opt <- nlminb(x0, f, gr, he, control=control)
    } else {
      opt <- nlminb(x0, f, gr, control=control)
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
    cat("Iteration:", i, "/", n, "\n")
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
  
  # return
  returnlist = list(
    profile.grid.and.nll = X,
    parameter.values = parlist,
    parameter.pairs = prof.pars,
    optimiser.results = opt.list,
    full.likelihood.optimum = fit$par.fixed
  )
  class(returnlist) <- "ctsmTMB.profile"
  
  return(returnlist)
}

#' #' Plot a profile likelihood ctsmTMB object
#' @param x a profile.ctstTMB object
#' @param y not in use
#' @param include.opt boolean which indicates whether or not to include the
#' total likelihood optimizer in the plot.
#' @param ... additional arguments
#' @export
plot.ctsmTMB.profile = function(x,y,include.opt=TRUE,...){
  list <- x
  l <- length(list$parameter.values)
  df <- list$profile.grid.and.nll
  par.names <- head(names(df),l)
  opt <- list$full.likelihood.optimum[par.names]
  
  # A simple plot is only needed if we are profiling one parameter
  if(l==1L){
    p <- ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x=df[[par.names[1]]],y=df$nll)) +
      {if(include.opt) ggplot2::geom_vline(ggplot2::aes(xintercept=opt,color="Full Likelihood Optimum"))} +
      ggplot2::labs(color="",
           title="Profile Likelihood",
           y = "Negative Log-Likelihood") +
      getggplot2theme()
  } else if (l==2L){
    p <- ggplot2::ggplot() + 
      ggplot2::stat_contour(ggplot2::aes(
        x=df[[par.names[1]]],
        y=df[[par.names[2]]],
        z=df$nll), bins=100) + {
      if(include.opt) ggplot2::geom_point(ggplot2::aes(x=opt[1], y=opt[2],fill="Full Likelihood Optimum"))
        } + 
      ggplot2::labs(fill="",
           title="Profile Likelihood") +
      getggplot2theme()
  } else {
    stop("Profile likelihood plotting for more than two parameters is not supported.")
  }
  
  return(p)
}
