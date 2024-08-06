#######################################################
# Print - S3 Method
#######################################################


#' Basic print of objects of class 'ctsmTMB'
#' @returns A huge amount of information
#' @export
print.ctsmTMB = function(object,...) {
  
  obj = object$print()
  #
  return(invisible(obj))
  
}

#' Basic print of objects of class 'ctsmTMB'
#' @returns A huge amount of information
#' @export
print.ctsmTMB.fit = function(fit) {
  
  mat = cbind(fit$par.fixed,fit$sd.fixed,fit$tvalue,fit$Pr.tvalue)
  colnames(mat) = c("Estimate","Std. Error","t value","Pr(>|t|)")
  cat("Coefficent Matrix \n")
  stats::printCoefmat(mat)
  
  return(invisible(mat))
  
}


#######################################################
# Summary - S3 Method
#######################################################

#' Basic summary of objects of class 'ctsmTMB'
#' @returns A huge amount of information
#' @export
summary.ctsmTMB = function(object,correlation=FALSE) {
  
  obj = object$summary(correlation)
  #
  return(invisible(obj))
  
}

#' @returns summary of fit object from \code{obj$estimate()}
#' @export
summary.ctsmTMB.fit = function(fit, correlation=FALSE) {
  
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
      legend.box = "vertical",
      legend.position = "top",
      plot.title = ggplot2::element_text(hjust=0.5)
    )
  return(mytheme)
}

getggplot2colors = function(n) {
  hues = seq(15, 375, length = n + 1)
  ggcolors = hcl(h = hues, l = 65, c = 100)[1:n]
  return(ggcolors)
}

#######################################################
# Plot - S3 Method
#######################################################

#' Basic summary of objects of class 'ctsmTMB.pred' from predict
#' @returns A huge amount of information
#' @export
plot.ctsmTMB.pred = function(pred.data,
                             n.ahead=0,
                             state.name=NULL,
                             ...) {
  
  # check n.ahead
  if(!any(n.ahead %in% unique(pred$n.ahead))){
    stop("n.ahead not found in the prediction data frame")
  }
  
  # set state name to plot
  if(is.null(state.name)){
    state.name = colnames(pred)[6]
  }
  
  # filter and plot
  bool = pred$n.ahead==n.ahead
  x = pred[bool,"t_{k+i}"]
  y = pred[bool,state.name]
  plot(x=x, y=y, type="l",...)
  
  
  return(invisible(NULL))
}


#' Basic summary of objects of class 'ctsmTMB'
#' @param plot.obs a vector to indicate which observations should be plotted for. If multiple
#' are chosen a list of plots for each observation is returned.
#' @param pacf logical to indicate whether or not the partial autocorrelations should be returned.
#' The default is FALSE in which case a histogram is returned instead.
#' @param extended logical. if TRUE additional information is printed
#' @param ggtheme ggplot2 theme to use for creating the ggplot.
#' @returns A huge amount of information
#' @export
plot.ctsmTMB = function(object,
                        plot.obs=1,
                        ggtheme=getggplot2theme()) {
  
  object$plot(plot.obs=plot.obs, ggtheme=ggtheme)
  
  return(invisible(NULL))
}

#' Basic summary of objects of class 'ctsmTMB'
#' @param plot.obs a vector to indicate which observations should be plotted for. If multiple
#' are chosen a list of plots for each observation is returned.
#' @param pacf logical to indicate whether or not the partial autocorrelations should be returned.
#' The default is FALSE in which case a histogram is returned instead.
#' @param extended logical. if TRUE additional information is printed
#' @param ggtheme ggplot2 theme to use for creating the ggplot.
#' @returns A list of plots
#' @export
#' @importFrom stats frequency fft spec.taper
plot.ctsmTMB.fit = function(fit,
                            plot.obs=1,
                            ggtheme=getggplot2theme()) {
  
  # get private fields from object
  private <- fit$private
  
  if (!(inherits(ggtheme,"theme") & inherits(ggtheme,"gg"))) {
    stop("ggtheme must be a ggplot2 theme")
  }
  if(!is.numeric(plot.obs)){
    stop("plot.obs must be a numeric (or integer) value")
  }
  if(plot.obs > private$number.of.states){
    message("Can't plot state ", plot.obs, " because there's only ", private$number.of.states, " state(s). Setting plot.obs = ",private$number.of.states)
    plot.obs = private$number.of.states
  }
  
  if(private$method == "laplace"){
    message("'plot' is not available for method 'laplace' yet.")
    return(NULL)
  }
  
  # retrieve user default parameter settings
  
  # use ggplot to plot
  mycolor = getggplot2colors(2)[2]
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
    
    # return plot list, and print one of the plots
    print(plots[[plot.obs]])
    return(invisible(plots))
  }
  
  # return
  return(invisible(NULL))
}

#' #' Full multi-dimensional profile likelihood calculations
#' @param parnames a named list of parameter values, one for each a vector of parameter names to be profiled.
#' @param initial values for all parameters
#' This is inspired by the TMB implementation at
#' https://github.com/kaskr/adcomp/blob/master/TMB/R/tmbprofile.R
#' @export
profile.ctsmTMB.fit = function(fit, 
                               parlist, 
                               par.vals=fit$par.fixed, 
                               silent=F, 
                               trace=0, 
                               control=list(trace=trace,iter.max=1e5,eval.max=1e5)
                               ){
  
  if(missing(fit)){
    stop("Please supply a fit from a ctsmTMB model.")
  }
  if(missing(parlist)){
    stop("Please supply a named list of parameter values")
  }
  
  # 1. names of parameters to be profiled
  # 2. initial values for optimizer
  
  # 1. Get likelihood function
  nll <- fit$private$nll
  
  X <- do.call(expand.grid, rev(parlist))
  names(X) <- rev(names(parlist))
  n <- nrow(X)
  prof.nll <- numeric(n)
  prof.pars <- vector("list",length=n)
  opt.list <- vector("list",length=n)
  
  # 2. Create parameter filter matrix that extracts only the non-profiled entries
  # which needs to be optimized over
  # Must map from length(par.fixed) to length(par.fixed) - length(parnames)
  parnames <- names(parlist)
  id <- names(fit$par.fixed) %in% parnames
  C <- diag(length(fit$par.fixed))[,!id,drop=F]
  
  f.optim <- function(x0){
    
    # print(par)
    
    # objective function
    f <- function(x){
      par0 <- par + as.numeric(C %*% x)
      nll$fn(par0)
    }
    
    # gradient
    gr <- function(x){
      par0 <- par + as.numeric(C %*% x)
      as.numeric(nll$gr(par0) %*% C)
    }
    
    # optimize
    opt <- nlminb(x0, f, gr, control=control)
    
    return(opt)
  }
  
  # Robustify f
  f <- function(x0){
    y <- try_withWarningRecovery(f.optim(x0), silent=TRUE)
    if(is(y, "try-error")) y <- NA
    return(y)
  }
  
  # Parameter values
  par <- par.vals
  x0 <- par[!id]
  par[!id] <- 0
  partemp <- par.vals
  for(i in 1:nrow(X)){
    par[id] <- unlist(X[i,])
    par.temp <- par
    # 
    if(!silent){
      cat("Iteration:", i, "/", n, "\n")
    }
    opt <- f(x0)
    # 
    opt.list[[i]] <- opt
    if(any(is.na(opt))){
      prof.nll[i] <- NA
      prof.pars[[i]] <- NA
      x0 <- fit$par.fixed[!id]
    } else {
      prof.nll[i] <- opt$objective
      par.temp[!id] <- opt$par
      prof.pars[[i]] <- par.temp
      if(opt$convergence==1) {
        x0 <- fit$par.fixed
      } else {
        x0 <- opt$par 
      }
    }
  }
  
  # test
  # x <- par.vals[!(names(par.vals) %in% parnames)]
  # new.f(x)
  # new.gr(x)
  
  # run optimizer
  # pars0 <- par.vals[!id]
  # opt <- nlminb(pars0, f, gr)
  returnlist = list(
    profile.grid = cbind(X, profile.nll = prof.nll),
    profile.pars = prof.pars,
    profile.opts = opt.list
  )
  
  # return
  # return(list(profileNegLogLikelihood = prof.nll, profileParameters = prof.pars, optLists=opt.list))
  return(returnlist)
}

#' #' Inner profile likelihood calculation
#' @param parnames a vector of parameter names to be profiled.
#' @param initial values for all parameters
#' This is inspired by the TMB implementation at
#' https://github.com/kaskr/adcomp/blob/master/TMB/R/tmbprofile.R
#' @export
profile0.ctsmTMB.fit = function(fit, parnames, par.vals=fit$par.fixed){
  
  if(missing(fit)){
    stop("Please supply a fit from a ctsmTMB model.")
  }
  if(missing(parnames)){
    stop("Please supply a parameter names")
  }
  
  # 1. names of parameters to be profiled
  # 2. initial values for optimizer
  
  # 1. Get likelihood function
  nll <- fit$private$nll
  
  # Parameter values
  par <- par.vals
  x0 <- par[!id]
  par[!id] <- 0
  
  # 2. Create parameter filter matrix that extracts only the non-profiled entries
  # which needs to be optimized over
  # Must map from length(par.fixed) to length(par.fixed) - length(parnames)
  id <- names(fit$par.fixed) %in% parnames
  C <- diag(length(fit$par.fixed))[,!id,drop=F]
  
  
  f.optim <- function(x0){
    
    # objective function
    f <- function(x){
      par0 <- par + as.numeric(C %*% x)
      nll$fn(par0)
    }
    
    # gradient
    gr <- function(x){
      par0 <- par + as.numeric(C %*% x)
      as.numeric(nll$gr(par0) %*% C)
    }
    
    # optimize
    opt <- nlminb(x0, f, gr, control=list(trace=1))
    opt$objective
  }
  
  # Create hessian for optimizer
  
  # Robustify f
  f <- function(x0){
    # y <- try(f.optim(x0), silent=TRUE)
    y <- try_withWarningRecovery(f.optim(x0), silent=TRUE)
    if(is(y, "try-error")) y <- NA
    y
  }
  
  # test
  # x <- par.vals[!(names(par.vals) %in% parnames)]
  # new.f(x)
  # new.gr(x)
  
  # run optimizer
  # pars0 <- par.vals[!id]
  # opt <- nlminb(pars0, f, gr)
  
  # return
  f(x0)
}
