#######################################################
# Print - S3 Method
#######################################################

#' Basic print of ctsmTMB objects
#' @param x  an object of class 'ctsmTMB'
#' @param ... additional arguments (not in use)
#'
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
#'
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
#' @param plot.series.i an integer determining with prediction series to plot, corresponding
#' to the index ".i" in the prediction data.frame.
#' @param type which series to plot (states, observations or both of those).
#' @param state.name the name of the state variable to plot
#' @param obs.name the name of the observation variable to plot
#' @param against.obs the name of the observations variable whose data series to plot against.
#' @param ggtheme a ggplot2 theme to use for the plots.
#' @param ... additional arguments
#' @examples
#' \dontrun{
#' data(Ornstein)
#' df <- Ornstein
#'
#' model <- create.Ornstein1D.model()
#'
#' # Given a prediction
#' pred <- model$predict(df, k.ahead=100)
#'
#' # plot the first series (starting from the posterior at i=0)
#' plot(pred)
#'
#' # We can also plot the observations instead, and compare with the data
#' plot(pred, type="obs", against.obs="y")
#'
#' # or plot a later series (startering from the posterior at i=50), for both obs and state
#' plot(pred, type="both", against.obs="y", plot.series.i=50)
#'
#' }
#' @returns A ggplot2 plot object.
#' @export
plot.ctsmTMB.pred = function(x,
                             y = NULL,
                             plot.series.i = 0,
                             type = c("states", "observations","both"),
                             state.name = attributes(pred)$state.names[1],
                             obs.name = attributes(pred)$obs.names[1],
                             against.obs = NULL,
                             ggtheme = getggplot2theme(),
                             ...) {

  type <- match.arg(type)
  pred <- x

  # Create plot
  colors <- getggplot2colors(3)
  colors <- colors[c(1,3,2)]
  names(colors) <- c(state.name, obs.name, paste0(against.obs,".data"))

  bool <- pred$states$mean_and_sd$i. %in% plot.series.i
  df1 <- df2 <- NULL
  if (type == "states" || type == "both") {
    df1 <- data.frame(
      t     = pred$states$mean_and_sd$t.j,
      x     = pred$states$mean_and_sd[[state.name]],
      x.sd  = pred$states$mean_and_sd[[paste0(state.name, ".sd")]],
      name = state.name
      # index = as.character(pred$states$mean_and_sd$i.)
    )[bool, ]
  }

  if (type == "observations" || type == "both") {
    df2 <- data.frame(
      t     = pred$observations$mean_and_sd$t.j,
      x     = pred$observations$mean_and_sd[[obs.name]],
      x.sd  = pred$observations$mean_and_sd[[paste0(obs.name, ".sd")]],
      name = obs.name
      # index = as.character(pred$observations$mean_and_sd$i.)
    )[bool, ]
  }
  df <- rbind(df1, df2)   # NULL entries are simply dropped

  # Add data
  obs.layer1 <- obs.layer2 <- NULL
  if (!is.null(against.obs)) {
    df.data <- data.frame(
      t = pred$data$t.j,
      y = pred$data[[against.obs]]
    )[bool, ]
    obs.layer1 <- ggplot2::geom_line(data=df.data,aes(x=t,y=y, color=paste0(against.obs,".data")),na.rm=T)
    obs.layer2 <- ggplot2::geom_point(data=df.data,aes(x=t,y=y, color=paste0(against.obs,".data")),na.rm=T)
  }

  p <- ggplot2::ggplot(data=df) +
    ggplot2::geom_ribbon(aes(x=t, ymin=x-2*x.sd, ymax=x+2*x.sd, fill=name, color=name), linetype="dashed", alpha=0.25) +
    ggplot2::geom_line(aes(x=t, y=x, color=name)) +
    ggplot2::geom_point(aes(x=t, y=x, color=name)) +
    obs.layer1 + obs.layer2 +
    scale_fill_manual(values=colors[-3],guide="none") +
    scale_color_manual(values=colors) +
    labs(color = NULL) +
    ggtheme

  # Show on plot window
  print(p)

  return(invisible(p))
}

#' This function creates residual plots for an estimated ctsmTMB object
#' @param x A R6 ctsmTMB fit object
#' @param print.plot an integer determining which of the produced plots to print in the plot window.
#' A value of 0 means not to plot anything.
#' @param type a character vector either 'residuals' or 'states' determining what
#' to plot.
#' @param state.type a character vector either 'prior' or 'posterior' determining what kind of states to plot.
#' @param against.obs name of an observation to plot state predictions against.
#' @param standardized.residuals logical whether or not the time-series and histogram plots show standardized residuals.
#' @param residual.burnin integer N - removes the first N entries in the residual vectors, including NAs.
#' @param residual.vs.obs.and.inputs the residual plots also include a new window
#' with time-series plots of residuals, associated observations and inputs
#' @param ggtheme ggplot2 theme to use for creating the ggplot.
#' @param time.series.ylims limits on the y-axis for residual time-series plot
#' @param hist.bins number of bins passed to \code{ggplot2::geom_histogram}
#' @param ... additional arguments
#' @examples
#' \dontrun{
#' # Given a fit from estimate i.e.
#' fit <- model$estimate(data)
#'
#' # We can generate a standard residual analysis plot by
#' plot(fit)
#'
#' # We can remove the first say 10 residuals of the series if they are off due
#' # to a bad initial prior state
#' plot(fit, residual.burnin=10)
#'
#' # We can also plot residuals against inputs and observations instead
#' plot(fit, residual.vs.obs.and.inputs=TRUE)
#'
#' # We can also plot the prior or posterior states using the 'type' and
#' # 'state.type' arguments.
#' plot(fit, type="states", state.type="prior")
#'
#' # We can also show observations on that plot
#' plot(fit, type="states", state.type="prior", against.obs="y")
#' }
#' @returns a list of residual analysis / state plots (ggplots) depending on \code{type}.
#' @export
plot.ctsmTMB.fit = function(x,
                            print.plot = 1,
                            type = c("residuals", "states"),
                            state.type = c("prior","posterior"),
                            against.obs=NULL,
                            standardized.residuals = TRUE,
                            residual.burnin = 0L,
                            residual.vs.obs.and.inputs = FALSE,
                            ggtheme=getggplot2theme(),
                            time.series.ylims = c(NA,NA),
                            hist.bins = 50,
                            ...) {

  type <- match.arg(type)
  state.type <- match.arg(state.type)

  fit <- x
  # remember fit$private has missing fields compared to general private
  private <- fit$private
  state.names <- private$names$states
  obs.names <- private$names$obs
  input.names <- private$names$inputs[-1] #remove 't' which is first input
  method <- private$method
  laplace.methods <- get.inference.algos("laplace")

  if (!ggplot2::is_theme(ggtheme)) {
    stop("'ggtheme' must be a ggplot2 theme, e.g. ggplot2::theme_minimal().")
  }

  if(!is.numeric(print.plot)){
    stop("print.plot must be a numeric value")
  }

  plots <- plots2 <- list()

  # residual plots
  ########################################
  if(type=="residuals"){

    # Are there any residuals?
    if(is.null(fit$residuals)){
      if(method %in% laplace.methods){
        stop("There are no residuals in the fit object. Did you remember to use 'laplace.residuals=TRUE' when you called 'estimate'?")
      } else {
        stop("There are no residuals in the fit object.")
      }
    }

    #check burnin
    stopifnot(is.numeric(residual.burnin), length(residual.burnin) == 1, residual.burnin >= 0)

    # color
    color <- "steelblue"
    ylab <- if(standardized.residuals) "Standardized Residuals" else "Raw Residuals"

    # define helper functions
    burn <- function(df, r.burn) if(r.burn == 0) df else df[-seq_len(r.burn),,drop=FALSE]

    correlation_plot <- function(e, title, type=c("acf","pacf")) {
      type <- match.arg(type)
      n.used <- sum(!is.na(e))
      if(type == "acf") {
        x <- stats::acf(e, na.action = stats::na.pass, plot = FALSE)
      } else {
        x <- stats::pacf(e, na.action = stats::na.pass, plot = FALSE)
      }
      lag <- if (type == "pacf") x$lag else x$lag[-1]   # acf includes lag 0, pacf does not
      val <- if (type == "pacf") x$acf else x$acf[-1]
      ci  <- 2 / sqrt(n.used)
      ggplot2::ggplot(data.frame(lag = c(lag), val = c(val)), ggplot2::aes(x = lag, ymin = 0, ymax = val)) +
        ggplot2::geom_linerange(color = color) +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::geom_hline(yintercept = c(-ci, ci), color = "blue", linetype = "dashed") +
        ggplot2::labs(title = title, x = "Lag", y = NULL) +
        ggtheme
    }

    if(residual.vs.obs.and.inputs && length(input.names) > 0) {
      # input vs time plot
      input.plots <- list()
      for(i in seq_along(input.names)) {
        nam <- input.names[i]
        df <- data.frame(
          t = private$data$t,
          u = private$data[[nam]]
        )
        input.plots[[i]] <-
          ggplot2::ggplot(data=df) +
          ggplot2::geom_line(ggplot2::aes(x=t, y=u), color=color) +
          ggplot2::labs(y=nam, x = "") +
          ggtheme
      }
      full.input.plot <- patchwork::wrap_plots(input.plots, ncol=1)
      # Add time to last plot (bottom most)
      full.input.plot[[length(full.input.plot)]]$labels$x <- "Time"
    }

    # RESIDUAL ANALYSIS PLOTS
    for (i in seq_along(obs.names)) {
      nam <- obs.names[i]

      # Create data frame
      df <- data.frame(
        t = fit$residuals$normalized$t,
        e = fit$residuals$residuals[[nam]],
        e.norm = fit$residuals$normalized[[nam]]
      )
      df$e0 <- if(standardized.residuals) df$e.norm else df$e
      # Burn away initial if requested
      df <- burn(df, residual.burnin)
      # Create another data.frame without NAs
      df.no.na <- df[!is.na(df$e), ]

      # residuals vs time
      # uses either standardized or not
      plot.res <- ggplot2::ggplot(data=df) +
        ggplot2::geom_line(ggplot2::aes(x=t, y=e0), color=color, na.rm=TRUE) +
        ggplot2::coord_cartesian(ylim=time.series.ylims) +
        ggplot2::labs(
          title = "Residuals" , subtitle = ylab, y = "", x = "Time"
        ) +
        ggtheme

      # histogram
      # uses either standardized or not
      plot.hist <-
        ggplot2::ggplot(data=df.no.na) +
        ggplot2::geom_histogram(ggplot2::aes(x=e0, y=ggplot2::after_stat(density)), bins=hist.bins, color="black", fill=color) +
        ggplot2::labs(title = "Histogram", subtitle=ylab, y = "",x = "") +
        ggtheme

      # QQ plot
      # always standardized
      plot.qq <-
        ggplot2::ggplot(data=df.no.na) +
        ggplot2::stat_qq(ggplot2::aes(sample=e.norm), color=color) +
        ggplot2::stat_qq_line(ggplot2::aes(sample=e.norm), lty="dashed") +
        ggplot2::labs(
          title = "Normal Q-Q Plot",
          y = "Sample Quantiles",
          x = "Theoretical Quantiles"
        ) +
        ggtheme

      # acf and pacf
      # always standardized
      plot.acf <- correlation_plot(df$e.norm,"Auto-Correlation", type="acf")
      plot.pacf <- correlation_plot(df$e.norm,"Partial Auto-Correlation", type="pacf")

      # cpgram
      # always standardized
      plot.cpgram <- ggfortify_ggcpgram(df$e.norm, colour=color) +
        ggplot2::labs(title = "Cumulative Periodogram", y = "", x = "Frequency") +
        ggtheme

      plots[[i]] <- patchwork::wrap_plots(plot.res,
                                          plot.hist,
                                          plot.qq,
                                          plot.cpgram,
                                          plot.acf,
                                          plot.pacf ,
                                          nrow=3) +
        patchwork::plot_annotation(title=paste("Residual Analysis for", nam),
                                   theme = ggtheme +
                                     ggplot2::theme(text=ggplot2::element_text(size=12, face="bold"))
        )

      # RESIDUAL, INPUTS AND OBSERVATIONS TIME SERIES PLOTS
      if(residual.vs.obs.and.inputs) {
        # Create data frame
        df <- data.frame(
          t = private$data$t,
          y = private$data[[nam]]
        )
        # time vs observations
        plot.obs <-
          ggplot2::ggplot(data=df) +
          ggplot2::geom_line(ggplot2::aes(x=t, y=y), color=color, na.rm=TRUE) +
          ggplot2::labs(y = nam, x = "") +
          ggtheme

        # save current residual vs input plot
        plot.res.no.label <- plot.res
        plot.res.no.label$labels$x = ""
        panels <- list(plot.res.no.label, plot.obs)
        if(length(input.names)>0) panels <- c(panels, list(full.input.plot))
        plots2[[i]] <- patchwork::wrap_plots(panels,
                                             ncol = 1,
                                             heights=c(1,1,length(input.names))
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
    ms <- fit$states[[state.type]]$mean_and_sd
    plot.obs <- !is.null(against.obs)
    if (plot.obs && !against.obs %in% obs.names) stop("'against.obs' must be one of: ", paste(obs.names, collapse = ", "))
    obs.layer <- NULL
    if (plot.obs) {
      obsdata <- data.frame(t = fit$private$data$t, y = fit$private$data[[against.obs]])
      obs.layer <- ggplot2::geom_point(data = obsdata,
                                       ggplot2::aes(x = t, y = y, color = against.obs),
                                       size = 0.5, na.rm=TRUE)
    }

    # We produce a ggplot for each state
    for (i in seq_len(private$dims$states)) {
      nam <- state.names[i]
      plot.title <- paste(capitalize_first(state.type),"estimates of", nam, if(plot.obs) paste("against",against.obs))
      # We must store the color in the data so that plots remember them on print
      plot.data <- data.frame(t=ms$t, x=ms[[nam]], x.sd=ms[[paste0(nam,".sd")]], col=nam)
      # Name color vector
      colors <- c("steelblue", if (plot.obs) "tomato")
      names(colors) <- c(nam, if (plot.obs) against.obs)

      # Create main plot
      plots[[i]] <- ggplot2::ggplot() +
        # Confidence interval
        ggplot2::geom_ribbon(data=plot.data,
                             ggplot2::aes(x=t, ymin=x-2*x.sd, ymax=x+2*x.sd),
                             fill="grey", alpha=0.6
        ) +
        # Mean estimate
        ggplot2::geom_line(data=plot.data,
                           ggplot2::aes(x=t, y=x, color=col)
        ) +
        # Observation points
        obs.layer +
        ggplot2::scale_color_manual(values=colors) +
        ggplot2::labs(title = plot.title, x="Time", color="", y = "") +
        ggtheme
    }
  }

  # Print to console
  if(print.plot != 0){
    temp.plot <- c(plots, plots2)[[print.plot]]
    print(temp.plot)
  }

  # return plot list
  return(invisible(c(plots, plots2)))
}

#' Performs full multi-dimensional profile likelihood calculations
#' @param fitted a ctsmTMB fit object
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
#' \dontrun{
#' # Given a model fit produced by
#' fit <- model$estimate(Ornstein)
#'
#' # We call profile with theta=NULL to use the standard range theta +- 3 std. devs.
#' prof <- profile(fit, parlist=list(theta=NULL))
#' }
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
    y <- try_with_warning_recovery(f.optim(x0), silent=TRUE)
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
#' \dontrun{
#' # Given a model fit produced by
#' fit <- model$estimate(Ornstein)
#'
#' # We calculate the profile likelihood
#' prof <- profile(fit, parlist=list(theta=NULL))
#'
#' # And then simply call plot on that object
#' plot(prof)
#' }
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

    # if (requireNamespace("geomtextpath", quietly = TRUE)) {
    #   # Use geomtextpath if available
    #   p <- ggplot2::ggplot() +
    #     geomtextpath::geom_textcontour(
    #       data = df,
    #       ggplot2::aes(
    #         x = .data[[par.names[1]]],
    #         y = .data[[par.names[2]]],
    #         z = .data$likelihood,
    #         label = labelfun(after_stat(level)),
    #         color = after_stat(level)
    #       ),
    #       breaks = chisquared_contours,
    #       linewidth = 0.6
    #     )
    # } else {
    # Fallback without labels if geomtextpath is not available
    # message("Please install 'geomtextpath' to get contour labels.")
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
    # }

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

