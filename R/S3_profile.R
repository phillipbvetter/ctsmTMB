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

