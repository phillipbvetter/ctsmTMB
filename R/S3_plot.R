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
#' @param state.type a character vector either 'prior', 'posterior' or 'smoothed' determining what kind of states to plot.
#' 'prior' and 'posterior' are available for the Kalman filter methods, while 'smoothed' is available for the
#' Laplace methods (and is used automatically for those fits).
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
#'
#' # For fits obtained with the Laplace methods the smoothed states are plotted
#' plot(fit, type="states", state.type="smoothed", against.obs="y")
#' }
#' @returns a list of residual analysis / state plots (ggplots) depending on \code{type}.
#' @export
plot.ctsmTMB.fit = function(x,
                            print.plot = 1,
                            type = c("residuals", "states"),
                            state.type = c("prior","posterior","smoothed"),
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

    # The Laplace methods produce smoothed state estimates, which are stored
    # in the fit object as fit$states$mean$smoothed and fit$states$sd$smoothed
    if (method %in% laplace.methods) {
      if (state.type != "smoothed") {
        message("The Laplace methods produce smoothed state estimates - using state.type='smoothed'.")
        state.type <- "smoothed"
      }
      mean.df <- fit$states$mean$smoothed
      sd.df <- fit$states$sd$smoothed
      if (is.null(mean.df) || is.null(sd.df)) {
        stop("No smoothed states found in the fit object. Did you run 'estimate' with 'report=TRUE'?")
      }
      # construct the same 'mean_and_sd' data.frame structure used by the Kalman filters
      ms <- mean.df
      for (nam in state.names) ms[[paste0(nam, ".sd")]] <- sd.df[[nam]]
    } else {
      if (state.type == "smoothed") {
        stop("Smoothed states are only available for the Laplace methods - use state.type='prior' or 'posterior'.")
      }
      ms <- fit$states[[state.type]]$mean_and_sd
      if (is.null(ms)) {
        stop("No ", state.type, " states found in the fit object. Did you run 'estimate' with 'report=TRUE'?")
      }
    }

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
