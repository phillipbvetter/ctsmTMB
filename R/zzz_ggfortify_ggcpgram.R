# -------------------------
# DISCLAIMER!!
# These functions are copy-pasted from ggfortify::ggcpgram !!!!
# is done in order to avoid the package dependency

# Link to function on Github:
# https://github.com/sinhrks/ggfortify/blob/33652f781a7046e55468599e3b992a01dfcc5165/R/tslib.R#L227
# Link to package CRAN package
# https://cran.r-project.org/web/packages/ggfortify/index.html
# -------------------------


################################################
# MAIN FUNCTION TO CALL FOR CUMULATIVE PERIDOGRAM
################################################
ggfortify_ggcpgram <- function (ts, taper = 0.1, colour = "#000000", linetype = "solid", 
                                conf.int = TRUE, conf.int.colour = "#0000FF", conf.int.linetype = "dashed", 
                                conf.int.fill = NULL, conf.int.alpha = 0.3) 
{
  is.univariate(ts)
  x <- as.vector(ts)
  x <- x[!is.na(x)]
  x <- stats::spec.taper(scale(x, TRUE, FALSE), p = taper)
  y <- Mod(stats::fft(x))^2/length(x)
  y[1L] <- 0
  n <- length(x)
  x <- (0:(n/2)) * stats::frequency(ts)/n
  if (length(x)%%2 == 0) {
    n <- length(x) - 1
    y <- y[1L:n]
    x <- x[1L:n]
  }
  else y <- y[seq_along(x)]
  xm <- stats::frequency(ts)/2
  mp <- length(x) - 1
  crit <- 1.358/(sqrt(mp) + 0.12 + 0.11/sqrt(mp))
  d <- data.frame(x = x, y = cumsum(y)/sum(y), upper = 1/xm * 
                    x + crit, lower = 1/xm * x - crit)
  p <- ggplot2::ggplot(data = d,
                       # mapping = ggplot2::aes_string(x = "x", y = "y")
                       mapping = ggplot2::aes(x=x, y=y)
                       ) + 
    geom_line(colour = colour, linetype = linetype) + 
    ggplot2::scale_x_continuous(name = "", limits = c(0, 
                                                      xm)) + ggplot2::scale_y_continuous(name = "", limits = c(0, 
                                                                                                               1))
  p <- plot_confint(p = p, data = d, conf.int = conf.int, conf.int.colour = conf.int.colour, 
                    conf.int.linetype = conf.int.linetype, conf.int.fill = conf.int.fill, 
                    conf.int.alpha = conf.int.alpha)
  p
}

################################################
# HELPER FUNCTION
################################################
is.univariate <- function (data, raise = TRUE) 
{
  if (ncol(as.matrix(data)) > 1) {
    if (raise) {
      stop("data must be univariate time series")
    }
    else {
      return(FALSE)
    }
  }
  return(TRUE)
}

################################################
# HELPER FUNCTION
################################################
plot_confint <- function (p, data = NULL, lower = "lower", upper = "upper", conf.int = TRUE, 
                          conf.int.geom = "line", conf.int.group = NULL, conf.int.colour = "#0000FF", 
                          conf.int.linetype = "none", conf.int.fill = "#000000", conf.int.alpha = 0.3) 
{
  if (missing(conf.int) && (!missing(conf.int.colour) || !missing(conf.int.linetype) || 
                            !missing(conf.int.fill) || !missing(conf.int.alpha))) {
    conf.int <- TRUE
  }
  if (is.null(data)) {
    stop("Internal Error: 'data' must be provided to plot_confint")
  }
  # if (conf.int.geom == "step") {
  # ribbon_func <- geom_confint
  # line_func <- geom_step
  # }
  # else {
  ribbon_func <- geom_ribbon
  line_func <- geom_line
  # }
  if (conf.int) {
    if (!is.null(conf.int.fill)) {
      p <- p + geom_factory(ribbon_func, data, ymin = lower, 
                            ymax = upper, group = conf.int.group, fill = conf.int.fill, 
                            alpha = conf.int.alpha, na.rm = TRUE)
    }
    if (conf.int.linetype != "none") {
      p <- p + geom_factory(line_func, data, y = lower, 
                            group = conf.int.group, colour = conf.int.colour, 
                            linetype = conf.int.linetype, na.rm = TRUE)
      p <- p + geom_factory(line_func, data, y = upper, 
                            group = conf.int.group, colour = conf.int.colour, 
                            linetype = conf.int.linetype, na.rm = TRUE)
    }
  }
  p
}

################################################
# HELPER FUNCTION
################################################
# geom_factory <- function (geomfunc, data = NULL, position = NULL, ...) 
# {
#   mapping <- list()
#   option <- list()
#   columns <- colnames(data)
#   for (key in names(list(...))) {
#     value <- list(...)[[key]]
#     if (is.null(value)) {
#     }
#     else if (!(is.vector(value) && length(value) > 1L) && 
#              value %in% columns) {
#       mapping[[key]] <- value
#     }
#     else {
#       option[[key]] <- value
#     }
#   }
#   if (!is.null(data)) {
#     option[["data"]] <- data
#   }
#   if (!is.null(position)) {
#     option[["position"]] <- position
#   }
#   option[["mapping"]] <- do.call(ggplot2::aes_string, mapping)
#   print(option[["mapping"]])
#   return(do.call(geomfunc, option))
# }

# re-write using rlang::sym to avoid depreciated aes_string
geom_factory <- function (geomfunc, data = NULL, position = NULL, ...) 
{
  mapping <- list()
  option <- list()
  columns <- colnames(data)
  for (key in names(list(...))) {
    value <- list(...)[[key]]
    if (is.null(value)) {
      # skip
    } else if (!(is.vector(value) && length(value) > 1L) && value %in% columns) {
      # capture column references properly with rlang::sym
      mapping[[key]] <- rlang::sym(value)
    } else {
      option[[key]] <- value
    }
  }
  if (!is.null(data)) {
    option[["data"]] <- data
  }
  if (!is.null(position)) {
    option[["position"]] <- position
  }
  option[["mapping"]] <- do.call(ggplot2::aes, mapping)   # use aes() not aes_string()
  return(do.call(geomfunc, option))
}
