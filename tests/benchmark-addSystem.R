## Benchmark: time to add each successive state/observation to a model.
##
## Each state follows an independent OU process:
##   dx_i = -x_i * dt + dw_i
## with a direct observation:
##   y_i ~ x_i,  Var = sigma^2
##
## Run from the package root:
##   Rscript tests/benchmark-addSystem.R

library(ctsmTMB)

N_MAX <- 30   # maximum number of states to add

# ------------------------------------------------------------------
# Time each successive addSystem / addObs / setVariance call
# ------------------------------------------------------------------

model <- ctsmTMB$new()

t_addSystem  <- numeric(N_MAX)
t_addObs     <- numeric(N_MAX)
t_setVar     <- numeric(N_MAX)

cat(sprintf("%-6s  %12s  %12s  %12s  %12s\n",
            "state", "addSystem", "addObs", "setVariance", "cumulative"))
cat(strrep("-", 62), "\n")

cumulative <- 0

for (i in seq_len(N_MAX)) {

  sys_expr <- as.formula(paste0("dx", i, " ~ -x", i, " * dt + dw", i))
  obs_expr <- as.formula(paste0("y",  i, " ~ x", i))
  var_expr <- as.formula(paste0("y",  i, " ~ sigma^2"))

  t_addSystem[i] <- system.time(model$addSystem(sys_expr))["elapsed"]
  t_addObs[i]    <- system.time(model$addObs(obs_expr))["elapsed"]
  t_setVar[i]    <- system.time(model$setVariance(var_expr))["elapsed"]

  cumulative <- cumulative + t_addSystem[i] + t_addObs[i] + t_setVar[i]

  cat(sprintf("%-6d  %12.4f  %12.4f  %12.4f  %12.4f\n",
              i, t_addSystem[i], t_addObs[i], t_setVar[i], cumulative))
}

# ------------------------------------------------------------------
# Summary table
# ------------------------------------------------------------------

results <- data.frame(
  n           = seq_len(N_MAX),
  addSystem   = t_addSystem,
  addObs      = t_addObs,
  setVariance = t_setVar,
  total_call  = t_addSystem + t_addObs + t_setVar
)
results$cumulative <- cumsum(results$total_call)

cat("\n--- Summary ---\n")
cat(sprintf("Total time to build %d-state model: %.2f s\n",
            N_MAX, tail(results$cumulative, 1)))
cat(sprintf("Time for last addSystem call:       %.4f s\n",
            tail(t_addSystem, 1)))
cat(sprintf("Time for first addSystem call:      %.4f s\n",
            t_addSystem[1]))
cat(sprintf("Slowdown (last vs first):           %.1fx\n",
            tail(t_addSystem, 1) / t_addSystem[1]))

# ------------------------------------------------------------------
# Plot (if ggplot2 available)
# ------------------------------------------------------------------

if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("tidyr",   quietly = TRUE)) {

  library(ggplot2)
  library(tidyr)

  long <- tidyr::pivot_longer(
    results[, c("n", "addSystem", "addObs", "setVariance")],
    cols      = -n,
    names_to  = "call",
    values_to = "seconds"
  )

  p1 <- ggplot2::ggplot(long, ggplot2::aes(n, seconds, colour = call)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::labs(
      title = "Time per call vs. number of states already in model",
      x     = "State index (i-th addition)",
      y     = "Elapsed time (s)",
      colour = NULL
    ) +
    ggplot2::theme_bw()

  p2 <- ggplot2::ggplot(results, ggplot2::aes(n, cumulative)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::labs(
      title = "Cumulative build time vs. number of states",
      x     = "Number of states",
      y     = "Cumulative elapsed time (s)"
    ) +
    ggplot2::theme_bw()

  out_file <- file.path("tests", "benchmark-addSystem.png")
  ggplot2::ggsave(out_file,
                  plot   = patchwork::wrap_plots(p1, p2, ncol = 1),
                  width  = 8,
                  height = 8)
  cat(sprintf("\nPlot saved to %s\n", out_file))

} else {
  cat("\n(install ggplot2 + tidyr to generate a plot)\n")
}
