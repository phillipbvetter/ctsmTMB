# # ---------------------------------------------------------------------------
# # Author-only smoke tests for the summer school exercise scripts
# #
# # These tests source the exercise scripts from maintenance/summer_school/ to
# # verify that they still run end-to-end with the current version of ctsmTMB.
# # They only run for the author, since the scripts are long running (many model
# # compilations and estimations) and they are written for interactive use in
# # RStudio (plots, dev.new(), setwd()). Note that maintenance/ is excluded
# # from the built package (see .Rbuildignore), so the scripts are additionally
# # skipped whenever they are not available.
# # ---------------------------------------------------------------------------
#
# # Graphics functions called by the exercise scripts and their helper
# # functions. They are replaced by no-ops while sourcing, so that nothing is
# # plotted and no graphics devices are opened or closed during the test.
# graphics_functions <- c("plot", "points", "lines", "par", "abline", "polygon",
#                         "legend", "title", "axis", "mtext", "acf", "cpgram",
#                         "curve", "layout", "dev.new", "dev.off")
#
# # Temporarily replace graphics functions with no-ops in the global environment
# # (the exercise scripts source their helper functions there, and those helpers
# # plot too) and restore the original bindings afterwards.
# without_plotting <- function(expr) {
#   existed <- graphics_functions %in% ls(globalenv())
#   old <- mget(graphics_functions[existed], envir = globalenv())
#   noop <- function(...) invisible(NULL)
#   for (fn in graphics_functions) assign(fn, noop, envir = globalenv())
#   on.exit({
#     rm(list = graphics_functions[!existed], envir = globalenv())
#     list2env(old, envir = globalenv())
#   }, add = TRUE)
#   force(expr)
# }
#
# # Helper to source a summer school exercise script outside of RStudio.
# #
# # The scripts begin by clearing the workspace (rm(list=ls())), setting the
# # working directory to the location of the script through
# # rstudioapi::documentPath(), and opening the documentation (?ctsmTMB), which
# # only works when the script is open in RStudio. Instead, the working
# # directory is set here, those interactive-only lines are removed from the
# # script, and the remaining code is evaluated in a fresh environment. Use
# # `extra` to pre-define variables for parts of the exercises that are left
# # open for students.
# #
# # All output is silenced: silent=TRUE is appended to the $estimate() calls,
# # graphics functions are replaced by no-ops (also for the helper functions
# # that the scripts source into the global environment), and the remaining
# # console output (prints, messages, optimizer traces) is captured.
# source_exercise <- function(script, extra = list()) {
#   script <- normalizePath(script, mustWork = TRUE)
#   oldwd <- setwd(dirname(script))
#   on.exit(setwd(oldwd), add = TRUE)
#
#   lines <- readLines(script)
#   is_init <- grepl("^\\s*rm\\(list\\s*=\\s*ls\\(\\)\\)", lines) |
#     grepl("rstudioapi::documentPath", lines) |
#     grepl("^\\s*\\?ctsmTMB", lines)
#   lines <- lines[!is_init]
#
#   # Silence the estimations (the $estimate() calls in the scripts do not set
#   # the argument themselves)
#   is_est <- grepl("\\$estimate\\(", lines) & !grepl("\\bsilent\\b", lines)
#   lines[is_est] <- sub("\\)\\s*$", ", silent=TRUE)", lines[is_est])
#
#   tmp <- tempfile(fileext = ".R")
#   writeLines(lines, tmp)
#   on.exit(unlink(tmp), add = TRUE)
#
#   env <- new.env(parent = globalenv())
#   list2env(extra, envir = env)
#   suppressMessages(capture.output(
#     without_plotting(source(tmp, local = env))
#   ))
#   invisible(env)
# }
#
# summer_school <- file.path("..", "..", "maintenance", "summer_school")
#
# testthat::test_that("ex_greybox: q1_fit_and_validate.R runs end-to-end", {
#   # Author only test
#   testthat::skip_if_not(Sys.getenv("USER") == "pbrve", "Author test - the summer school exercise scripts are long running, interactive, and only available in the repository")
#   testthat::skip_if_not(dir.exists(summer_school), "Summer school exercise scripts not available")
#
#   env <- source_exercise(file.path(summer_school, "ex_greybox", "r", "q1_fit_and_validate.R"))
#   testthat::expect_true(is.finite(env$fit.rk4$nll))
# })
#
# testthat::test_that("ex_greybox: q2_q3_model_selection.R runs end-to-end", {
#   # Author only test
#   testthat::skip_if_not(Sys.getenv("USER") == "pbrve", "Author test - the summer school exercise scripts are long running, interactive, and only available in the repository")
#   testthat::skip_if_not(dir.exists(summer_school), "Summer school exercise scripts not available")
#
#   env <- source_exercise(file.path(summer_school, "ex_greybox", "r", "q2_q3_model_selection.R"))
#   testthat::expect_true(all(vapply(env$L, function(fit) is.finite(fit$nll), logical(1))))
# })
#
# testthat::test_that("ex_greybox_part_2: q1_system_noise_levels.R runs end-to-end", {
#   # Author only test
#   testthat::skip_if_not(Sys.getenv("USER") == "pbrve", "Author test - the summer school exercise scripts are long running, interactive, and only available in the repository")
#   testthat::skip_if_not(dir.exists(summer_school), "Summer school exercise scripts not available")
#
#   env <- source_exercise(file.path(summer_school, "ex_greybox_part_2", "r", "q1_system_noise_levels.R"))
#   testthat::expect_true(is.finite(env$fit2$nll))
# })
#
# testthat::test_that("ex_greybox_part_2: q2_splined_ga_value.R runs end-to-end", {
#   # Author only test
#   testthat::skip_if_not(Sys.getenv("USER") == "pbrve", "Author test - the summer school exercise scripts are long running, interactive, and only available in the repository")
#   testthat::skip_if_not(dir.exists(summer_school), "Summer school exercise scripts not available")
#
#   # The script leaves the boundary angles for the gA-spline as an exercise
#   # for students (azimuth_boundary <- c(?, ?)). Set them here to the values
#   # implied by the plotting code at the end of the script (Boundary.knots=c(1,5)
#   # in radians, since aoiSunAzimuthDeg() returns radians).
#   env <- source_exercise(
#     file.path(summer_school, "ex_greybox_part_2", "r", "q2_splined_ga_value.R"),
#     extra = list(azimuth_boundary = c(1, 5))
#   )
#   testthat::expect_true(is.finite(env$fit2$nll))
# })
