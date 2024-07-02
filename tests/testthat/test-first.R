
# This is a testthat script for automatically testing that the functions in
# sdeTMB are working as intended.

obj = sdeTMB$new()
testthat::expect_s3_class(obj,class=c("sdeTMB","R6"))

# system equations
testthat::expect_error(obj$add_systems(x ~ y))
testthat::expect_error(obj$add_systems(x123 ~ y))
testthat::expect_error(obj$add_systems(dx ~ 1 * dt + dw))
testthat::expect_error(obj$add_systems(dx ~ y123))


testthat::expect_no_error(obj$add_systems(dx ~ dw))
testthat::expect_no_error(check_system_eqs(dx ~ y*dt + dw))
