## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

### Author's Comments:

**Date: 2/4/2025:**

Hello CRAN team. This is the first submission of the 'ctsmTMB' package. The package has been repeatedly testing as follows, to make sure we follow the guidelines:

  1) Github Actions testing "using R-latest on Linux, Mac, and Windows, and using R-devel and R-oldrel on Linux."
      with the `usethis::use_github_action` command
      
  2) using  the `devtools` package via the commands: 
    `devtools::check()`
    `devtools::check_win_devel()`
    `devtools::check_win_release()`
    `devtools::check_mac_release()`

  3) Using a CRAN-like test via
  `withr::with_options(list(repos = c(CRAN = "https://cloud.r-project.org/")), {callr::default_repos() rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran")) })`

The testing has revealed the following NOTE, which have been addressed:

  1) *Package suggested but not available for checking: ‘RTMBode’*.
  The *Suggested* 'RTMBode' package is not on **CRAN**. We have added to the **DESCRIPTION** file the r-universe repository in the *Additional_repositories* field. 
  We have also added an `.onLoad()` function which checks for the package via `requireNamespace()`. In addition, the functionalities in 'ctsmTMB' that require the 'RTMBode' package are similarly wrapped in a `requireNamespace()` throwing an error if the package is not installed.


                        
