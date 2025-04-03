## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

### Author's Comments:

**Date: 2/4/2025:**

Hello CRAN team. This is the first submission of the 'ctsmTMB' package. The package has been repeatedly testing as follows, to make sure we follow the guidelines:

  1) Github Actions testing "using R-latest on Linux, Mac, and Windows, and using R-devel and R-oldrel on Linux."
      with the `usethis::use_github_action` command
      
  2) using  the `devtools` package via the commands: 
  - `devtools::check()`
  - `devtools::check_win_devel()`
  - `devtools::check_win_release()`
  - `devtools::check_mac_release()`

  3) Using a CRAN-like test via
  `withr::with_options(list(repos = c(CRAN = "https://cloud.r-project.org/")), {callr::default_repos() rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran")) })`

  4) Spellchecking and URL checking via
  `spelling::spell_check_package()` and `urlchecker::url_check()`

The testing has revealed the following NOTE, which have been addressed:

  1)  The *Suggested* 'RTMBode' package is not on **CRAN**. We have added to the **DESCRIPTION** file the r-universe repository in the *Additional_repositories* field. 
  The functionalities in 'ctsmTMB' that require the 'RTMBode' package have been wrapped in a `requireNamespace()` throwing an error if the package is not installed.
