
.onLoad <- function(lib, pkg) {
  
  # Check if RTMBode is installed:
  check <- requireNamespace("RTMBode",quietly=TRUE)
  if(!check){
    message(
      "Note: The RTMBode package for adaptive ODE-solvers is not installed.
The package can be installed from source using
  install.packages('RTMBode', repos = c('https://kaskr.r-universe.dev', 'https://cloud.r-project.org'))
     ") 
  }
  NULL
}
