// This include folder is not really needed but is used to allow compiling drift functions etc with cppXPtr when using
// devtools::load_all, since it looks for headers in sdeTMB/include.
// When installing the package normally this is taken care by with inst/include folder.

// Allows model formulations with the 'invlogit' function
double invlogit(double x){
  return 1/(1+exp(-x));
};
