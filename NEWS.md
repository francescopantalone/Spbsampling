# Spbsampling 1.1.0

* Added a `NEWS.md` file to track changes to the package.
* In the function pwd() and swd() the parameter nrepl is now set to default at 1.
* In the function pwd() and swd() the parameter niter is now set to default at 10.
* In the function stprod() and stsum() the parameter differ is now set to default at 1e-15.
* In the function stprod() and stsum() the parameter niter is now set to default at 1000.
* In the function stprod() and stsum() there is no print for each convergence step anymore. 
* The function stprod() is now implemented in C++ through the Armadillo library.
* The function stsum() is now implemented in C++ through the Armadillo library.
* The function hpwd() is now implemented in C++ through the Armadillo library.
* In all the functions there are now checks about the correct inputs.
* Correction of some typos along the guides.

