# LTB-model_cpp
A c++ implementation of the concepts demonstrated in the python repository LTB-model. The purpose of this repository is to have a generic tool that solves the background and geodesic equations in the inhomogeneous Lemaître–Tolman  cosmological model. There are a few libraries written in Fortran and python that do exactly the above, but they rely on  exact analytical forms for the Mass function M(r), and Energy function E(r) which is a handicap.

To perform parameter constraints, the code relies on the python implementation (see http://dan.iel.fm/emcee/current/) of the affine invariant algorithm, see http://msp.org/camcos/2010/5-1/p04.xhtml. The algorithm performs much better then traditional samplers such as Metropolis-Hastings. For parameter estimation the speed bottleneck will not be because of the python interface but rather how fast one can calculate the likelihood function that one is maximizing or minimizing. You will also need Daniel et al.'s https://github.com/dfm/triangle.py to make the triangle plots.

A three parameter toy model for the dipole magnitude is used to test the code.

### Compiling and running
On linux do the following:
  * python setup.py build_ext --inplace

  * place cmb.o and scipy_tools.o produced from the first step in the working directory containing the rest of the code

  * compile the Makefile with the make command

### Defects and known issues
  * compiling with OpenMP for parallel programing crashes (while performing MCMC) on some systems.
  * error and exception handling needs to be implemented to conform to c++ standard

  
