========
ClassMC
========
:ClassMC:	A Cosmological Parameter Estimation Code based on CLASS and Affine-Invariant MCMC algrithm
:Author:	Youhua Xu

Description
=============
This code is based on:
1) CLASS (--written by Julien Lesgourgues, https://github.com/lesgourg/class_public)
2) imcmc (--written by Youhua Xu, https://github.com/xyh-cosmo/imcmc/)

Prerequisite
=============
1) Planck Likelihood Code
2) imcmc (depends on OpenMPI & GSL)
3) Armadillo (a C++ linear algebra library, http://arma.sourceforge.net/)
4) WMAP7 / WMAP9 / Planck likelihoods 
5) If Multinest or PolyChord is used, then the multinest or polychord library should be installed properly.

Examples (added @2019-08-02)
=========
We provide several examples to demonstrate the usage of ClassMC:

1) LCDM_JD16_no_err: use mock JD16* samples with random error;
2) LCDM_JD16_1: use mock JD16* samples to constrain {H0, MB, omegabh2 omegach2}
3) LCDM_JDF_1: use mock JDF samples to constrain {H0, MB, omegabh2 omegach2}

The following are for the reconstruction of w(z):

4) DDE_JD16_no_err
5) DDE_JD16_1
6) DDE_JDF_1


The following commands (in the master folder) can be used to run the above 6 examples (correspondingly):

1) mpirun -np 4 bin/ClassMC LCDM_JD16_no_err/input.ini

2) mpirun -np 4 bin/ClassMC LCDM_JD16_1/input_1.ini

3) mpirun -np 4 bin/ClassMC LCDM_JDF_1/input_1.ini

4) mpirun -np 4 bin/ClassMC DDE_JD16_no_err/input.ini

5) mpirun -np 4 bin/ClassMC LCDM_JD16_1/input_1.ini

6) mpirun -np 4 bin/ClassMC LCDM_JDF_1/input_1.ini

Notes:
==========
We will kept updating this README and comments inside the example *ini files to make the useage of ClassMC more easier.
