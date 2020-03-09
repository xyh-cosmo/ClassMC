# ClassMC

This is a Cosmological Parameter Estimation Code based on CLASS and Affine-Invariant MCMC algorithm.

# Dependence

This code is based on:

1) CLASS (written by Julien Lesgourgues, https://github.com/lesgourg/class_public)

2) imcmc (written by Youhua Xu, https://github.com/xyh-cosmo/imcmc/)

# Other prerequisites
1) Armadillo (a C++ linear algebra library, http://arma.sourceforge.net/);
2) WMAP7 / WMAP9 likelihoods (https://lambda.gsfc.nasa.gov/product/map/current/likelihood_get.cfm/ ; we have included these likelihood codes inside ClassMC); 
3) Planck likelihood (https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/CMB_spectrum_%26_Likelihood_Code#2018_Likelihood);
4) If one wants to use the Multinest or PolyChord sampler, then he/she needs to install the multinest (https://github.com/farhanferoz/MultiNest) or polychord (https://github.com/PolyChord/PolyChordLite) library properly (we have provided two interfaces to Multinest and PolyChord, see the C++ source files inside the main/ folder).

# Examples (added @ 2019-08-02)

To demonstrate the usage of ClassMC, we provide several examples (extracted from the works in https://arxiv.org/abs/1907.13298), each of which is grouped into a single folder:

1) LCDM_JD16_no_err/: use mock JD16* samples with random error;
2) LCDM_JD16_1/: use mock JD16* samples to constrain {H0, MB, omegabh2 omegach2};
3) LCDM_JDF_1/: use mock JDF samples to constrain {H0, MB, omegabh2 omegach2};

and the following are for the reconstruction of w(z) from JD16* without/with random errors and from  JDF16:

4) DDE_JD16_no_err/
5) DDE_JD16_1/
6) DDE_JDF_1/


The following commands (in the master folder) can be used to run the above 6 examples (correspondingly):
```bash
mpirun -np 4 bin/ClassMC LCDM_JD16_no_err/input.ini
```
```bash
mpirun -np 4 bin/ClassMC LCDM_JD16_1/input_1.ini
```
```bash
mpirun -np 4 bin/ClassMC LCDM_JDF_1/input_1.ini

```
```bash
mpirun -np 4 bin/ClassMC DDE_JD16_no_err/input.ini
```
```bash
mpirun -np 4 bin/ClassMC DDE_JD16_1/input_1.ini
```
```bash
mpirun -np 4 bin/ClassMC DDE_JDF_1/input_1.ini
```

# TODO

- [ ] Updating this README as well as comments inside the example ini files to make the useage of ClassMC easier;
- [ ] Writing a detailed document about the code structure of ClassMC as well as modifications to CLASS;
- [ ] Adding support for Planck 2018 likelihoods (https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/CMB_spectrum_%26_Likelihood_Code#2018_Likelihood);
- [ ] Updating to newer version of CLASS.

# Problems & Solutions

- [x] Recently, I eccountered the following compilation problem after upgrading GSL to the latest version (v2.6):
```
In file included from ../source/datasets/Data_SNE_JLA.cc:21:
/usr/include/x86_64-linux-gnu/cblas.h:69:13: error: conflicting declaration of C function ‘int cblas_izamax(int, const void*, int)’
   69 | CBLAS_INDEX cblas_izamax(const int N, const void   *X, const int incX);
      |             ^~~~~~~~~~~~
In file included from /usr/local/include/gsl/gsl_blas_types.h:28,
                 from /usr/local/include/gsl/gsl_matrix_complex_long_double.h:29,
                 from /usr/local/include/gsl/gsl_matrix.h:4,
                 from /usr/local/include/gsl/gsl_randist.h:24,
                 from /usr/local/include/imcmc/imcmc.hpp:29,
                 from ../include/jla.h:15,
                 from ../source/datasets/Data_SNE_JLA.cc:3:
/usr/local/include/gsl/gsl_cblas.h:102:13: note: previous declaration ‘size_t cblas_izamax(int, const void*, int)’
  102 | CBLAS_INDEX cblas_izamax(const int N, const void   *X, const int incX);
      |             ^~~~~~~~~~~~

```
To resolve this problem, I changed a little bit inside Data_SNE_JLA.cc in line 22:
```
//#include <cblas.h>
#include <gsl/gsl_cblas.h>
```
