# ClassMC

This is a Cosmological Parameter Estimation Code based on CLASS and Affine-Invariant MCMC algorithm.

# Dependence

This code is based on:

1) CLASS (written by Julien Lesgourgues, https://github.com/lesgourg/class_public)

2) imcmc (written by Youhua Xu, https://github.com/xyh-cosmo/imcmc/)

# Other prerequisites
1) Armadillo (a C++ linear algebra library, http://arma.sourceforge.net/);
2) WMAP7 / WMAP9 (https://lambda.gsfc.nasa.gov/product/map/current/likelihood_get.cfm/); 
3) Planck likelihood (https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/CMB_spectrum_%26_Likelihood_Code#2018_Likelihood);
4) If one wants to use the Multinest or PolyChord sampler, then he/she needs to install the multinest (https://github.com/farhanferoz/MultiNest) or polychord (https://github.com/PolyChord/PolyChordLite) library properly (we have provided two interfaces to Multinest and PolyChord, see the C++ source files inside the main/ folder).

# Examples (added @ 2019-08-02)

To demonstrate the usage of ClassMC, we provide several examples (extracted from the works in https://arxiv.org/abs/1907.13298),

each of which is grouped into a folder:

1) LCDM_JD16_no_err/: use mock JD16* samples with random error;
2) LCDM_JD16_1/: use mock JD16* samples to constrain {H0, MB, omegabh2 omegach2}
3) LCDM_JDF_1/: use mock JDF samples to constrain {H0, MB, omegabh2 omegach2}

The following are for the reconstruction of w(z):

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
