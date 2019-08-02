Brief summary of the mock samples:

1) 	BAO
	.../6dF									[has conjugate, for JD16*]
	.../Dr12								[has conjugate, for JD16*]
	.../Lya									[has conjugate, for JD16*]
	.../MGS									[has conjugate, for JD16*]
	.../forecast							[no conjugate, for JDF*]
	.../forecasted_BAO_icovmat.txt			[for JDF*]
	.../forecasted_fractional_errs.txt

Note: the inv-covrariance matrix of Dr12 is "data/BAO_Zhao2016/DR12_tomoBAO_Zhao_invcov.txt"

2)	HST										[has conjugate, for JD16* & JDF*]

3)	Hz 										[has conjugate, for JD16*]

4)	Hz_noLya								[has conjugate, for JDF*]
This is basically the same as Hz, but without H(z=2.34) which is derived from the Lya-forest BAO measurement.

5)	JLA										[has conjugate, for JD16* and JDF*]
	jla_covmat.txt:							covariance matrix with light-cureve parameters fixed to their best fit values (assuming LCDM)
	jla_covmat_inv.txt						inverse of the above covariance matrix

6)	PLK15_dist_prior						[has conjugate, for JD16* and JDF*]

7)	WFIRST									[no conjugate, for JDF*]


##################################################################################

JD16* = BAO(Dr12 + Lya + MGS + 6dF) + HST + Hz + JLA + PLK15_dist_prior

JDF*  = BAO/forecast + HST + Hz_noLya + JLA + PLK15_dist_prior + WFIRST


##################################################################################
Each of the following mock samples has an extra realization that with error(s) set to zero:
BAO/6dF
BAO/Dr12
BAO/Lya
BAO/MGS

HST/
Hz/, Hz_noLya/
JLA/
PLK15_dist_prior/




