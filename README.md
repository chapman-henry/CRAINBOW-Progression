Code and data files for "Nonlinear progression during the occult transition established cancer lethality"
Using mouse survey data, we estimate parameter values of an ODE model, plot the trajectories of the model with error bars,
and forward simulate different synthetic scenarios to determine difference between 2 HER2+ isoforms 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
Data files:

FSL_bootstrap_exper_FbyOLS.mat -- output of param_est_script.m for 3 compartment FSL model  Includes 4 objects of class compartmentalTumorProgression
	mouse_model_d16_2C_15cutoff (obj)
	mouse_model_d16_2C_90cutoff (obj)
	mouse_model_p95_2C_15cutoff (obj)
	mouse_model_p95_2C_90cutoff (obj)

FOSL_bootstrap_exper_FbyOLS.mat -- output of param_est_script. NOTE! the objects of class compartmentalTumorProgression have the same names as above
	mouse_model_d16_2C_15cutoff (obj)
	mouse_model_d16_2C_90cutoff (obj)
	mouse_model_p95_2C_15cutoff (obj)
	mouse_model_p95_2C_90cutoff (obj)
	Nsamp (int) -- number of samples for bootstrap == 1000
	old_cutoff (dbl) -- size threshold for S / L tumors == 1.5
	palp_cutoff (dbl) -- size threshold for S / L tumors == 9
	Pd (struct) --- expect value of d16 field parameters, 
	taken from OLS estimation of logistic equation from experimental data
	Fields:
	--- F0 == 26080
	--- K == 233000
	--- r == 0.7313
	Pp (struct) --- expect value of p95 field parameters, 
	taken from OLS estimation of logistic equation from experimental data
	Fields:
	--- F0 == 12000
	--- K == 26000
	--- r == 0.1346
 
bootstrap_noOccult_Jan2025.mat --- parameter estimations for 1000 bootstrap samples using 
both d16 and p95 specific estimates of F0, K, and r using FSL model
	param_results (struct) with 2 fields
	--- d16_params (table) 1000 x 7, cols K, T0, r, P1, P2, LSQ_boot, LSQ_pop 
	--- p95_params (table) same


MouseSurveyData_Formatted.mat --- eponymous table. Josh Ginzels mouse survey data.  1917 x 7 with columns:
	--- ID (categorical) unique id of mouse
	--- Age (dbl) weeks of age at sacrifice
	--- Gland (categorical) region of origin, 3 on each side of butterflied mouse
	--- Genotype : isoform, wt, d16, p95 
	--- Area (int) tumor area in sq. micron
	--- Area (dbl) tumor area in mm^2
	--- Mm2 (boolean) S = 0 , L = 1, with 1.5 mm^2 cutoff

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
Class definition:

compartmentalTumorProgression.m methods for estimating parameters for compartmental ODE models describing tumor progression in CRAINBOW HER2+ mouse survey data

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
Scripts:

param_est_script.m --- instantiates and calls methods of class compartmentalTumorProgression(). 
Produces FSL_bootstrap_exper_FbyOLS.mat and FOSL_bootstrap_exper_FbyOLS.mat

mouse_survey_publication_figs.m --- makes figures 3 and 4 

FSL_expected_tumors.m --- forward simulates FSL model as a function of K == F0, data for figure 5
