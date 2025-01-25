load('MouseSurveyData_Formatted.mat');

% Aproximate Bayesian Computation estimates:
% Pp.F0 = 13000;
% Pp.K = 26000; 
% Pp.r = .407;

% d16 estimates:
% Pd.F0 = 29500;
% Pd.K = 255000;
% Pd.r = 1.98;

%%%%%%%%%%%%%%%%%%%%%
% OLS estimates (using t = 3, 4, and 6 weeks)
Pp.F0 = 12000; %lb
Pp.K = 26000; %ub
Pp.r = 0.1346;

Pd.F0 = 26080;
Pd.K = 233000;
Pd.r = 0.7313;

palp_cutoff = 9;
old_cutoff = 1.5;

%%

mouse_model_p95_2C_15cutoff = compartmentalTumorProgression(myDataJan,Pd,Pp,old_cutoff);
mouse_model_p95_2C_90cutoff = compartmentalTumorProgression(myDataJan,Pd,Pp,palp_cutoff);
mouse_model_d16_2C_15cutoff = compartmentalTumorProgression(myDataJan,Pd,Pp,old_cutoff);
mouse_model_d16_2C_90cutoff = compartmentalTumorProgression(myDataJan,Pd,Pp,palp_cutoff);


%%
Nsamp = 1000;
tic
mouse_model_p95_2C_15cutoff.estimate_bootstrap_parameters('p95','2C',Nsamp);
%isplay(toc)
mouse_model_p95_2C_90cutoff.estimate_bootstrap_parameters('p95','2C',Nsamp);
%display(toc)
mouse_model_d16_2C_15cutoff.estimate_bootstrap_parameters('d16','2C',Nsamp);
%display(toc)
mouse_model_d16_2C_90cutoff.estimate_bootstrap_parameters('d16','2C',Nsamp);
toc

%%

mouse_model_p95_3C_15cutoff = compartmentalTumorProgression(myDataJan,Pd,Pp,old_cutoff);
mouse_model_p95_3C_90cutoff = compartmentalTumorProgression(myDataJan,Pd,Pp,palp_cutoff);
mouse_model_d16_3C_15cutoff = compartmentalTumorProgression(myDataJan,Pd,Pp,old_cutoff);
mouse_model_d16_3C_90cutoff = compartmentalTumorProgression(myDataJan,Pd,Pp,palp_cutoff);


%%
Nsamp = 1000;
tic
mouse_model_p95_3C_15cutoff.estimate_bootstrap_parameters('p95','3C_partialObs',Nsamp);
display(toc)
mouse_model_p95_3C_90cutoff.estimate_bootstrap_parameters('p95','3C_partialObs',Nsamp);
display(toc)
mouse_model_d16_3C_15cutoff.estimate_bootstrap_parameters('d16','3C_partialObs',Nsamp);
display(toc)
mouse_model_d16_3C_90cutoff.estimate_bootstrap_parameters('d16','3C_partialObs',Nsamp);
toc