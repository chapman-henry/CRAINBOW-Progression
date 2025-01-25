% Figures
% 20 Mar 2024

load('FSL_bootstrap_exper_FbyOLS.mat')

% Variables of interest
% mouse_model_d16_2C_15cutoff, mouse_model_p95_2C_15cutoff,
% mouse_model_d16_2C_90cutoff, mouse_model_p95_2C_90cutoff
%   all of class compartmentalTumorProgression
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % FSL Model trajectories with error bars
% Fig 3 C,D,E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tsteps = 3:0.1:26;
T = length(tsteps);
N = 1000;
state_dim = 3;
boot_traj_p95 = zeros(N,T,state_dim);
boot_traj_d16 = zeros(N,T,state_dim);

for n = 1:N
    Pd.K = mouse_model_d16_2C_15cutoff.param_estimate(n,1);
    Pd.F0 = mouse_model_d16_2C_15cutoff.param_estimate(n,2);
    Pd.r = mouse_model_d16_2C_15cutoff.param_estimate(n,3);
    Pd.p1 = mouse_model_d16_2C_15cutoff.param_estimate(n,4);
    Pd.p2 = mouse_model_d16_2C_15cutoff.param_estimate(n,5);

    Pp.K = mouse_model_p95_2C_15cutoff.param_estimate(n,1);
    Pp.F0 = mouse_model_p95_2C_15cutoff.param_estimate(n,2);
    Pp.r = mouse_model_p95_2C_15cutoff.param_estimate(n,3);
    Pp.p1 = mouse_model_p95_2C_15cutoff.param_estimate(n,4);
    Pp.p2 = mouse_model_p95_2C_15cutoff.param_estimate(n,5);


    boot_traj_d16(n,:,:) = mouse_model_d16_2C_15cutoff.F2C_all(tsteps,Pd);
    boot_traj_p95(n,:,:) =  mouse_model_d16_2C_15cutoff.F2C_all(tsteps,Pp);
end
%%
S_dm = reshape(median(boot_traj_d16),[T,state_dim]);
S_pm = reshape(median(boot_traj_p95),[T,state_dim]);

S_d05 = reshape(prctile(boot_traj_d16,5),[T,state_dim]);
S_p05 = reshape(prctile(boot_traj_p95,5),[T,state_dim]);

S_d95 = reshape(prctile(boot_traj_d16,95),[T,state_dim]);
S_p95 = reshape(prctile(boot_traj_p95,95),[T,state_dim]);


%%
%  F(t) with error
t2 = [tsteps, fliplr(tsteps)];

% % % % 
f = figure; % 
hold on
plot(tsteps,S_dm(:,1),'k','LineWidth',1)
plot(tsteps,S_pm(:,1),'k','LineWidth',1)
SdCI = [S_d05(:,1); flipud(S_d95(:,1))];
SpCI = [S_p05(:,1); flipud(S_p95(:,1))];
fill(t2,SdCI,'y','FaceAlpha', 0.4)
fill(t2,SpCI,'m','FaceAlpha', 0.4)
xlim([0 25])
title('F, c=1.5','FontSize',16)
xlabel('Time (weeks)')
ylabel('count')
exportgraphics(f,'figures\FSL_FvsTime_c15_Fosl.pdf','Resolution',600)


% % % %
% S(t) with error
f = figure;
hold on
plot(tsteps,S_dm(:,2),'k','LineWidth',1)
plot(tsteps,S_pm(:,2),'k','LineWidth',1)
SdCI = [S_d05(:,2); flipud(S_d95(:,2))];
SpCI = [S_p05(:,2); flipud(S_p95(:,2))];
fill(t2,SdCI,'y','FaceAlpha', 0.4)
fill(t2,SpCI,'m','FaceAlpha', 0.4)
xlim([0 25])
title('S, c=1.5','FontSize',16)
xlabel('Time (weeks)')
ylabel('count')
exportgraphics(f,'figures\FSL_SvsTime_c15_Fosl.pdf','Resolution',600)


% % % % % %
% L(t) with error
f = figure;
hold on
plot(tsteps,S_dm(:,3),'k','LineWidth',1)
plot(tsteps,S_pm(:,3),'k','LineWidth',1)
SdCI = [S_d05(:,3); flipud(S_d95(:,3))];
SpCI = [S_p05(:,3); flipud(S_p95(:,3))];
fill(t2,SdCI,'y','FaceAlpha', 0.4)
fill(t2,SpCI,'m','FaceAlpha', 0.4)
xlim([0 25])
title('L, c=1.5','FontSize',16)
xlabel('Time (weeks)')
ylabel('count')
exportgraphics(f,'figures\FSL_LvsTime_c15_Fosl.pdf','Resolution',600)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter estimates
% % % Bootstrap histograms
% Fig 3 G,H 
%  FSL palpable large tumors p1 bootstrap, p2 bootstrap 
plot_params_d16 = sort(mouse_model_d16_2C_90cutoff.param_estimate(:,4:5));
plot_params_p95 = sort(mouse_model_p95_2C_90cutoff.param_estimate(:,4:5));


f = figure;
hold on
histogram(plot_params_d16(1:999,1),20,'FaceColor',[1 1 0],'Normalization','probability')
histogram(plot_params_p95(:,1),20,'FaceColor',[1 0 1],'Normalization','probability')
legend('d16','p95')
xlabel('Bootstrap estimates of p_1')
ylabel('Frequency')
title('FSL c = 9')
%exportgraphics(f,'figures\FSL_p1bootstrap_c90_Fosl.pdf','Resolution',600)

f = figure;
hold on
histogram(plot_params_d16(1:999,2),20,'FaceColor',[1 1 0],'Normalization','probability')
histogram(mouse_model_p95_2C_90cutoff.param_estimate(:,5),20,'FaceColor',[1 0 1],'Normalization','probability')
legend('d16','p95')
xlabel('Bootstrap estimates of p_2')
ylabel('Frequency')
title('FSL c = 9')
%exportgraphics(f,'figures\FSL_p2bootstrap_c90_Fosl.pdf','Resolution',600)


% Cutoff > 1.5
plot_params_d16 = sort(mouse_model_d16_2C_15cutoff.param_estimate(:,4:5));
plot_params_p95 = (mouse_model_p95_2C_15cutoff.param_estimate(:,4:5));

f = figure;
hold on
histogram(plot_params_d16(1:999,1),20,'FaceColor',[1 1 0],'Normalization','probability')
histogram(plot_params_p95(:,1),20,'FaceColor',[1 0 1],'Normalization','probability')
legend('d16','p95')
xlabel('Bootstrap estimates of p_1')
ylabel('Frequency')
title('FSL c = 1.5')
%exportgraphics(f,'figures\FSL_p1bootstrap_c15_Fosl.pdf','Resolution',600)


r = figure;
hold on
histogram(plot_params_d16(1:999,2),20,'FaceColor',[1 1 0],'Normalization','probability')
histogram(mouse_model_p95_2C_15cutoff.param_estimate(:,5),20,'FaceColor',[1 0 1],'Normalization','probability')
legend('d16','p95')
xlabel('Bootstrap estimates of p_2')
ylabel('Frequency')
title('FSL c = 1.5')
%exportgraphics(f,'figures\FSL_p2bootstrap_c15_Fosl.pdf','Resolution',600)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %  FOSL model, trajectories with error bars
% Fig 4 G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('FOSL_bootstrap_exper_FbyOLS.mat')

%%
tsteps = 3:0.1:26;
T = length(tsteps);
N = 1000;
state_dim = 4;
boot_traj_p95 = zeros(N,T,state_dim);
boot_traj_d16 = zeros(N,T,state_dim);

for n = 1:N
    Pd.K = mouse_model_d16_3C_90cutoff.param_estimate(n,1);
    Pd.F0 = mouse_model_d16_3C_90cutoff.param_estimate(n,2);
    Pd.r = mouse_model_d16_3C_90cutoff.param_estimate(n,3);
    Pd.p0 = mouse_model_d16_3C_90cutoff.param_estimate(n,4);
    Pd.p1 = mouse_model_d16_3C_90cutoff.param_estimate(n,5);
    Pd.p2 = mouse_model_d16_3C_90cutoff.param_estimate(n,6);

    Pp.K = mouse_model_p95_3C_90cutoff.param_estimate(n,1);
    Pp.F0 = mouse_model_p95_3C_90cutoff.param_estimate(n,2);
    Pp.r = mouse_model_p95_3C_90cutoff.param_estimate(n,3);
    Pp.p0 = mouse_model_p95_3C_90cutoff.param_estimate(n,4);
    Pp.p1 = mouse_model_p95_3C_90cutoff.param_estimate(n,5);
    Pp.p2 = mouse_model_p95_3C_90cutoff.param_estimate(n,6);


    boot_traj_d16(n,:,:) = mouse_model_d16_3C_90cutoff.F3C_all(tsteps,Pd);
    boot_traj_p95(n,:,:) = mouse_model_d16_3C_90cutoff.F3C_all(tsteps,Pp);
end
%%
S_dm = reshape(median(boot_traj_d16),[T,state_dim]);
S_pm = reshape(median(boot_traj_p95),[T,state_dim]);

S_d05 = reshape(prctile(boot_traj_d16,5),[T,state_dim]);
S_p05 = reshape(prctile(boot_traj_p95,5),[T,state_dim]);

S_d95 = reshape(prctile(boot_traj_d16,95),[T,state_dim]);
S_p95 = reshape(prctile(boot_traj_p95,95),[T,state_dim]);
%%
%  F(t) with error
t2 = [tsteps, fliplr(tsteps)];

% % % % 
f = figure; % 
hold on
plot(tsteps,S_dm(:,1),'k','LineWidth',1)
plot(tsteps,S_pm(:,1),'k','LineWidth',1)
SdCI = [S_d05(:,1); flipud(S_d95(:,1))];
SpCI = [S_p05(:,1); flipud(S_p95(:,1))];
fill(t2,SdCI,'y','FaceAlpha', 0.4)
fill(t2,SpCI,'m','FaceAlpha', 0.4)
xlim([0 25])
title('F, FOSL, c = 9','FontSize',16)
xlabel('Time (weeks)')
ylabel('count')
%exportgraphics(f,'figures\FOSL_FvsTime_c90_Fols.pdf','Resolution',600)

%  O(t) with error
f = figure;
hold on
plot(tsteps,S_dm(:,2),'k','LineWidth',1)
plot(tsteps,S_pm(:,2),'k','LineWidth',1)
SdCI = [S_d05(:,2); flipud(S_d95(:,2))];
SpCI = [S_p05(:,2); flipud(S_p95(:,2))];
fill(t2,SdCI,'y','FaceAlpha', 0.4)
fill(t2,SpCI,'m','FaceAlpha', 0.4)
ylim([0 30])
xlim([0 25])
title('O, FOSL, c = 9','FontSize',16)
xlabel('Time (weeks)')
ylabel('count')
%exportgraphics(f,'figures\FOSL_OvsTime_c90_Fols.pdf','Resolution',600)


%  S(t) with error
f = figure;
hold on
plot(tsteps,S_dm(:,3),'k','LineWidth',1)
plot(tsteps,S_pm(:,3),'k','LineWidth',1)
SdCI = [S_d05(:,3); flipud(S_d95(:,3))];
SpCI = [S_p05(:,3); flipud(S_p95(:,3))];
fill(t2,SdCI,'y','FaceAlpha', 0.4)
fill(t2,SpCI,'m','FaceAlpha', 0.4)
xlim([0 25])
title('S, FOSL, c = 9','FontSize',16)
xlabel('Time (weeks)')
ylabel('count')
%exportgraphics(f,'figures\FOSL_SvsTime_c90_Fols.pdf','Resolution',600)



%  L(t) with error
f = figure;
hold on
plot(tsteps,S_dm(:,4),'k','LineWidth',1)
plot(tsteps,S_pm(:,4),'k','LineWidth',1)
SdCI = [S_d05(:,4); flipud(S_d95(:,4))];
SpCI = [S_p05(:,4); flipud(S_p95(:,4))];
fill(t2,SdCI,'y','FaceAlpha', 0.4)
fill(t2,SpCI,'m','FaceAlpha', 0.4)
xlim([0 25])
title('L, FOSL, c = 9','FontSize',16)
xlabel('Time (weeks)')
ylabel('count')
%exportgraphics(f,'figures\FOSL_LvsTime_c90_Fols.pdf','Resolution',600)

%% FOSL Parameter Estimates, bootstrap histograms
% Fig 4 D,E,F
%  

plot_params_d16 = sort(mouse_model_d16_3C_90cutoff.param_estimate(:,4:6));
plot_params_p95 = sort(mouse_model_p95_3C_90cutoff.param_estimate(:,4:6));


% p0
f = figure;
hold on
histogram(plot_params_d16(26:975,1),20,'FaceColor',[1 1 0],'Normalization','probability')
histogram(plot_params_p95(26:975,1),20,'FaceColor',[1 0 1],'Normalization','probability')
legend('d16','p95')
xlabel('Bootstrap estimates of p_0')
ylabel('Frequency')
title('FOSL c = 9')
%exportgraphics(f,'figures\FOSL_p0bootstrap_c90_Fols.pdf','Resolution',600)


% p1
f = figure;
hold on
histogram(plot_params_d16(26:975,2),20,'FaceColor',[1 1 0],'Normalization','probability')
histogram(plot_params_p95(26:975,2),20,'FaceColor',[1 0 1],'Normalization','probability')
legend('d16','p95')
xlabel('Bootstrap estimates of p_1')
ylabel('Frequency')
title('FOSL c = 9')
%exportgraphics(f,'figures\FOSL_p1bootstrap_c90_Fols.pdf','Resolution',600)


% p2
f = figure;
hold on
histogram(plot_params_d16(1:995,3),20,'FaceColor',[1 1 0],'Normalization','probability')
histogram(plot_params_p95(:,3),20,'FaceColor',[1 0 1],'Normalization','probability')
legend('d16','p95')
xlabel('Bootstrap estimates of p_2')
ylabel('Frequency')
title('FOSL c = 9')
%exportgraphics(f,'figures\FOSL_p2bootstrap_c90_Fols.pdf','Resolution',600)


% % % % % % % %
%% Histograms FOSL S / L cuttoff = 1.5
% S1 B,C,D
plot_params_d16 = sort(mouse_model_d16_3C_15cutoff.param_estimate(:,4:6));
plot_params_p95 = sort(mouse_model_p95_3C_15cutoff.param_estimate(:,4:6));


% p0
f = figure;
hold on
histogram(plot_params_d16(26:975,1),20,'FaceColor',[1 1 0],'Normalization','probability')
histogram(plot_params_p95(26:975,1),20,'FaceColor',[1 0 1],'Normalization','probability')
legend('d16','p95')
xlabel('Bootstrap estimates of p_0')
ylabel('Frequency')
title('FOSL c = 1.5')
%exportgraphics(f,'figures\FOSL_p0bootstrap_c15_Fols.pdf','Resolution',600)


% p1
f = figure;
hold on
histogram(plot_params_d16(26:975,2),20,'FaceColor',[1 1 0],'Normalization','probability')
histogram(plot_params_p95(26:975,2),20,'FaceColor',[1 0 1],'Normalization','probability')
legend('d16','p95')
xlabel('Bootstrap estimates of p_1')
ylabel('Frequency')
title('FSL c = 1.5')
%exportgraphics(f,'figures\FOSL_p1bootstrap_c15_Fols.pdf','Resolution',600)

% p2
f = figure;
hold on
histogram(plot_params_d16(:,3),20,'FaceColor',[1 1 0],'Normalization','probability')
histogram(plot_params_p95(:,3),20,'FaceColor',[1 0 1],'Normalization','probability')
legend('d16','p95')
xlabel('Bootstrap estimates of p_2')
ylabel('Frequency')
title('FSL c = 1.5')
%exportgraphics(f,'figures\FOSL_p2bootstrap_c15_Fols.pdf','Resolution',600)

     


