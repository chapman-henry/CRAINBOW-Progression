% FSL expected tumors
% July 2024
% % FIGURE 5C/D
% note --- data save('FSL_expected_tumor_sims.mat','n_large_d','n_large_p','N_k','t_b')
% contains data used to make FIG 5C/D

load('bootstrap_NoOccult_JAN2024.mat')
%%
mean_d16table = mean(param_results.d16_params);
mean_p95table = mean(param_results.p95_params);

% parameters for FSL model
Pd.p1 = mean_d16table.P1;
Pd.p2 = mean_d16table.P2;
Pd.r = mean_d16table.r;
Pd.K = mean_d16table.K;
Pd.T0 = mean_d16table.T0;

Pp.p1 = mean_p95table.P1;
Pp.p2 = mean_p95table.P2;
Pp.r = mean_p95table.r;
Pp.K = mean_p95table.K;
Pp.T0 = mean_p95table.T0;

%% Loop over Carrying capacity (k) and observation time (t_s)

%Log2 carrying capacity values
N_k = 8:0.2:18;
% time grid for solving ODE
t_b = 0:0.5:100;

n_T = length(t_b);
n_K = length(N_k);

% allocating memory for output: number of large tumors at end point
large_d = zeros(n_K,n_T);
large_p = zeros(n_K,n_T);

tic;
for k = 1:n_K
    for t_s = 1:n_T
        
        %setting time, starting at 
        t = t_b(t_s):0.1:104;

        % setting number of observed transformed cells
        % .T0 is the observation
        % and we suppose the number of transformed cells stays constant
        % over life of the sim mouse
        % same values for p95 and d16
        Pp.K = 2^(N_k(k));
        Pp.T0 = Pp.K;
        Pd.K = 2^(N_k(k));
        Pd.T0 = Pd.K;
        
        %solve eqns
        Sd = Kinetics_FSL(t,Pd);
        Sp = Kinetics_FSL(t,Pp);

        large_d(k,t_s) = Sd(end,3);
        large_p(k,t_s) = Sp(end,3);
        
    end
    %disp(k)
end
toc

%save('FSL_expected_tumor_sims.mat','n_large_d','n_large_p','N_k','t_b')
%% adjusting variables for plots
K_plot = N_k;
for k = 1:n_K
    K_plot(k) = 2^K_plot(k);
end

large_d_plot = mean(large_d,3);
large_d_plot = log10(flipud(large_d_plot));

large_p_plot = mean(large_p,3);
large_p_plot = log10(flipud(large_p_plot));

%%
% d16
figure;
%imagesc(t_b,K_plot,large_d_plot)
imagesc(large_d_plot)
colorbar;
xlabel('Age of mouse at detection(weeks)')
ylabel('Number of transformed cells (log 10)')
title('Expected log large d16 tumors at 104 weeks')
%%
% p95
figure;
%imagesc(t_b,K_plot,large_p_plot)
imagesc(large_p_plot)
colorbar;
xlabel('Age of mouse at detection(weeks)')
ylabel('Number of transformed cells (log 10)')
title('Expected log large p95 tumors at 104 weeks')