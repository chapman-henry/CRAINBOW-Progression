% Morris screen sensitivity analysis for the CRAINBOW FOSL model
% The inputs are as follows:
%   r: scalar denoting the number of samples to take of each parameter
%   p: parameter structure containing one value for each of the following
%   parameters: p0, p1, p2, F0, r, and k
%   genotype: string denoting the genotype being examined (used in making
%   graph titles)
%   output_name: string containing the output file name to save the EET
%   results and parameter samples. If a specific location is desired,
%   provide the full path. Otherwise, the file will be saved to the current
%   working directory
function CRAINBOW_field_SA_FOSL(r, p, genotype, output_name)
%addpath 'C:\Users\htc9\Documents\MATLAB\SAFE-matlab-main\sampling';
%addpath 'C:\Users\htc9\Documents\MATLAB\SAFE-matlab-main\EET';
%addpath 'C:\Users\htc9\Documents\MATLAB\SAFE-matlab-main\VBSA';
addpath '/hpc/home/htc9/EET';
addpath '/hpc/home/htc9/sampling';

%% parameters
t_span = 3:.1:25;
params = horzcat(p.p0, p.p1, p.p2, p.r, p.k, p.F0);
pars_sup = .9*params;
pars_inf = 1.1*params;
X_labels = {'P_0', 'P_1', 'P_2', 'r', 'k', 'F0'};
empStruct = struct('P_0', {}, 'P_1', {}, 'P_2', {}, 'r', {}, 'k', {}, 'F0', {});
A = struct('F_mu', empStruct, 'F_sig', empStruct, ...
    'H_mu', empStruct, 'H_sig', empStruct, ...
    'S_mu', empStruct, 'S_sig', empStruct, ...
    'L_mu', empStruct, 'L_sig', empStruct);
field_names = fieldnames(empStruct);
L = length(X_labels);
DistrPar = cell(1,L);
for i=1:L; DistrPar{i} = [pars_sup(i) pars_inf(i)]; end

DistrFun = cell(1,L);
for i=1:L; DistrFun{i} = 'unif'; end

SampStrategy = 'lhs';
design_type = 'radial';

X = OAT_sampling(r, L, DistrFun, DistrPar, SampStrategy, design_type);



% max_viral_load = zeros(M,height(epsilon_samples));
%figure

% myCluster = parcluster('local');
% myCluster.NumWorkers = 24;  % 'Modified' property now TRUE
% saveProfile(myCluster);
% parpool('local',16);


[F_traj, H_traj, S_traj, L_traj] = deal(zeros(length(X), length(t_span)));
tic
%fprintf('\nCounter: ')
parfor ii = 1:length(X)
    trajs = KineticsX_wrap(X(ii,:), t_span);
    [F_traj(ii,:), H_traj(ii,:), S_traj(ii,:), L_traj(ii,:)] = ...
        deal(trajs.F, trajs.H, trajs.S, trajs.L);
end
T = width(F_traj);
[F_mu, F_sig, H_mu, H_sig, S_mu, S_sig, L_mu, L_sig] = ...
        deal(zeros(T,L));

for jj = 1:T
    [F_mu(jj,:), F_sig(jj,:)] = EET_indices(r, pars_sup, pars_inf, X, F_traj(:,jj), design_type);
    [H_mu(jj,:), H_sig(jj,:)] = EET_indices(r, pars_sup, pars_inf, X, H_traj(:,jj), design_type);
    [S_mu(jj,:), S_sig(jj,:)] = EET_indices(r, pars_sup, pars_inf, X, S_traj(:,jj), design_type);
    [L_mu(jj,:), L_sig(jj,:)] = EET_indices(r, pars_sup, pars_inf, X, L_traj(:,jj), design_type);
end
for k = 1:length(field_names)
    A.F_mu(1).(field_names{k}) = F_mu(:,k); A.F_sig(1).(field_names{k}) = F_sig(:,k);
    A.H_mu(1).(field_names{k}) = H_mu(:,k); A.H_sig(1).(field_names{k}) = H_sig(:,k);
    A.S_mu(1).(field_names{k}) = S_mu(:,k); A.S_sig(1).(field_names{k}) = S_sig(:,k);
    A.L_mu(1).(field_names{k}) = L_mu(:,k); A.L_sig(1).(field_names{k}) = L_sig(:,k);
end
    
parameter_samples = X;
EET_results = A;
save(output_name, "parameter_samples", "EET_results");

toc
figure()
morris_mu = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
set(gcf, 'DefaultAxesColorOrder', [turbo(8)])
nexttile(1)
plot(t_span, F_mu, 'LineWidth', 1.5)
grid on
title('Field', 'FontSize', 16)
ylabel('\mu^*', 'FontSize', 16)
xlabel('Time (weeks)', 'FontSize', 16)


nexttile(2)
plot(t_span, H_mu, 'LineWidth', 1.5)
grid on
title('Occult', 'FontSize', 16)
ylabel('\mu^*', 'FontSize', 16)
xlabel('Time (weeks)', 'FontSize', 16)

nexttile(3)
plot(t_span, S_mu, 'LineWidth', 1.5)
grid on
title('Small', 'FontSize', 16)
ylabel('\mu^*', 'FontSize', 16)
xlabel('Time (weeks)', 'FontSize', 16)

nexttile(4)
plot(t_span, L_mu, 'LineWidth', 1.5)
grid on
title('Large', 'FontSize', 16)
ylabel('\mu^*', 'FontSize', 16)
xlabel('Time (weeks)', 'FontSize', 16)

hobj = legend(X_labels, 'Orientation', 'vertical', 'FontSize', 14);
%     h1 = findobj(hobj, 'type', 'line');
hobj.Layout.Tile = 'east';
%     set(h1, 'LineWidth', 3);

title(morris_mu, {genotype, ' Tumor Counts, Morris Method, \mu^*'}, 'FontSize', 20)
fig_name = [output_name, '_fig.pdf'];
exportgraphics(morris_mu, fig_name, 'Resolution',600)
% For the file destination change to match desired file. File path is not
% automatically created
exportgraphics(morris_mu, 'CRAINBOW_SA_results_2/all_figs.pdf', 'Resolution', 600, 'Append', true)
% 
% figure()
% morris_sig = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
% 
% nexttile(1)
% plot(t_span, F_sig, 'LineWidth', 1.5)
% grid on
% title('Field', 'FontSize', 16)
% ylabel('\sigma^*', 'FontSize', 16)
% xlabel('Time (weeks)', 'FontSize', 16)
% set(gcf, 'DefaultAxesColorOrder', [turbo(7)])
% 
% nexttile(2)
% plot(t_span, H_sig, 'LineWidth', 1.5)
% grid on
% title('Occult', 'FontSize', 16)
% ylabel('\sigma^*', 'FontSize', 16)
% xlabel('Time (weeks)', 'FontSize', 16)
% set(gcf, 'DefaultAxesColorOrder', [turbo(7)])
% 
% nexttile(3)
% plot(t_span, S_sig, 'LineWidth', 1.5)
% grid on
% title('Small', 'FontSize', 16)
% ylabel('\sigma^*', 'FontSize', 16)
% xlabel('Time (weeks)', 'FontSize', 16)
% set(gcf, 'DefaultAxesColorOrder', [turbo(7)])
% 
% nexttile(4)
% plot(t_span, L_sig, 'LineWidth', 1.5)
% grid on
% title('Large', 'FontSize', 16)
% ylabel('\sigma^*', 'FontSize', 16)
% xlabel('Time (weeks)', 'FontSize', 16)
% set(gcf, 'DefaultAxesColorOrder', [turbo(7)])
% 
% hobj = legend(X_labels, 'Orientation', 'vertical', 'FontSize', 14);
% %     h1 = findobj(hobj, 'type', 'line');
% hobj.Layout.Tile = 'east';
% %     set(h1, 'LineWidth', 3);
% 
% title(morris_sig, [genotype, ' Tumor Counts, Morris Method, \sigma^*'], 'FontSize', 20)
end
%%

function trajs = KineticsX_wrap(x, t_span)
    par_cell = num2cell(x);

    [P.p0, P.p1, P.p2, P.r, P.k, P.F0] = deal(par_cell{:});
    y0 = [P.F0, 0, 0, 0];
    [~,trajs_out] = ode45(@(t,y) growth(t,y,P),t_span,y0);


    [trajs.F, trajs.H, trajs.S, trajs.L] = ...
        deal(trajs_out(:,1)', trajs_out(:,2)', trajs_out(:,3)', trajs_out(:,4)');

    
	function dydt = growth(t,y,P)
		
		dydt = zeros(4, 1);
        F = y(1);
        H = y(2);
        S = y(3);
        L = y(4);
		dydt(1) = P.r*F*(1-F/(P.k))-P.p0*F; %F
		dydt(2) = P.p0*F - P.p1*H; %Unobservable
		dydt(3) = P.p1*H - P.p2*S; %Small
		dydt(4) = P.p2*S;          %Large
    end

end


%%