% SensitivityKchannel.m
%
% Sensitivity analysis: the potassium-channel gain functions.
%
% The shear-gated K+ conductance and the resulting membrane potential are:
%     gK(tau) = ksens*|tau| + gK0                         [S/m^2]
%     Em(tau) = (gNa*ENa + gLeak*ELeak + gK*EK) / gm       [V],  gm = gNa+gLeak+gK
% Em is the drive for adaptation (rdot ~ Vref - Em). ksens is the gain: it sets
% how strongly shear hyperpolarises the vessel.
%
% Part A: plot the gain functions gK(tau) and Em(tau) for several ksens (formula
%         only, no simulation).
% Part B: sweep ksens and measure loop survival (ArtCoupling, FlowReg).
%
% Output: FinalResults/SensitivityKchannel.png + _Results.mat

rng(1234);
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
close all;
PROJ = fileparts(mfilename('fullpath'));
addpath(PROJ); addpath(fullfile(PROJ,'ElectricModel'));

% K+ model constants (defaults from defaultPars.m)
gNa=0.2; gLeak=0.2; gK0=0.5; ENa=60e-3; ELeak=0; EK=-90e-3; Vref=-0.05;

%% Part A: gain functions (formula only)
tau = linspace(0, 3, 200);                 % wall shear stress [Pa]
ksens_curves = [0.1 0.4 1.0];              % gain values to illustrate
cmap = lines(numel(ksens_curves));

figure('Color','w','Position',[120 200 1000 640]);

subplot(2,2,1); hold on;
for k = 1:numel(ksens_curves)
    gK = ksens_curves(k)*abs(tau) + gK0;
    plot(tau, gK, 'LineWidth',1.8, 'Color',cmap(k,:), 'DisplayName',sprintf('ksens=%.1f',ksens_curves(k)));
end
xlabel('wall shear stress \tau [Pa]'); ylabel('g_K [S/m^2]');
title('K^+ conductance gain  g_K(\tau) = ksens\cdot|\tau| + g_{K0}');
legend('Location','northwest'); grid on; box off;

subplot(2,2,2); hold on;
for k = 1:numel(ksens_curves)
    gK = ksens_curves(k)*abs(tau) + gK0;
    gm = gNa + gLeak + gK;
    Em = (gNa*ENa + gLeak*ELeak + gK.*EK)./gm;
    plot(tau, Em*1e3, 'LineWidth',1.8, 'Color',cmap(k,:), 'DisplayName',sprintf('ksens=%.1f',ksens_curves(k)));
end
yline(Vref*1e3,'k--','V_{ref}','HandleVisibility','off');
xlabel('wall shear stress \tau [Pa]'); ylabel('E_m [mV]');
title('membrane potential  E_m(\tau)   (drive = V_{ref} - E_m)');
legend('Location','northeast'); grid on; box off;

%% Part B: ksens sweep -> loops
r_dead = 5e-6;  P_op = 86;  MAX_CHUNKS = 150;  density = 3.35;  ECr = 0.1;
fprintf('Building network...\n');
[S_ref, pial_idx, ~, Qreg, nloops_total, A_mm2, npen] = sens_build(2, 2e-3, density, 1234);
fprintf('  area=%.1f mm^2  n_pen=%d  n_loops_total=%d\n\n', A_mm2, npen, nloops_total);

ksens_list = [0.6 0.7 0.8 0.9 1.0 1.1 1.2];   % focused 0.6-1.2 (FlowReg basin transition)
models = {'ArtCoupling','FullCouplingFlowReg'};
nM = numel(models);
loops = nan(numel(ksens_list), nM);

for ik = 1:numel(ksens_list)
    S_use = S_ref;  S_use.ksens = ksens_list(ik);       % override the K+ gain
    fprintf('ksens=%.2f\n', ksens_list(ik));
    for im = 1:nM
        t0 = tic;
        [nl,~,cv] = sens_ss(S_use, models{im}, Qreg, ECr, P_op*133, pial_idx, r_dead, MAX_CHUNKS, []);
        loops(ik,im) = nl;
        fprintf('   %-22s loops=%2d/%2d  [%.0fs conv=%d]\n', models{im}, nl, nloops_total, toc(t0), cv);
    end
end

% table
fprintf('\n%s\nLOOPS vs ksens (ECresistivity=%.2g)\n%s\n', repmat('=',1,52), ECr, repmat('=',1,52));
fprintf('%-8s', 'ksens'); for im=1:nM, fprintf(' %22s', models{im}); end; fprintf('\n');
for ik=1:numel(ksens_list)
    fprintf('%-8.2f', ksens_list(ik));
    for im=1:nM, fprintf(' %16d/%-5d', loops(ik,im), nloops_total); end
    fprintf('\n');
end

subplot(2,2,[3 4]); hold on;
mc = [0.85 0.45 0.05; 0.05 0.60 0.20];
for im=1:nM
    plot(ksens_list, loops(:,im), '-o','LineWidth',1.9,'Color',mc(im,:), ...
        'MarkerFaceColor',mc(im,:),'DisplayName',models{im});
end
xline(0.4,'k:','default','HandleVisibility','off');
yline(nloops_total,'k:','max','HandleVisibility','off');
xlabel('K^+ shear sensitivity  ksens'); ylabel('surviving loops');
title('Loop survival vs K^+ channel gain (ksens)'); grid on; box off; legend('Location','southeast');

sgtitle('Potassium-channel gain: functions and loop sensitivity','FontWeight','bold');
saveas(gcf, fullfile(PROJ,'FinalResults','SensitivityKchannel.png'));
save(fullfile(PROJ,'FinalResults','SensitivityKchannel_Results.mat'), ...
    'ksens_list','loops','models','nloops_total','ksens_curves','ECr', ...
    'gNa','gLeak','gK0','ENa','ELeak','EK','Vref');
fprintf('\nSaved: FinalResults/SensitivityKchannel.png and _Results.mat\n');
