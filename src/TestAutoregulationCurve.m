% TestAutoregulationCurve.m
%
% Generates the autoregulation curve:
%
%   Y-axis: Q_total / Q_0   (normalised to each model's own operating flow)
%   X-axis: arterial pressure P (mmHg), range 50-120 mmHg
%
% All adapted curves pass through (P_0=86, 1) by construction.
% Passive reference: straight line Q_norm = P/P_0 through origin
%   (exact for fixed-conductance network -- no simulation needed).
% Perfect autoregulation: horizontal line at Q_norm = 1.
%
% Autoregulation index (Lassen 1959):
%   AI = 1 - (dQ_norm/dP) / (1/P_0)  =  1 - P_0 * (dQ_norm/dP)
% where dQ_norm/dP is the local slope of the normalised adapted curve at P_0,
% estimated from the adjacent pressure points (80 and 90 mmHg).
%
% Why curves are NOT straight lines:
%   Vessels re-adapt from r_ss_86 at each new pressure. At low P some
%   pial segments collapse (r < r_dead), sharply reducing network
%   conductance and bending the curve downward. At high P over-dilation
%   causes the curve to bend upward. This nonlinearity is only visible
%   with a wide pressure range and enough sample points.

rng(1234);
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
close all;

PROJ = fileparts(mfilename('fullpath')); addpath(PROJ);
addpath(fullfile(PROJ,'ElectricModel'));

r_dead = 5e-6;
P_op   = 86;    % mmHg, operating pressure
KSENSQ = 3.0;   % oxygen-sensing gain (peak of the loop-stability inverted-U; only affects FlowReg)

%% Network geometry (identical to TestPhysiologicalPenetrators)
ElementLength = 2e-3;
ncombsx       = 2;
ncombsy       = 4;
n_hex         = ncombsx * ncombsy;
A_hex         = (3*sqrt(3)/2) * ElementLength^2;
A_total_mm2   = n_hex * A_hex * 1e6;
pen_density   = 3.35;    % penetrators/mm^2 (measured: 23 / 6.86 mm^2 crop)
n_pen_target  = round(pen_density * A_total_mm2);

fprintf('Building network  (%.1f mm^2, %d penetrators)...\n', A_total_mm2, n_pen_target);
Sin = struct;
Sin.npenetrator         = 1;
Sin.ncombsx             = ncombsx;
Sin.branchlevel         = 4;
Sin.balance_Gsin_Gsout  = 10;
Sin.Vref                = -0.05;
Sin.gainRadius2Vref     = 10;
Sin.ksens               = 0.4;
Sin.kReg                = 0.01;
Sin.CapillaryPotential  = -0.06;
Sin.CapillaryConductance = 2e-7;
Sin.kSensQ              = 0.4;
Sin.gCapBase            = 0.2;
Sin.ECapBase            = -0.04;
Sin.DoCapillaryCurrent  = 0;
Sin.ECresistivity       = 0.1;
Sin.KeepDynamics        = 0;
Sin.GsSink              = 1.6e-16;

S_ref = defaultPars('leptomeningeal_geometry', Sin);

pen   = load_penetrator_data('penetrators.csv');
sp    = pp_sample_with_replacement(pen, n_pen_target);
S_ref = integrate_sampled_penetrators(S_ref, sp);
n_new = S_ref.nSE - length(S_ref.sources);
S_ref.sources = [S_ref.sources, zeros(1, n_new)];

pial_idx = find([S_ref.IE.pial] == 1);
sink_idx = find(~S_ref.sources);
fprintf('  n_pial=%d (post-split), n_penetrators=%d\n\n', length(pial_idx), length(sink_idx));

%% Pressure sweep points
pressures_mmHg = [50, 60, 70, 80, 86, 90, 100, 110, 120];
nP      = length(pressures_mmHg);
P_op_idx = find(pressures_mmHg == P_op);   % index of operating point

%% Calibrate Qref from NoCoupling SS at 86 mmHg
fprintf('Calibrating Qreg from NoCoupling SS...\n');
S86_cal = ac_setP(S_ref, P_op*133);
[R_cal, S_cal, ~] = ac_run(S86_cal, 0, 0.1, 'NoCoupling', []);
r_cal = R_cal.X(end,:);
S_tmp = S_cal; [S_tmp.IE.r] = vout(r_cal);
S_tmp = calcConductance2025(S_tmp); S_tmp = solvehemodyn2025(S_tmp);
Qreg = mean(abs([S_tmp.SE(sink_idx).Qs]));   % regulation target = mean penetrator flow of NoCoupling SS
fprintf('  Qreg (target) = %.3f fL/s\n\n', Qreg*1e15);

%% Models
models  = {'NoCoupling','ArtCoupling','FullCoupling','FullCouplingFlowReg'};
labels  = {'No coupling','Art. coupling','Full coupling','Full + flow reg.'};
colors  = [0.15 0.35 0.75;   % blue
           0.85 0.45 0.05;   % orange
           0.55 0.10 0.55;   % purple
           0.05 0.60 0.20];  % green
nM = length(models);

Q_norm_all = nan(nM, nP);   % normalised adapted flow [model × pressure]
AI_all     = nan(1, nM);
n_surv_all = nan(nM, nP);   % surviving pial segments at each pressure point

%% Main loop
for imod = 1:nM
    model_name = models{imod};
    fprintf('%s\n', repmat('=',1,65));
    fprintf('Model: %s\n', model_name);
    fprintf('%s\n', repmat('=',1,65));

    % Step 1: SS at operating pressure (fresh start)
    S86m = ac_setP(S_ref, P_op*133);
    S86m.QrefCapillary = Qreg;
    S86m.gCapBase      = 0.05;
    S86m.kSensQ        = KSENSQ;
    S86m.GscCapillary  = 2e-5;
    S86m.kSmCapillary  = 7e-6;
    [R86, S86_out, conv86] = ac_run(S86m, Qreg, 0.1, model_name, []);
    r_ss86 = R86.X(end,:);
    fprintf('  86 mmHg SS: %d chunks, converged=%d\n', conv86.n_chunks, conv86.converged);

    % Q_0 at operating point
    S_tmp = S86_out; [S_tmp.IE.r] = vout(r_ss86);
    S_tmp = calcConductance2025(S_tmp); S_tmp = solvehemodyn2025(S_tmp);
    Q0 = mean(abs([S_tmp.SE(sink_idx).Qs]));
    fprintf('  Q_0 = %.3f fL/s,  n_surv_86 = %d/%d\n\n', ...
        Q0*1e15, sum(r_ss86(pial_idx) > r_dead), length(pial_idx));

    % Step 2: Pressure sweep re-adapting from r_ss86
    for iP = 1:nP
        P_mmHg = pressures_mmHg(iP);
        Spm = ac_setP(S_ref, P_mmHg*133);
        Spm.QrefCapillary = Qreg;
        Spm.gCapBase      = 0.05;
        Spm.kSensQ        = KSENSQ;
        Spm.GscCapillary  = 2e-5;
        Spm.kSmCapillary  = 7e-6;
        [Rpm, Spm_out, conv_pm] = ac_run(Spm, Qreg, 0.1, model_name, r_ss86);
        r_fin = Rpm.X(end,:);

        S_tmp2 = Spm_out; [S_tmp2.IE.r] = vout(r_fin);
        S_tmp2 = calcConductance2025(S_tmp2); S_tmp2 = solvehemodyn2025(S_tmp2);
        Q_ad = mean(abs([S_tmp2.SE(sink_idx).Qs]));

        Q_norm_all(imod, iP) = Q_ad / Q0;
        n_surv_all(imod, iP) = sum(r_fin(pial_idx) > r_dead);
        fprintf('  P=%3d mmHg: Q_norm=%.4f  n_surv=%d  chunks=%d  conv=%d\n', ...
            P_mmHg, Q_ad/Q0, n_surv_all(imod,iP), conv_pm.n_chunks, conv_pm.converged);
    end

    % AI from local slope at P_op using adjacent points (80 & 90 mmHg)
    % Indices for 80 and 90 mmHg
    iLo = find(pressures_mmHg == 80);
    iHi = find(pressures_mmHg == 90);
    if ~isempty(iLo) && ~isempty(iHi)
        dQn_dP = (Q_norm_all(imod,iHi) - Q_norm_all(imod,iLo)) / ...
                 (pressures_mmHg(iHi) - pressures_mmHg(iLo));
    else
        % Fallback: use full range slope
        dQn_dP = (Q_norm_all(imod,end) - Q_norm_all(imod,1)) / ...
                 (pressures_mmHg(end) - pressures_mmHg(1));
    end
    % Lassen AI = 1 - (dQ/dP) / (Q_0/P_0)
    % In normalised coords: Q_norm = Q/Q_0 -> dQn/dP = dQ/(Q_0 dP)
    % Reference slope = (Q_0/P_0)/Q_0 = 1/P_0
    AI_all(imod) = 1 - dQn_dP * P_op;
    fprintf('\n  AI = %.3f  (local slope at 86: dQ_norm/dP = %.5f /mmHg)\n\n', ...
        AI_all(imod), dQn_dP);
end

%% Figure
figure('Position', [100 80 750 520], 'Color', 'w');
hold on;

% Passive reference: Q_norm = P / P_op  (exact for fixed conductance)
P_line = linspace(0, 130, 300);
h_pass = plot(P_line, P_line/P_op, 'k--', 'LineWidth', 1.5);

% Perfect autoregulation: horizontal at 1
h_perf = yline(1, 'k:', 'LineWidth', 1.2);

% Model curves
h_mod = gobjects(nM,1);
for imod = 1:nM
    h_mod(imod) = plot(pressures_mmHg, Q_norm_all(imod,:), '-o', ...
        'Color', colors(imod,:), 'LineWidth', 2.2, 'MarkerSize', 6, ...
        'MarkerFaceColor', colors(imod,:));
end

% Operating point marker
plot(P_op, 1, 'ko', 'MarkerSize', 11, 'MarkerFaceColor', [0.9 0.9 0.9], ...
    'LineWidth', 1.5, 'HandleVisibility', 'off');
text(P_op + 1.5, 1.03, sprintf('P_0 = %d mmHg', P_op), 'FontSize', 9.5);

% Legend
legend_entries = [h_pass; h_mod(:)];
legend_labels  = ['Passive (no adapt.)'; ...
    cellfun(@(lbl, ai) sprintf('%s  (AI = %.3f)', lbl, ai), ...
            labels', num2cell(AI_all)', 'UniformOutput', false)];
legend(legend_entries, legend_labels, 'Location', 'northwest', 'FontSize', 10);

% Axes
xlabel('Arterial pressure  (mmHg)', 'FontSize', 13);
ylabel('Q_{total} / Q_0   (normalised flow)', 'FontSize', 13);
title('Autoregulation curve — effect of endothelial coupling', 'FontSize', 13, 'FontWeight', 'bold');
xlim([35 130]);
ylim([0 2.0]);
grid on; box on;
set(gca, 'FontSize', 11, 'LineWidth', 1);

% Interpretation note
annotation('textbox', [0.60 0.12 0.28 0.14], ...
    'String', {'AI = 0 : no regulation', 'AI = 1 : perfect regulation', 'AI < 0 : WSS amplification'}, ...
    'FitBoxToText', 'on', 'BackgroundColor', [0.97 0.97 0.97], ...
    'EdgeColor', [0.7 0.7 0.7], 'FontSize', 9);

saveas(gcf, fullfile(PROJ,'FinalResults','AutoregulationCurve.png'));
fprintf('Figure saved to AutoregulationCurve.png\n');

%% Summary
fprintf('\n%s\n', repmat('=',1,70));
fprintf('SUMMARY\n');
fprintf('%s\n', repmat('=',1,70));
fprintf('%-25s %8s %8s %8s\n', 'Model', 'AI', 'n_surv@86', 'n_surv@50');
for imod = 1:nM
    fprintf('%-25s %8.3f %8d %8d\n', models{imod}, AI_all(imod), ...
        n_surv_all(imod, P_op_idx), n_surv_all(imod, 1));
end
fprintf('\nPassive reference: Q_norm = P / %d  (analytical, no simulation)\n', P_op);

save(fullfile(PROJ,'FinalResults','AutoregulationCurve_Results.mat'), 'Q_norm_all', 'n_surv_all', ...
    'pressures_mmHg', 'AI_all', 'models', 'labels', 'P_op', 'Qreg', 'KSENSQ');
fprintf('Results saved to AutoregulationCurve_Results.mat\n');


%% ════════════════════════════════════════════════════════════════════════
%% LOCAL HELPERS
%% ════════════════════════════════════════════════════════════════════════

function sp = pp_sample_with_replacement(pen, n)
    n_avail = length(pen.r);
    idx = randi(n_avail, n, 1);
    sp.r = pen.r(idx);
    sp.l = pen.l(idx);
    sp.R = pen.R(idx);
    sp.G = 1 ./ pen.R(idx);
end

function [R, S_out, conv_info] = ac_run(S_geo, QrefCap, ECr, model_name, r0_init)
    chunk_size = 1e5;
    tol_drdt   = 1e-12;
    max_chunks = 150;

    S = S_geo;
    S.ECresistivity = ECr;
    S.QrefCapillary = QrefCap;

    pen_ie = find(~[S.IE.pial]);
    S.kReg = S.kReg * ones(S.nIE, 1);
    S.kReg(pen_ie) = 0;

    if nargin >= 5 && ~isempty(r0_init)
        S.r0 = r0_init(:)';
    else
        S.r0 = 100e-6 * ones(1, S.nIE);
        S.r0(pen_ie) = 15e-6;
    end
    S.KeepDynamics = 0;
    needs_caps = any(strcmp(model_name, {'FullCoupling','FullCouplingFlowReg'}));
    S.DoCapillaryCurrent = needs_caps;
    S = MakeElectricCircuit(S, 0);
    S.DoElectricCoupling       = ~strcmp(model_name, 'NoCoupling');
    S.DoCapillaryCurrent       = needs_caps;
    S.FlowRegulatedCapillaries = strcmp(model_name, 'FullCouplingFlowReg');
    S.annot = model_name;

    opts  = odeset('Events', @ac_detectNaN, 'RelTol', 1e-6, 'AbsTol', 1e-9);
    t_now = 0;  y_now = S.r0(:);  converged = false;  n_chunks = 0;  drdt_last = NaN;

    while n_chunks < max_chunks
        [t_c, X_c, tErr] = ode45(@(t,X) rdotfunCoupling(t, X, S), ...
                                  [t_now, t_now + chunk_size], y_now, opts);
        n_chunks = n_chunks + 1;
        if ~isempty(tErr), y_now = X_c(end,:)';  break, end
        drdt_last = max(abs(X_c(end,:) - X_c(1,:))) / (t_c(end) - t_c(1));
        if drdt_last < tol_drdt, converged = true;  y_now = X_c(end,:)';  break, end
        t_now = t_c(end);  y_now = X_c(end,:)';
    end

    R.X = y_now';
    R.t = t_now + chunk_size * converged;
    conv_info.converged  = converged;
    conv_info.n_chunks   = n_chunks;
    conv_info.drdt_final = drdt_last;
    conv_info.t_total    = R.t;
    S_out = S;
end

function [value, isterminal, direction] = ac_detectNaN(~, X)
    value = all(~isnan(X));  isterminal = 1;  direction = 0;
end

function S_out = ac_setP(S, Psrc)
    S_out = S;
    S_out.sourceP = Psrc;
    for k = find(S_out.sources)
        S_out.SE(k).Ps = Psrc;
    end
end
