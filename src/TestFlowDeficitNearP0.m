% TestFlowDeficitNearP0.m
%
% Fine pressure sweep below the operating point (75-86 mmHg, 1 mmHg steps).
% Runs ArtCoupling and FullCouplingFlowReg only. At each pressure it records
% the normalised flow Q/Q_0, the surviving pial segment count, and:
%     adapted deficit   D_ad   = 1 - Q_norm
%     passive deficit   D_pass = 1 - P/P_0
%     compensation      C      = Q_norm - P/P_0      (>0 above passive)
% plus a local autoregulation index on each 1 mmHg bracket:
%     AI_loc = 1 - P_0 * dQ_norm/dP
% This resolves where the flow-regulated model sits relative to the passive
% reference across the range.

rng(1234);
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
close all;

PROJ = fileparts(mfilename('fullpath')); addpath(PROJ);
addpath(fullfile(PROJ,'ElectricModel'));

r_dead = 5e-6;
P_op   = 86;    % mmHg, operating pressure
KSENSQ = 3.0;   % oxygen-sensing gain (peak of the loop-stability inverted-U; only affects FlowReg)

%% Network geometry (identical to TestAutoregulationCurve)
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

%% Fine pressure sweep points (below and up to the set point)
pressures_mmHg = 75:1:86;      % 1 mmHg resolution, 75 -> 86
nP       = length(pressures_mmHg);
P_op_idx = find(pressures_mmHg == P_op);

%% Calibrate Qref from NoCoupling SS at 86 mmHg
fprintf('Calibrating Qreg from NoCoupling SS...\n');
S86_cal = ac_setP(S_ref, P_op*133);
[R_cal, S_cal, ~] = ac_run(S86_cal, 0, 0.1, 'NoCoupling', []);
r_cal = R_cal.X(end,:);
S_tmp = S_cal; [S_tmp.IE.r] = vout(r_cal);
S_tmp = calcConductance2025(S_tmp); S_tmp = solvehemodyn2025(S_tmp);
Qreg = mean(abs([S_tmp.SE(sink_idx).Qs]));
fprintf('  Qreg (target) = %.3f fL/s\n\n', Qreg*1e15);

%% Models (only the two informative ones)
models  = {'ArtCoupling','FullCouplingFlowReg'};
labels  = {'Art. coupling','Full + flow reg.'};
colors  = [0.85 0.45 0.05;   % orange
           0.05 0.60 0.20];  % green
nM = length(models);

Q_norm_all = nan(nM, nP);
n_surv_all = nan(nM, nP);
conv_all   = false(nM, nP);

%% Main loop
for imod = 1:nM
    model_name = models{imod};
    fprintf('%s\n', repmat('=',1,65));
    fprintf('Model: %s\n', model_name);
    fprintf('%s\n', repmat('=',1,65));

    % SS at operating pressure (fresh start)
    S86m = ac_setP(S_ref, P_op*133);
    S86m.QrefCapillary = Qreg;
    S86m.gCapBase      = 0.05;
    S86m.kSensQ        = KSENSQ;
    S86m.GscCapillary  = 2e-5;
    S86m.kSmCapillary  = 7e-6;
    [R86, S86_out, conv86] = ac_run(S86m, Qreg, 0.1, model_name, []);
    r_ss86 = R86.X(end,:);

    S_tmp = S86_out; [S_tmp.IE.r] = vout(r_ss86);
    S_tmp = calcConductance2025(S_tmp); S_tmp = solvehemodyn2025(S_tmp);
    Q0 = mean(abs([S_tmp.SE(sink_idx).Qs]));
    fprintf('  86 mmHg SS: Q_0=%.3f fL/s, n_surv_86=%d/%d, conv=%d\n\n', ...
        Q0*1e15, sum(r_ss86(pial_idx) > r_dead), length(pial_idx), conv86.converged);

    % Fine pressure sweep re-adapting from r_ss86
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
        conv_all(imod, iP)   = conv_pm.converged;
        fprintf('  P=%3d mmHg: Q_norm=%.4f  D_ad=%.4f  n_surv=%d  conv=%d\n', ...
            P_mmHg, Q_ad/Q0, 1 - Q_ad/Q0, n_surv_all(imod,iP), conv_pm.converged);
    end
    fprintf('\n');
end

%% Deficit descriptors and local AI
Q_pass    = pressures_mmHg / P_op;              % passive reference (row)
D_pass    = 1 - Q_pass;                          % passive deficit
D_ad      = 1 - Q_norm_all;                      % adapted deficit  [model x P]
Compens   = Q_norm_all - repmat(Q_pass, nM, 1);  % >0 regulating, <0 amplifying

Pmid   = (pressures_mmHg(1:end-1) + pressures_mmHg(2:end)) / 2;
AI_loc = nan(nM, numel(Pmid));
for imod = 1:nM
    dQ = diff(Q_norm_all(imod,:));
    dP = diff(pressures_mmHg);
    AI_loc(imod,:) = 1 - P_op * (dQ ./ dP);
end

%% Figure
figure('Color','w','Position',[100 120 1150 380]);

% (A) normalised flow vs pressure, zoomed
subplot(1,3,1); hold on;
plot(pressures_mmHg, Q_pass, 'k--', 'LineWidth',1.3);           % passive
yline(1,'k:','LineWidth',1.1);                                   % perfect
for imod = 1:nM
    plot(pressures_mmHg, Q_norm_all(imod,:), '-o','LineWidth',2, ...
        'Color',colors(imod,:),'MarkerFaceColor',colors(imod,:),'DisplayName',labels{imod});
end
xline(P_op,'k:','HandleVisibility','off');
xlabel('Arterial pressure (mmHg)'); ylabel('Q / Q_0');
title('Normalised flow (75–86 mmHg)'); grid on; box off;
legend([{'Passive'} labels],'Location','southeast','FontSize',8);

% (B) local AI per 1 mmHg bracket - the regulation notch
subplot(1,3,2); hold on;
for imod = 1:nM
    plot(Pmid, AI_loc(imod,:), '-o','LineWidth',2, ...
        'Color',colors(imod,:),'MarkerFaceColor',colors(imod,:),'DisplayName',labels{imod});
end
yline(0,'k-','LineWidth',1,'HandleVisibility','off');
xline(P_op,'k:','Label','P_0=86','HandleVisibility','off');
xlabel('Arterial pressure (mmHg)'); ylabel('local AI (per bracket)');
title('Regulation strength vs distance from P_0'); grid on; box off;
legend('Location','best','FontSize',8);

% (C) compensation relative to passive
subplot(1,3,3); hold on;
for imod = 1:nM
    plot(pressures_mmHg, Compens(imod,:), '-o','LineWidth',2, ...
        'Color',colors(imod,:),'MarkerFaceColor',colors(imod,:),'DisplayName',labels{imod});
end
yline(0,'k-','LineWidth',1,'HandleVisibility','off');
xline(P_op,'k:','HandleVisibility','off');
xlabel('Arterial pressure (mmHg)');
ylabel('compensation  Q_{norm} - P/P_0');
title('Deficit compensation (>0 regulating)'); grid on; box off;
legend('Location','best','FontSize',8);

sgtitle('Flow-deficit stabilisation near the operating point','FontWeight','bold');
saveas(gcf,fullfile(PROJ,'FinalResults','FlowDeficitNearP0.png'));
fprintf('Saved FlowDeficitNearP0.png\n');

%% Summary tables
fprintf('\n%s\n', repmat('=',1,78));
fprintf('SUMMARY — flow deficit near P_0 = %d mmHg\n', P_op);
fprintf('%s\n', repmat('=',1,78));
for imod = 1:nM
    fprintf('\n%s\n', labels{imod});
    fprintf('  %-6s %8s %8s %8s %10s %8s %6s\n', ...
        'P', 'Q_norm', 'D_ad', 'D_pass', 'compens', 'n_surv', 'conv');
    for iP = 1:nP
        fprintf('  %-6d %8.4f %8.4f %8.4f %10.4f %8d %6d\n', ...
            pressures_mmHg(iP), Q_norm_all(imod,iP), D_ad(imod,iP), ...
            D_pass(iP), Compens(imod,iP), n_surv_all(imod,iP), conv_all(imod,iP));
    end
end

fprintf('\nLocal AI per 1 mmHg bracket:\n');
fprintf('%-22s', 'bracket (mmHg)');
for k = 1:numel(Pmid), fprintf('%8d-%-2d', pressures_mmHg(k), pressures_mmHg(k+1)); end
fprintf('\n');
for imod = 1:nM
    fprintf('%-22s', labels{imod});
    for k = 1:numel(Pmid), fprintf('%11.2f', AI_loc(imod,k)); end
    fprintf('\n');
end

save(fullfile(PROJ,'FinalResults','FlowDeficitNearP0_Results.mat'), 'Q_norm_all','n_surv_all','conv_all', ...
    'pressures_mmHg','Pmid','AI_loc','D_ad','D_pass','Compens', ...
    'models','labels','P_op','Qreg','KSENSQ');
fprintf('\nResults saved to FlowDeficitNearP0_Results.mat\n');


%% ════════════════════════════════════════════════════════════════════════
%% LOCAL HELPERS  (copied verbatim from TestAutoregulationCurve.m)
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
