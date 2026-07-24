% TestPhysiologicalPenetrators.m
%
% Tests all four coupling models with a physiologically correct penetrator
% count, derived from the actual network area.
%
% Network area calculation:
%   branchlevel=4 -> ncombsy=4, ncombsx=2 -> 8 hexagons
%   ElementLength=2mm -> hexagon side = 2mm
%   Area per hexagon = (3*sqrt(3)/2) * (2e-3)^2 = 10.39 mm^2
%   Total area = 8 * 10.39 = 83.1 mm^2
%   At 2.5 penetrators/mm^2 -> ~208 penetrators needed
%   Currently: 20 penetrators -> 0.24/mm^2 (10x too sparse)
%
% Fix: sample penetrators WITH REPLACEMENT from the 23 available in
% penetrators.csv to reach the physiological count.
%
% Runs all 4 models at 86mmHg and reports:
%   - n_surv, Q per penetrator, convergence
%   - Quick 70/86/100 mmHg Q-P sweep for each model

rng(1234);
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
close all;

PROJ = fileparts(mfilename('fullpath')); addpath(PROJ);
addpath(fullfile(PROJ,'ElectricModel'));

r_dead = 5e-6;

%% Compute physiological penetrator count
ElementLength = 2e-3;    % m -- from defaultPars.m
ncombsx       = 2;
ncombsy       = 4;       % 2^(branchlevel-2) = 2^2
n_hex         = ncombsx * ncombsy;
A_hex         = (3*sqrt(3)/2) * ElementLength^2;
A_total_mm2   = n_hex * A_hex * 1e6;
pen_density   = 3.35;    % penetrators/mm^2 (measured: 23 / 6.86 mm^2 crop)
n_pen_target  = round(pen_density * A_total_mm2);

fprintf('\n%s\n', repmat('=',1,60));
fprintf('Network geometry\n');
fprintf('%s\n', repmat('=',1,60));
fprintf('  Total area     = %.1f mm^2  (%.1fx%.1f mm equivalent)\n', ...
    A_total_mm2, sqrt(A_total_mm2), sqrt(A_total_mm2));
fprintf('  Density target = %.1f pen/mm^2\n', pen_density);
fprintf('  n_pen_target   = %d\n', n_pen_target);
fprintf('  Previous count = 20  (implied density = %.2f/mm^2)\n\n', 20/A_total_mm2);

%% Build network
fprintf('Building network...\n');
Sin = struct;
Sin.npenetrator        = 1;
Sin.ncombsx            = ncombsx;
Sin.branchlevel        = 4;
Sin.balance_Gsin_Gsout = 10;
Sin.Vref               = -0.05;
Sin.gainRadius2Vref    = 10;
Sin.ksens              = 0.4;
Sin.kReg               = 0.01;
Sin.CapillaryPotential  = -0.06;
Sin.CapillaryConductance = 2e-7;
Sin.kSensQ              = 0.4;
Sin.gCapBase            = 0.2;
Sin.ECapBase            = -0.04;
Sin.DoCapillaryCurrent  = 0;
Sin.ECresistivity       = 0.1;
Sin.KeepDynamics        = 0;
Sin.GsSink              = 1.6e-16;  % 100x smaller than default; keeps R_pial and R_sink comparable

S_ref = defaultPars('leptomeningeal_geometry', Sin);

% Sample WITH REPLACEMENT to reach physiological count
pen  = load_penetrator_data('penetrators.csv');
n_available = length(pen.r);
fprintf('  Penetrators in CSV: %d\n', n_available);
fprintf('  Sampling %d with replacement (bootstrap from %d anatomical measurements)\n\n', ...
    n_pen_target, n_available);

%% Loop count: formula m - n + c (pre-penetrator topology)
% Cycle rank = number of independent loops in the pial graph.
% m = pial segments, n = pial nodes, c = connected components.
% Exact: no enumeration needed, works for any loop length.
%
% The same formula applied to the SURVIVING subgraph gives alive loops.
% integrate_sampled_penetrators splits pial segments, so run on original topology.

sp = pp_sample_with_replacement(pen, n_pen_target);
S_ref = integrate_sampled_penetrators(S_ref, sp);
n_new = S_ref.nSE - length(S_ref.sources);
S_ref.sources = [S_ref.sources, zeros(1, n_new)];

pial_idx = find([S_ref.IE.pial] == 1);
sink_idx = find(~S_ref.sources);
n_pial   = length(pial_idx);
n_sinks  = length(sink_idx);
[n_loops_total, ~] = pial_cycle_rank(S_ref, pial_idx);

fprintf('\n%s\n', repmat('-',1,60));
fprintf('Post-penetrator pial graph: n_loops (m-n+c) = %d\n', n_loops_total);
fprintf('%s\n\n', repmat('-',1,60));

fprintf('Network built: nIE=%d, n_pial=%d (post-split), n_penetrators=%d\n\n', ...
    S_ref.nIE, n_pial, n_sinks);

%% Qref calibration
fprintf('Calibrating Qref from NoCoupling SS at 86 mmHg...\n');
S86 = pp_setP(S_ref, 86*133);
[R_cal, S_cal, conv_cal] = pp_run(S86, 0, 0.1, 'NoCoupling', []);
fprintf('  Calibration: %d chunk(s), max|dr/dt|=%.2e m/s, converged=%d\n', ...
    conv_cal.n_chunks, conv_cal.drdt_final, conv_cal.converged);
r_ss_cal  = R_cal.X(end,:);
S_f = S_cal; [S_f.IE.r] = vout(r_ss_cal);
S_f = calcConductance2025(S_f); S_f = solvehemodyn2025(S_f);
Q_per_pen = abs([S_f.SE(sink_idx).Qs]);
Q_op      = mean(Q_per_pen);
Qreg      = Q_op;                % flow-regulation target = actual mean flow
Qref      = Q_op / (1 - 0.5);   % fd0=0.5, used for adequacy calc only

fprintf('  Q_op per penetrator = %.4f fL/s  (was ~%.0f fL/s with 20 pens)\n', ...
    Q_op*1e15, 131000/n_sinks);
fprintf('  Qref = %.4f fL/s\n\n', Qref*1e15);

%% Run all 4 models at 86 mmHg
models = {'NoCoupling','ArtCoupling','FullCoupling','FullCouplingFlowReg'};
pressures_mmHg = [70 86 100];
nP = length(pressures_mmHg);

fprintf('%s\n', repmat('=',1,65));
fprintf('Running all 4 models\n');
fprintf('%s\n', repmat('=',1,65));

results = struct('model',{},'converged',{},'n_surv_86',{},'n_loops_86',{},'Q_mean',{},'Q_std',{},'Q_cv',{}, ...
                 'slope_passive',{},'slope_adapt',{},'ratio',{},'AI',{},'Qp_pass',{},'Qp_ad',{});

for imod = 1:4
    model_name = models{imod};
    fprintf('\n--- %s ---\n', model_name);

    % 86 mmHg SS
    S86m = pp_setP(S_ref, 86*133);
    S86m.gCapBase      = 0.05;
    S86m.kSensQ        = 0.6;
    S86m.GscCapillary  = 2e-5;
    S86m.kSmCapillary  = 7e-6;
    try
        [R86, S86_out, conv86] = pp_run(S86m, Qreg, 0.1, model_name, []);
        r_ss = R86.X(end,:);
        results(imod).r_ss = r_ss;
        converged = conv86.converged;   % reflect real convergence, not just "did not error"
        fprintf('  Convergence: %d chunk(s), max|dr/dt|=%.2e m/s, converged=%d\n', ...
            conv86.n_chunks, conv86.drdt_final, conv86.converged);
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        results(imod).model = model_name;
        results(imod).converged = 0;
        continue;
    end

    n_surv_86 = sum(r_ss(pial_idx) > r_dead);
    alive_post = pial_idx(r_ss(pial_idx) > r_dead);
    [n_loops_86, loop_info] = pial_cycle_rank(S_ref, alive_post);
    S_f = S86_out; [S_f.IE.r] = vout(r_ss);
    S_f = calcConductance2025(S_f); S_f = solvehemodyn2025(S_f);
    Q_i = abs([S_f.SE(sink_idx).Qs]);
    fprintf('  n_surv pial segs (excl. penetrators) = %d/%d\n', n_surv_86, n_pial);
    fprintf('  n_loops alive = %d/%d  (m=%d, n=%d, c=%d)\n', ...
        n_loops_86, n_loops_total, loop_info.m, loop_info.n, loop_info.c);
    fprintf('  Q per pen: mean=%.4f fL/s  std=%.4f  CV=%.3f\n', ...
        mean(Q_i)*1e15, std(Q_i)*1e15, std(Q_i)/mean(Q_i));

    % Passive reference
    Qp_pass = nan(1,nP);
    for iP = 1:nP
        Qp_pass(iP) = mean(pp_passiveQ(S_ref, r_ss, pressures_mmHg(iP)*133, sink_idx))*1e15;
    end
    slope_pass = (Qp_pass(3)-Qp_pass(1))/(pressures_mmHg(3)-pressures_mmHg(1));

    % Pressure sweep
    Qp_ad = nan(1,nP);
    for iP = 1:nP
        Spm = pp_setP(S_ref, pressures_mmHg(iP)*133);
        Spm.QrefCapillary = Qreg;
        Spm.gCapBase      = 0.05;
        Spm.kSensQ        = 0.6;
        Spm.GscCapillary  = 2e-5;
        Spm.kSmCapillary  = 7e-6;
        try
            [Rpm, ~, conv_pm] = pp_run(Spm, Qreg, 0.1, model_name, r_ss);
            r_fin = Rpm.X(end,:);
            dr_pial = r_fin(pial_idx) - r_ss(pial_idx);
            fprintf('    P=%d mmHg: mean|Δr|=%.3f µm  max|Δr|=%.3f µm\n', ...
                pressures_mmHg(iP), mean(abs(dr_pial))*1e6, max(abs(dr_pial))*1e6);
            S_f2 = Spm; [S_f2.IE.r] = vout(r_fin);
            S_f2 = calcConductance2025(S_f2); S_f2 = solvehemodyn2025(S_f2);
            Qp_ad(iP) = mean(abs([S_f2.SE(sink_idx).Qs]))*1e15;
            if ~conv_pm.converged
                fprintf('  P=%d mmHg: WARNING not converged after %d chunks (max|dr/dt|=%.2e m/s)\n', ...
                    pressures_mmHg(iP), conv_pm.n_chunks, conv_pm.drdt_final);
            end
        catch ME
            fprintf('  P=%d ERROR: %s\n', pressures_mmHg(iP), ME.message);
        end
    end

    slope_ad = (Qp_ad(3)-Qp_ad(1))/(pressures_mmHg(3)-pressures_mmHg(1));
    ratio    = slope_ad / slope_pass;
    AI       = 1 - ratio;   % 0=no regulation, 1=perfect, <0=WSS amplification
    fprintf('  Qp_pass [70 86 100]: %s fL/s\n', mat2str(round(Qp_pass,1)));
    fprintf('  Qp_ad   [70 86 100]: %s fL/s\n', mat2str(round(Qp_ad,1)));
    fprintf('  slope_pass=%.2f  slope_adapt=%.2f  ratio=%.3f  AI=%.3f\n', slope_pass, slope_ad, ratio, AI);

    results(imod).model        = model_name;
    results(imod).converged    = converged;
    results(imod).n_surv_86    = n_surv_86;
    results(imod).n_loops_86   = n_loops_86;
    results(imod).Q_mean       = mean(Q_i)*1e15;
    results(imod).Q_std        = std(Q_i)*1e15;
    results(imod).Q_cv         = std(Q_i)/mean(Q_i);
    results(imod).slope_passive = slope_pass;
    results(imod).slope_adapt  = slope_ad;
    results(imod).ratio        = ratio;
    results(imod).AI           = AI;
    results(imod).Qp_pass      = Qp_pass;
    results(imod).Qp_ad        = Qp_ad;
end

%% Summary
fprintf('\n%s\n', repmat('=',1,85));
fprintf('SUMMARY  (n_pial=%d [pial segs only, excl. penetrators], n_pens=%d, area=%.1f mm^2)\n', ...
    n_pial, n_sinks, A_total_mm2);
fprintf('         n_loops_total=%d  (formula: m-n+c on pre-penetrator topology)\n', n_loops_total);
fprintf('%s\n', repmat('=',1,85));
fprintf('%-25s %5s %8s %8s %8s %8s %8s %8s %8s\n', ...
    'Model','OK?','n_surv','n_loops','Q_mean','CV','slope','ratio','AI');
for i = 1:length(results)
    if results(i).converged
        fprintf('%-25s %5s %8d %8d %8.4f %8.3f %8.2f %8.3f %8.3f\n', ...
            results(i).model, 'YES', results(i).n_surv_86, results(i).n_loops_86, ...
            results(i).Q_mean, results(i).Q_cv, ...
            results(i).slope_adapt, results(i).ratio, results(i).AI);
    else
        fprintf('%-25s %5s  FAILED\n', results(i).model, 'NO');
    end
end

fprintf('\nKey check: if all 4 models converge and n_surv is reasonable,\n');
fprintf('the physiological penetrator count is compatible with the model.\n');

save(fullfile(PROJ,'FinalResults','PhysiologicalPenetrators_Results.mat'), 'results', 'n_pial', 'n_sinks', ...
    'A_total_mm2', 'pen_density', 'Qref', 'Q_op', 'pressures_mmHg', 'n_loops_total', 'pial_idx');
fprintf('\nSaved to PhysiologicalPenetrators_Results.mat\n');


%% ════════════════════════════════════════════════════════════════════════
%% LOCAL HELPERS
%% ════════════════════════════════════════════════════════════════════════

function [cycle_rank, info] = pial_cycle_rank(S, seg_list)
% Compute the cycle rank (= number of independent loops) of the subgraph
% defined by seg_list (a vector of IE indices, all with IE.pial==1).
%
% Formula: cycle_rank = m - n + c
%   m = number of segments in seg_list
%   n = number of unique nodes those segments touch
%   c = number of connected components in that subgraph
%
% Works for any loop length -- no enumeration needed.
% Use on the pre-penetrator pial topology for total loops.
% Use on the surviving subset (r > r_dead) for alive loops.

    m = length(seg_list);
    if m == 0
        cycle_rank = 0;
        info = struct('m',0,'n',0,'c',0);
        return
    end

    all_nodes = unique([[S.IE(seg_list).nodes]]);
    n         = length(all_nodes);

    % BFS to count connected components
    max_nid  = max(all_nodes);
    node_map = zeros(max_nid, 1);
    node_map(all_nodes) = 1:n;

    adj = cell(n, 1);
    for k = 1:m
        ns = S.IE(seg_list(k)).nodes;
        a  = node_map(ns(1));  b = node_map(ns(2));
        adj{a}(end+1) = b;
        adj{b}(end+1) = a;
    end

    visited = false(n, 1);
    c = 0;
    for s = 1:n
        if visited(s), continue; end
        c = c + 1;
        queue = s;  visited(s) = true;
        while ~isempty(queue)
            v = queue(1);  queue(1) = [];
            for nb = adj{v}
                if ~visited(nb)
                    visited(nb) = true;
                    queue(end+1) = nb; %#ok<AGROW>
                end
            end
        end
    end

    cycle_rank = m - n + c;
    info = struct('m', m, 'n', n, 'c', c);
end

function sp = pp_sample_with_replacement(pen, n)
% Sample n penetrators with replacement from pen struct
    n_avail = length(pen.r);
    idx = randi(n_avail, n, 1);
    sp.r = pen.r(idx);
    sp.l = pen.l(idx);
    sp.R = pen.R(idx);
    sp.G = 1 ./ pen.R(idx);
end

function [R, S_out, conv_info] = pp_run(S_geo, QrefCap, ECr, model_name, r0_init)
% Runs the simulation using a chunk loop (from CompareODE45vsChunk.m).
% Integrates in successive 1e5-s windows; stops when the maximum radius
% change across a window falls below tol_drdt = 1e-12 m/s.
% conv_info: struct with fields converged, n_chunks, drdt_final, t_total.
    chunk_size = 1e5;     % [s] per chunk
    tol_drdt   = 1e-12;   % [m/s] tightened from 1e-11 for publication-quality convergence
    max_chunks = 150;     % safety limit; raised from 50 so FullCouplingFlowReg converges

    S = S_geo;
    S.ECresistivity = ECr;
    S.QrefCapillary = QrefCap;

    % Penetrators (pial=0) must not adapt: fixed radius at CSV value.
    % Make kReg a per-element column vector so penetrators get kReg=0.
    pen_ie = find(~[S.IE.pial]);
    S.kReg = S.kReg * ones(S.nIE, 1);
    S.kReg(pen_ie) = 0;

    if nargin >= 5 && ~isempty(r0_init)
        S.r0 = r0_init(:)';
    else
        S.r0 = 100e-6 * ones(1, S.nIE);
        S.r0(pen_ie) = 15e-6;   % penetrators: typical arteriole radius, stays fixed (kReg=0)
    end
    S.KeepDynamics = 0;
    needs_caps = any(strcmp(model_name, {'FullCoupling','FullCouplingFlowReg'}));
    S.DoCapillaryCurrent = needs_caps;
    S = MakeElectricCircuit(S, 0);
    S.DoElectricCoupling       = ~strcmp(model_name, 'NoCoupling');
    S.DoCapillaryCurrent       = needs_caps;
    S.FlowRegulatedCapillaries = strcmp(model_name, 'FullCouplingFlowReg');
    S.annot = model_name;

    opts = odeset('Events', @pp_detectNaN, 'RelTol', 1e-6, 'AbsTol', 1e-9);

    t_now     = 0;
    y_now     = S.r0(:);
    converged = false;
    n_chunks  = 0;
    drdt_last = NaN;

    while n_chunks < max_chunks
        [t_c, X_c, tErr] = ode45(@(t,X) rdotfunCoupling(t, X, S), ...
                                  [t_now, t_now + chunk_size], y_now, opts);
        n_chunks = n_chunks + 1;

        if ~isempty(tErr)
            fprintf('    [chunk %d] NaN detected — stopping early\n', n_chunks);
            y_now = X_c(end,:)';
            break
        end

        dr_chunk   = max(abs(X_c(end,:) - X_c(1,:)));
        drdt_last  = dr_chunk / (t_c(end) - t_c(1));

        if drdt_last < tol_drdt
            converged = true;
            y_now = X_c(end,:)';
            break
        end
        t_now = t_c(end);
        y_now = X_c(end,:)';
    end

    R.X = y_now';          % 1 x nIE, matches RunSimulation output format
    R.t = t_now + chunk_size * converged;

    conv_info.converged  = converged;
    conv_info.n_chunks   = n_chunks;
    conv_info.drdt_final = drdt_last;   % [m/s]
    conv_info.t_total    = t_now + chunk_size * converged;

    S_out = S;
end

function [value, isterminal, direction] = pp_detectNaN(~, X)
    value      = all(~isnan(X));
    isterminal = 1;
    direction  = 0;
end

function S_out = pp_setP(S, Psrc)
    S_out = S;
    S_out.sourceP = Psrc;
    for k = find(S_out.sources)
        S_out.SE(k).Ps = Psrc;
    end
end

function Q = pp_passiveQ(S, r_fixed, Psrc, sink_idx)
    S = pp_setP(S, Psrc);
    S.DoElectricCoupling = 0;
    S.DoCapillaryCurrent = 0;
    S.FlowRegulatedCapillaries = 0;
    [S.IE.r] = vout(r_fixed);
    S = calcConductance2025(S);
    S = solvehemodyn2025(S);
    Q = abs([S.SE(sink_idx).Qs]);
end
