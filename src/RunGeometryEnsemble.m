% RunGeometryEnsemble.m
%
% Stochastic geometry ensemble: test whether the loop-
% survival result is statistically robust across DIFFERENT network geometries,
% not just the one canonical network.
%
% For each realization:
%   1. pick a network size (ncombsx varies -> different surface area)
%   2. compute the pial surface area of THAT geometry
%   3. penetrator count = round(2.5 /mm^2 * area)         (area x density rule)
%   4. sample that many penetrators (fresh seed -> stochastic sampling/placement)
%   5. run all 4 coupling models to steady state
%   6. record loop count, loop FRACTION (n_loops / n_loops_total), n_surv, CV(Q)
%
% Because network size (and thus the maximum possible loops) varies between
% realizations, the cross-geometry comparison uses the loop FRACTION.
%
% Aggregates over N realizations: mean +/- SD per model, and a paired
% significance test (Wilcoxon signed-rank across the shared geometries) for the
% two key contrasts:  NoCoupling vs ArtCoupling  and  ArtCoupling vs FlowReg.
%
% Outputs (into FinalResults/):
%   GeometryEnsemble_Results.mat   - all raw + aggregated numbers
%   GeometryEnsemble.png           - loop-fraction distribution per model
%
% RUNTIME: N x 4 models; larger ncombsx networks are slower and FlowReg is the
% bottleneck. Budget a few hours for N=10. Progress is printed per realization.

rng(1234);
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
close all;

PROJ = fileparts(mfilename('fullpath')); addpath(PROJ);
addpath(fullfile(PROJ,'ElectricModel'));

OUTDIR     = 'FinalResults';
r_dead     = 5e-6;
P_op       = 86;
MAX_CHUNKS = 150;
pen_density = 3.35;                        % penetrators / mm^2 (measured: 23 / 6.86 mm^2 crop)

% Ensemble design: 10 realizations, network SIZE varies via ncombsx
% ncombsx rotates through {1,2,3} (different surface areas); each realization
% gets its own random seed so penetrator sampling/placement also varies.
N            = 10;
ncombsx_list = [1 2 3 1 2 3 1 2 3 2];       % size per realization
seeds        = 100 + (1:N);
models = {'NoCoupling','ArtCoupling','FullCoupling','FullCouplingFlowReg'};
nM = numel(models);

% storage [N x nM] -- loops only (per request)
loops_abs  = nan(N,nM);
loops_frac = nan(N,nM);
surv_abs   = nan(N,nM);
area_mm2   = nan(N,1);
npen_used  = nan(N,1);
loops_max  = nan(N,1);

ElementLength = 2e-3;
A_hex = (3*sqrt(3)/2) * ElementLength^2;    % area of one hexagon [m^2]

t_start = tic;
for ir = 1:N
    ncombsx = ncombsx_list(ir);
    rng(seeds(ir));
    branchlevel = 4;
    ncombsy = 2^(branchlevel-2);            % derived as in defaultPars
    n_hex   = ncombsx * ncombsy;
    A_tot   = n_hex * A_hex * 1e6;           % mm^2
    n_pen   = round(pen_density * A_tot);
    area_mm2(ir) = A_tot;  npen_used(ir) = n_pen;

    fprintf('\n%s\n[Realization %d/%d]  ncombsx=%d  area=%.1f mm^2  n_pen=%d  seed=%d\n%s\n', ...
        repmat('=',1,70), ir, N, ncombsx, A_tot, n_pen, seeds(ir), repmat('=',1,70));

    % build this geometry
    Sin = base_Sin(ncombsx, branchlevel);
    S_ref = defaultPars('leptomeningeal_geometry', Sin);
    pen   = load_penetrator_data('penetrators.csv');
    sp    = tr_sample(pen, n_pen);
    S_ref = integrate_sampled_penetrators(S_ref, sp);
    n_new = S_ref.nSE - length(S_ref.sources);
    S_ref.sources = [S_ref.sources, zeros(1, n_new)];

    pial_idx = find([S_ref.IE.pial] == 1);
    sink_idx = find(~S_ref.sources);
    [nloops_total,~] = pial_cycle_rank(S_ref, pial_idx);
    loops_max(ir) = nloops_total;
    fprintf('  n_pial=%d  n_pens=%d  n_loops_total=%d\n', ...
        numel(pial_idx), numel(sink_idx), nloops_total);

    % calibrate Qreg from this geometry's NoCoupling SS
    [Rc, Sc, ~] = tr_run(tr_setP(S_ref, P_op*133), 0, 0.1, 'NoCoupling', [], 50);
    rc = Rc.X(end,:); St = Sc; [St.IE.r] = vout(rc);
    St = calcConductance2025(St); St = solvehemodyn2025(St);
    Qreg = mean(abs([St.SE(sink_idx).Qs]));

    % run the 4 models
    for im = 1:nM
        t0 = tic;
        M = run_ss(S_ref, models{im}, Qreg, P_op*133, pial_idx, sink_idx, r_dead, MAX_CHUNKS);
        [nl,~] = pial_cycle_rank(S_ref, M.alive_segs);
        loops_abs(ir,im)  = nl;
        loops_frac(ir,im) = nl / max(nloops_total,1);
        surv_abs(ir,im)   = numel(M.alive_segs);
        fprintf('   %-20s loops=%2d/%2d (%.2f)  surv=%3d  [%.0fs conv=%d]\n', ...
            models{im}, nl, nloops_total, loops_frac(ir,im), surv_abs(ir,im), ...
            toc(t0), M.converged);
    end
end
fprintf('\nEnsemble finished in %.1f min\n', toc(t_start)/60);

%% Aggregate
mean_frac = mean(loops_frac,1);   sd_frac = std(loops_frac,0,1);
mean_abs  = mean(loops_abs,1);    sd_abs  = std(loops_abs,0,1);

fprintf('\n%s\nGEOMETRY ENSEMBLE SUMMARY  (N=%d realizations, size-varying)\n%s\n', ...
    repmat('=',1,70), N, repmat('=',1,70));
fprintf('%-22s %16s %16s\n','Model','loops (abs)','loop fraction');
for im = 1:nM
    fprintf('%-22s %7.1f +/- %-5.1f %8.2f +/- %-5.2f\n', ...
        models{im}, mean_abs(im), sd_abs(im), mean_frac(im), sd_frac(im));
end

% paired significance tests (across the shared geometries)
fprintf('\nPaired significance (Wilcoxon signed-rank on loop fraction, N=%d):\n', N);
sig_pair('NoCoupling vs ArtCoupling', loops_frac(:,1), loops_frac(:,2));
sig_pair('ArtCoupling vs FullCoupling', loops_frac(:,2), loops_frac(:,3));
sig_pair('ArtCoupling vs FlowReg',     loops_frac(:,2), loops_frac(:,4));

%% Figure: loop-fraction distribution per model
colors = [0.15 0.35 0.75; 0.85 0.45 0.05; 0.55 0.10 0.55; 0.05 0.60 0.20];
figure('Color','w','Position',[150 200 720 420]); hold on;
for im = 1:nM
    xj = im + 0.12*(rand(N,1)-0.5);
    scatter(xj, loops_frac(:,im), 32, colors(im,:), 'filled', ...
        'MarkerFaceAlpha',0.7);                                  % individual realizations
    bar(im, mean_frac(im), 0.5, 'FaceColor', colors(im,:), 'FaceAlpha', 0.25, 'EdgeColor','none');
    errorbar(im, mean_frac(im), sd_frac(im), 'k', 'LineWidth', 1.3, 'CapSize', 10);
end
set(gca,'XTick',1:nM,'XTickLabel',models,'XTickLabelRotation',15);
ylabel('surviving loop fraction  (n\_loops / n\_loops\_total)');
title(sprintf('Loop survival across %d random geometries (size-varying)', N));
ylim([0 1.05]); box off; grid on;
saveas(gcf, fullfile(OUTDIR,'GeometryEnsemble.png'));

save(fullfile(OUTDIR,'GeometryEnsemble_Results.mat'), ...
    'loops_abs','loops_frac','surv_abs','area_mm2','npen_used','loops_max', ...
    'models','ncombsx_list','seeds','N','mean_frac','sd_frac','mean_abs','sd_abs');
fprintf('\nSaved: %s/GeometryEnsemble.png  and  %s/GeometryEnsemble_Results.mat\n', OUTDIR, OUTDIR);

%% ════════════════════════════════════════════════════════════════════════
%% LOCAL HELPERS
%% ════════════════════════════════════════════════════════════════════════
function Sin = base_Sin(ncombsx, branchlevel)
    Sin = struct;
    Sin.npenetrator          = 1;
    Sin.ncombsx              = ncombsx;
    Sin.branchlevel          = branchlevel;
    Sin.balance_Gsin_Gsout   = 10;
    Sin.Vref                 = -0.05;
    Sin.gainRadius2Vref      = 10;
    Sin.ksens                = 0.4;
    Sin.kReg                 = 0.01;
    Sin.CapillaryPotential   = -0.06;
    Sin.CapillaryConductance = 2e-7;
    Sin.kSensQ               = 0.4;
    Sin.gCapBase             = 0.2;
    Sin.ECapBase             = -0.04;
    Sin.DoCapillaryCurrent   = 0;
    Sin.ECresistivity        = 0.1;
    Sin.KeepDynamics         = 0;
    Sin.GsSink               = 1.6e-16;
end

function sig_pair(label, a, b)
    d = a - b;
    if all(d==0)
        fprintf('  %-28s  identical (all diffs 0)\n', label); return;
    end
    if exist('signrank','file')
        p = signrank(a,b);
        fprintf('  %-28s  mean diff=%+.3f   p=%.4f %s\n', label, mean(d), p, star(p));
    else
        % fallback: paired t-test by hand
        t = mean(d)/(std(d)/sqrt(numel(d)));
        fprintf('  %-28s  mean diff=%+.3f   t=%.2f (install stats tb for p)\n', label, mean(d), t);
    end
end

function s = star(p)
    if p<0.001, s='***'; elseif p<0.01, s='**'; elseif p<0.05, s='*'; else, s='(n.s.)'; end
end

function M = run_ss(S_ref, model_name, Qreg, Psrc, pial_idx, sink_idx, r_dead, maxchunks) %#ok<INUSD>
    % loops only: sink_idx kept in the signature for call-site symmetry, unused here
    S = tr_setP(S_ref, Psrc);
    S.gCapBase = 0.05;  S.kSensQ = 0.6;
    S.GscCapillary = 2e-5;  S.kSmCapillary = 7e-6;
    [R, ~, conv] = tr_run(S, Qreg, 0.1, model_name, [], maxchunks);
    r_ss = R.X(end,:);
    M.alive_segs = pial_idx(r_ss(pial_idx) > r_dead);
    M.converged  = conv.converged;
end

function sp = tr_sample(pen, n)
    idx = randi(length(pen.r), n, 1);
    sp.r = pen.r(idx);  sp.l = pen.l(idx);
    sp.R = pen.R(idx);  sp.G = 1./pen.R(idx);
end

function [cycle_rank, info] = pial_cycle_rank(S, seg_list)
    m = length(seg_list);
    if m == 0, cycle_rank=0; info=struct('m',0,'n',0,'c',0); return; end
    all_nodes = unique([[S.IE(seg_list).nodes]]);
    n = length(all_nodes);
    node_map = zeros(max(all_nodes),1);
    node_map(all_nodes) = 1:n;
    adj = cell(n,1);
    for k = 1:m
        ns = S.IE(seg_list(k)).nodes;
        a = node_map(ns(1)); b = node_map(ns(2));
        adj{a}(end+1)=b; adj{b}(end+1)=a;
    end
    visited=false(n,1); c=0;
    for s=1:n
        if visited(s), continue; end
        c=c+1; q=s; visited(s)=true;
        while ~isempty(q)
            v=q(1); q(1)=[];
            for nb=adj{v}
                if ~visited(nb), visited(nb)=true; q(end+1)=nb; end %#ok<AGROW>
            end
        end
    end
    cycle_rank = m-n+c;
    info = struct('m',m,'n',n,'c',c);
end

function [R, S_out, conv] = tr_run(S_geo, QrefCap, ECr, model_name, r0_init, max_chunks_in)
    chunk_size = 1e5;  tol_drdt = 1e-12;
    max_chunks = 150;
    if nargin>=6 && ~isempty(max_chunks_in), max_chunks=max_chunks_in; end
    S = S_geo;
    S.ECresistivity = ECr;  S.QrefCapillary = QrefCap;
    pen_ie = find(~[S.IE.pial]);
    S.kReg = S.kReg * ones(S.nIE,1);  S.kReg(pen_ie) = 0;
    if nargin>=5 && ~isempty(r0_init)
        S.r0 = r0_init(:)';
    else
        S.r0 = 100e-6*ones(1,S.nIE);  S.r0(pen_ie) = 15e-6;
    end
    S.KeepDynamics = 0;
    nc = any(strcmp(model_name,{'FullCoupling','FullCouplingFlowReg'}));
    S.DoCapillaryCurrent = nc;
    S = MakeElectricCircuit(S,0);
    S.DoElectricCoupling       = ~strcmp(model_name,'NoCoupling');
    S.DoCapillaryCurrent       = nc;
    S.FlowRegulatedCapillaries = strcmp(model_name,'FullCouplingFlowReg');
    S.annot = model_name;
    opts = odeset('Events',@tr_nan,'RelTol',1e-6,'AbsTol',1e-9);
    t=0; y=S.r0(:); converged=false; n_ch=0; drdt=NaN;
    while n_ch < max_chunks
        [tc,Xc,tE] = ode45(@(t,X) rdotfunCoupling(t,X,S),[t,t+chunk_size],y,opts);
        n_ch = n_ch+1;
        if ~isempty(tE), y=Xc(end,:)'; break; end
        drdt = max(abs(Xc(end,:)-Xc(1,:)))/(tc(end)-tc(1));
        if drdt<tol_drdt, converged=true; y=Xc(end,:)'; break; end
        t=tc(end); y=Xc(end,:)';
    end
    R.X = y';  R.t = t;
    S_out = S;
    conv = struct('converged',converged,'n_chunks',n_ch,'drdt_final',drdt);
end

function [v,i,d] = tr_nan(~,X)
    v=all(~isnan(X)); i=1; d=0;
end

function S_out = tr_setP(S, Psrc)
    S_out=S; S_out.sourceP=Psrc;
    for k=find(S_out.sources), S_out.SE(k).Ps=Psrc; end
end
