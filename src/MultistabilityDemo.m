% MultistabilityDemo.m
%
% Multistability analysis: are there multiple stable states?
%
% The adaptation system is nonlinear. This asks whether the FINAL steady state
% depends on the INITIAL vessel diameters. Different initial pial diameters are
% run to steady state; we record BOTH the loop count AND the full final pial
% radius vector. If different starts settle to different final diameters (and/or
% loop counts), the system is multistable (initial conditions govern the
% outcome) -- deterministic multiple basins, not chaos.
%
% NOTE: distinct loop counts imply distinct diameters, but the same loop count
% can still hide different diameter configurations -- so we cluster the final
% states in DIAMETER space (not just by loop count).
%
% Run for BOTH ArtCoupling (expected ~single deep attractor) and
% FullCouplingFlowReg (expected multistable) to contrast the two.
%
% Output: FinalResults/MultistabilityDemo.png + _Results.mat

rng(1234);
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
close all;
PROJ = fileparts(mfilename('fullpath'));
addpath(PROJ); addpath(fullfile(PROJ,'ElectricModel'));

r_dead = 5e-6;  P_op = 86;  MAX_CHUNKS = 150;  density = 3.35;  ECr = 0.1;
STATE_TOL_UM = 1.0;   % two final states counted as "the same" if max|dr| < this

fprintf('Building network...\n');
[S_ref, pial_idx, ~, Qreg, nloops_total, A_mm2, npen] = sens_build(2, 2e-3, density, 1234);
np = numel(pial_idx);
fprintf('  area=%.1f mm^2  n_pen=%d  n_pial=%d  n_loops_total=%d\n\n', A_mm2, npen, np, nloops_total);

% initial-condition sets (pial diameters)
uniform_um = [100 60 30 15];
n_random   = 4;
init_labels = {};  init_r0 = {};
for u = uniform_um
    init_labels{end+1} = sprintf('uniform %dum', u);   %#ok<SAGROW>
    init_r0{end+1}     = u*1e-6;
end
for k = 1:n_random
    rng(2000+k);
    init_labels{end+1} = sprintf('random #%d', k);      %#ok<SAGROW>
    init_r0{end+1}     = (10 + 140*rand(1,np))*1e-6;
end
nI = numel(init_r0);

models = {'ArtCoupling','FullCouplingFlowReg'};
nM = numel(models);
loops  = nan(nI, nM);
r_mean = nan(nI, nM);
r_max  = nan(nI, nM);
R_pial = nan(nI, np, nM);     % full final pial radii [um]
conv_f = nan(nI, nM);

for im = 1:nM
    fprintf('%s\nModel: %s\n%s\n', repmat('=',1,64), models{im}, repmat('=',1,64));
    for ii = 1:nI
        t0 = tic;
        [nl,~,cv,r_ss] = sens_ss(S_ref, models{im}, Qreg, ECr, P_op*133, pial_idx, r_dead, MAX_CHUNKS, init_r0{ii});
        rp = r_ss(pial_idx)*1e6;                       % final pial radii [um]
        loops(ii,im)=nl; r_mean(ii,im)=mean(rp); r_max(ii,im)=max(rp);
        R_pial(ii,:,im)=rp; conv_f(ii,im)=cv;
        fprintf('  %-14s -> loops=%2d/%2d  r_mean=%5.1f  r_max=%5.1f um  [%.0fs conv=%d]\n', ...
            init_labels{ii}, nl, nloops_total, mean(rp), max(rp), toc(t0), cv);
    end

    % cluster final states in DIAMETER space (greedy, tolerance STATE_TOL_UM)
    reps = [];                                         % indices of cluster representatives
    cid  = zeros(nI,1);
    for ii = 1:nI
        placed = false;
        for c = 1:numel(reps)
            if max(abs(R_pial(ii,:,im) - R_pial(reps(c),:,im))) < STATE_TOL_UM
                cid(ii)=c; placed=true; break;
            end
        end
        if ~placed, reps(end+1)=ii; cid(ii)=numel(reps); end %#ok<SAGROW>
    end
    nStates = numel(reps);
    fprintf('  distinct final DIAMETER states: %d  (loop counts of reps: %s)  -> %s\n\n', ...
        nStates, mat2str(loops(reps,im)'), ternary(nStates>1,'MULTISTABLE','single attractor'));
    clusters{im} = cid; %#ok<SAGROW>
    nstates_all(im) = nStates; %#ok<SAGROW>
end

% summary table
fprintf('%s\nMULTISTABILITY SUMMARY (loops / r_mean / r_max per initial condition)\n%s\n', ...
    repmat('=',1,78), repmat('=',1,78));
for im=1:nM
    fprintf('\n%s  (%d distinct diameter states)\n', models{im}, nstates_all(im));
    fprintf('  %-14s %8s %10s %10s %8s\n','init','loops','r_mean','r_max','state#');
    for ii=1:nI
        fprintf('  %-14s %5d/%-2d %9.1f %10.1f %6d\n', init_labels{ii}, loops(ii,im), ...
            nloops_total, r_mean(ii,im), r_max(ii,im), clusters{im}(ii));
    end
end

% figures
mc = [0.85 0.45 0.05; 0.05 0.60 0.20];
figure('Color','w','Position',[100 150 1100 430]);

% (1) loops vs init
subplot(1,3,1); hold on;
for im=1:nM
    plot(1:nI, loops(:,im), 'o-','LineWidth',1.6,'Color',mc(im,:), ...
        'MarkerFaceColor',mc(im,:),'MarkerSize',7,'DisplayName',models{im});
end
yline(nloops_total,'k:','max','HandleVisibility','off');
set(gca,'XTick',1:nI,'XTickLabel',init_labels,'XTickLabelRotation',30);
ylabel('final loops'); title('Loops vs initial condition'); legend('Location','best');
grid on; box off; ylim([-1 nloops_total+2]);

% (2,3) final pial radius profiles per init, one panel per model
%       (sorted descending; overlapping curves = same attractor)
for im=1:nM
    subplot(1,3,1+im); hold on;
    cmap = turbo(nI);
    for ii=1:nI
        plot(sort(R_pial(ii,:,im),'descend'), 'LineWidth',1.2, 'Color',cmap(ii,:), ...
            'DisplayName',init_labels{ii});
    end
    yline(r_dead*1e6,'k--','r_{dead}','HandleVisibility','off');
    xlabel('pial segment (sorted)'); ylabel('final radius [\mum]');
    title(sprintf('%s: final diameters (%d states)', models{im}, nstates_all(im)));
    grid on; box off;
    if im==nM, legend('Location','northeast','FontSize',7); end
end
sgtitle('Multistability: initial diameters \rightarrow final state','FontWeight','bold');
saveas(gcf, fullfile(PROJ,'FinalResults','MultistabilityDemo.png'));
save(fullfile(PROJ,'FinalResults','MultistabilityDemo_Results.mat'), ...
    'loops','r_mean','r_max','R_pial','conv_f','clusters','nstates_all', ...
    'init_labels','models','nloops_total','uniform_um','n_random','ECr','STATE_TOL_UM','pial_idx');
fprintf('\nSaved: FinalResults/MultistabilityDemo.png and _Results.mat\n');

function s = ternary(c,a,b), if c, s=a; else, s=b; end, end
