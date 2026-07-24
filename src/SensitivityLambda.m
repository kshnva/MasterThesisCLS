% SensitivityLambda.m
%
% Sensitivity analysis: vary the electrotonic length constant
% lambda and see how loop survival responds.
%
% lambda is set by the EC resistivity:
%     lambda = sqrt( ECheight / (2*(gNa+gLeak+gK0)*ECresistivity) )
% so sweeping ECresistivity sweeps lambda. Pial segment length L = ElementLength.
% Phase-1 found loops stabilise when lambda/L >= 1; here we re-test that on the
% corrected Phase-2 network (physiological penetrators, 3.35/mm^2).
%
% Models: ArtCoupling and FullCouplingFlowReg (lambda acts through EC coupling).
% Output: FinalResults/SensitivityLambda.png + printed table.

rng(1234);
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
close all;
PROJ = fileparts(mfilename('fullpath'));
addpath(PROJ); addpath(fullfile(PROJ,'ElectricModel'));

r_dead     = 5e-6;
P_op       = 86;
MAX_CHUNKS = 150;
density    = 3.35;
ElementLength = 2e-3;                 % pial segment length L [m]

% build the network once (lambda does not change geometry)
fprintf('Building network...\n');
[S_ref, pial_idx, sink_idx, Qreg, nloops_total, A_mm2, npen] = ...
    sens_build(2, ElementLength, density, 1234);
fprintf('  area=%.1f mm^2  n_pen=%d  n_pial=%d  n_loops_total=%d\n\n', ...
    A_mm2, npen, numel(pial_idx), nloops_total);

% ECresistivity sweep -> lambda
ECr_list = [3 1 0.5 0.3 0.2 0.15 0.1 0.05 0.03 0.01];
gm0 = 0.2 + 0.2 + 0.5;                                  % gNa+gLeak+gK0 (resting)
lam = sqrt(1e-6 ./ (2*gm0*ECr_list));                  % lambda [m], ECheight=1e-6
lam_mm = lam*1e3;  lamL = lam/ElementLength;

models = {'ArtCoupling','FullCouplingFlowReg'};
nM = numel(models);
loops = nan(numel(ECr_list), nM);

for ie = 1:numel(ECr_list)
    ECr = ECr_list(ie);
    fprintf('ECresistivity=%.3g  lambda=%.2f mm  lambda/L=%.2f\n', ECr, lam_mm(ie), lamL(ie));
    for im = 1:nM
        t0 = tic;
        [nl,~,cv] = sens_ss(S_ref, models{im}, Qreg, ECr, P_op*133, pial_idx, r_dead, MAX_CHUNKS, []);
        loops(ie,im) = nl;
        fprintf('   %-22s loops=%2d/%2d  [%.0fs conv=%d]\n', models{im}, nl, nloops_total, toc(t0), cv);
    end
end

% table
fprintf('\n%s\nLOOPS vs LAMBDA (L=%.1f mm)\n%s\n', repmat('=',1,64), ElementLength*1e3, repmat('=',1,64));
fprintf('%-8s %-9s %-8s', 'ECresist','lambda(mm)','lambda/L');
for im=1:nM, fprintf(' %20s', models{im}); end
fprintf('\n');
for ie=1:numel(ECr_list)
    fprintf('%-8.3g %-9.2f %-8.2f', ECr_list(ie), lam_mm(ie), lamL(ie));
    for im=1:nM, fprintf(' %14d/%-5d', loops(ie,im), nloops_total); end
    fprintf('\n');
end

% figure: loops vs lambda/L
colors = [0.85 0.45 0.05; 0.05 0.60 0.20];
figure('Color','w','Position',[150 200 720 430]); hold on;
for im=1:nM
    plot(lamL, loops(:,im), '-o', 'LineWidth',1.9, 'Color',colors(im,:), ...
        'MarkerFaceColor',colors(im,:), 'DisplayName',models{im});
end
xline(1,'k--','\lambda/L = 1','LabelVerticalAlignment','bottom','HandleVisibility','off');
yline(nloops_total,'k:','max','HandleVisibility','off');
xlabel('\lambda / L   (electrotonic length / segment length)');
ylabel('surviving loops'); title('Loop survival vs electrotonic length constant \lambda');
set(gca,'XScale','log'); grid on; box off; legend('Location','southeast');
saveas(gcf, fullfile(PROJ,'FinalResults','SensitivityLambda.png'));
save(fullfile(PROJ,'FinalResults','SensitivityLambda_Results.mat'), ...
    'ECr_list','lam','lam_mm','lamL','loops','models','nloops_total','ElementLength', ...
    'A_mm2','npen');
fprintf('\nSaved: FinalResults/SensitivityLambda.png and SensitivityLambda_Results.mat\n');
