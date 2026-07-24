% SensitivityLengthLambda.m
%
% Sensitivity analysis: length vs lambda.
%
% Tests whether the dimensionless ratio lambda/L (electrotonic length constant /
% segment length) is the governing stability parameter, by varying BOTH the
% segment length L (ElementLength) and lambda (via ECresistivity) and checking
% whether loop survival collapses onto a single curve when plotted against
% lambda/L.
%
%   lambda = sqrt( ECheight / (2*(gNa+gLeak+gK0)*rho') )
%   L      = ElementLength (pial arcade edge)
%
% Changing L rebuilds the network (area, and thus penetrator count at 3.35/mm^2,
% change), so loop survival is compared as a FRACTION (n_loops / n_loops_total).
% Models: ArtCoupling and FullCouplingFlowReg.
%
% Output: FinalResults/SensitivityLengthLambda.png + _Results.mat
%
% RUNTIME: (#lengths) x (#ECr) x (#models) steady states + one build per length.

rng(1234);
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
close all;
PROJ = fileparts(mfilename('fullpath'));
addpath(PROJ); addpath(fullfile(PROJ,'ElectricModel'));

r_dead = 5e-6;  P_op = 86;  MAX_CHUNKS = 150;  density = 3.35;

L_list   = [1e-3 2e-3 3e-3];               % segment lengths [m]
ECr_list = [1 0.3 0.1 0.03];               % -> lambda
gm0 = 0.2 + 0.2 + 0.5;                      % gNa+gLeak+gK0 (resting)
models = {'ArtCoupling','FullCouplingFlowReg'};
nL = numel(L_list); nE = numel(ECr_list); nM = numel(models);

% storage
frac    = nan(nL, nE, nM);
loopsA  = nan(nL, nE, nM);
lamL    = nan(nL, nE);
maxloop = nan(nL,1);  areaL = nan(nL,1);  penL = nan(nL,1);

for il = 1:nL
    L = L_list(il);
    fprintf('%s\n[Length %d/%d]  L = %.1f mm\n%s\n', repmat('=',1,60), il, nL, L*1e3, repmat('=',1,60));
    [S_ref, pial_idx, ~, Qreg, nloops_total, A_mm2, npen] = sens_build(2, L, density, 1234);
    maxloop(il)=nloops_total; areaL(il)=A_mm2; penL(il)=npen;
    fprintf('  area=%.1f mm^2  n_pen=%d  n_loops_total=%d\n', A_mm2, npen, nloops_total);
    for ie = 1:nE
        ECr = ECr_list(ie);
        lam = sqrt(1e-6/(2*gm0*ECr));      % [m]
        lamL(il,ie) = lam / L;
        fprintf(' ECr=%.3g  lambda=%.2fmm  lambda/L=%.2f\n', ECr, lam*1e3, lamL(il,ie));
        for im = 1:nM
            t0=tic;
            [nl,~,cv] = sens_ss(S_ref, models{im}, Qreg, ECr, P_op*133, pial_idx, r_dead, MAX_CHUNKS, []);
            loopsA(il,ie,im) = nl;
            frac(il,ie,im)   = nl / max(nloops_total,1);
            fprintf('   %-22s loops=%2d/%2d frac=%.2f [%.0fs conv=%d]\n', ...
                models{im}, nl, nloops_total, frac(il,ie,im), toc(t0), cv);
        end
    end
end

%% figure: loop fraction vs lambda/L, per model, colored by L
Lcolors = lines(nL);
figure('Color','w','Position',[130 180 980 430]);
for im=1:nM
    subplot(1,nM,im); hold on;
    for il=1:nL
        x = lamL(il,:);  y = squeeze(frac(il,:,im));
        [x,ord]=sort(x); y=y(ord);
        plot(x, y, '-o', 'LineWidth',1.7, 'Color',Lcolors(il,:), ...
            'MarkerFaceColor',Lcolors(il,:), 'DisplayName',sprintf('L=%.0fmm',L_list(il)*1e3));
    end
    xline(1,'k--','\lambda/L=1','HandleVisibility','off');
    xlabel('\lambda / L'); ylabel('surviving loop fraction');
    title(models{im}); ylim([-0.02 1.05]); grid on; box off;
    if im==1, legend('Location','southeast'); end
end
sgtitle('Length vs \lambda : does loop survival collapse onto \lambda/L ?','FontWeight','bold');
saveas(gcf, fullfile(PROJ,'FinalResults','SensitivityLengthLambda.png'));
save(fullfile(PROJ,'FinalResults','SensitivityLengthLambda_Results.mat'), ...
    'L_list','ECr_list','frac','loopsA','lamL','maxloop','areaL','penL','models');
fprintf('\nSaved: FinalResults/SensitivityLengthLambda.png and _Results.mat\n');
