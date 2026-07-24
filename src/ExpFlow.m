% ExpFlow.m
%
% Experiment: penetrator perfusion at the operating point AND under occlusion,
% for the 4 coupling models. Perfusion and occlusion share the same baseline
% steady state, so both are derived from a single SS solve per model.
% (Replaces the separate ExpPerfusion.m + ExpOcclusion.m; density configurable.)
%
% Per model:
%   baseline SS at 86 mmHg
%     -> perfusion : CV of penetrator flow, flow distribution, radius distribution
%   occlude one feeding source, RE-ADAPT from baseline radii
%     -> occlusion : % total flow maintained, % penetrators >= 50% of own baseline
%
% Set n_src > 1 to occlude several different sources and average (robustness).
%
% Output: FinalResults/Exp_Flow.png + _Results.mat

rng(1234);
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
close all;
PROJ = fileparts(mfilename('fullpath'));
addpath(PROJ); addpath(fullfile(PROJ,'ElectricModel'));

r_dead = 5e-6;  P_op = 86;  MAX_CHUNKS = 150;
DENSITY = 3.35;
n_src   = 1;                 % # sources to occlude (robustness); 1 = original single-source

fprintf('Building network (density=%.2f /mm^2)...\n', DENSITY);
[S_ref, pial_idx, sink_idx, Qreg, nloops_total, A_mm2, npen] = sens_build(2, 2e-3, DENSITY, 1234);
Qref = Qreg / (1 - 0.5);                  % = 2*Q_mean, adequacy reference only
src_idx = find(S_ref.sources);
occ_srcs = src_idx(1:min(n_src, numel(src_idx)));
fprintf('  area=%.1f mm^2  n_pen=%d  n_sources=%d (occluding %d)  Qreg=%.3f fL/s\n\n', ...
    A_mm2, npen, numel(src_idx), numel(occ_srcs), Qreg*1e15);

models = {'NoCoupling','ArtCoupling','FullCoupling','FullCouplingFlowReg'};
labels = {'No coupling','Art. coupling','Full coupling','Full + flow reg.'};
colors = [0.15 0.35 0.75; 0.85 0.45 0.05; 0.55 0.10 0.55; 0.05 0.60 0.20];
nM = numel(models);

CVq=nan(1,nM); Qmean=nan(1,nM); Qdist=cell(1,nM); r_pial=cell(1,nM);
flow_maint=nan(numel(occ_srcs),nM); adequacy=nan(numel(occ_srcs),nM); Qocc_dist=cell(1,nM);

for im = 1:nM
    fprintf('--- %s ---\n', models{im});
    % baseline SS (shared by perfusion + occlusion)
    [~,~,cvb,r_ss] = sens_ss(S_ref, models{im}, Qreg, 0.1, P_op*133, pial_idx, r_dead, MAX_CHUNKS, []);
    Q_base = solve_flows(S_ref, r_ss, P_op*133, sink_idx);
    % perfusion
    Qdist{im}=Q_base*1e15; Qmean(im)=mean(Q_base)*1e15; CVq(im)=std(Q_base)/mean(Q_base);
    r_pial{im}=r_ss(pial_idx)*1e6;
    fprintf('   perfusion : Q_mean=%.1f fL/s  CV(Q)=%.3f  r_mean=%.1f um  [conv=%d]\n', ...
        Qmean(im), CVq(im), mean(r_pial{im}), cvb);
    % occlusion (one or more sources)
    for is = 1:numel(occ_srcs)
        S_occ = occlude(S_ref, occ_srcs(is));
        [~,~,cvo,r_occ] = sens_ss(S_occ, models{im}, Qreg, 0.1, P_op*133, pial_idx, r_dead, MAX_CHUNKS, r_ss(pial_idx));
        Q_occ = solve_flows(S_occ, r_occ, P_op*133, sink_idx);
        flow_maint(is,im) = 100*sum(Q_occ)/sum(Q_base);
        adequacy(is,im)   = 100*mean(Q_occ >= 0.5*Q_base);
        if is==1, Qocc_dist{im}=Q_occ*1e15; end
        fprintf('   occlusion src %d: flow maintained=%.1f%%  adequacy=%.1f%%  [conv=%d]\n', ...
            occ_srcs(is), flow_maint(is,im), adequacy(is,im), cvo);
    end
end
fm = mean(flow_maint,1);  ad = mean(adequacy,1);

%% summary
fprintf('\n%s\nFLOW SUMMARY (density=%.2f /mm^2, n_pen=%d)\n%s\n', repmat('=',1,66), DENSITY, npen, repmat('=',1,66));
fprintf('%-22s %8s %8s %14s %14s\n','Model','Q_mean','CV(Q)','flow maint.','adequacy');
for im=1:nM
    fprintf('%-22s %8.1f %8.3f %12.1f%% %12.1f%%\n', models{im}, Qmean(im), CVq(im), fm(im), ad(im));
end

%% figure (2x2)
figure('Color','w','Position',[80 120 980 720]);

subplot(2,2,1);
b=bar(1:nM, CVq,'FaceColor','flat','EdgeColor','none'); b.CData=colors;
set(gca,'XTick',1:nM,'XTickLabel',labels,'XTickLabelRotation',20);
ylabel('CV of penetrator flow'); title('Perfusion heterogeneity (lower = more uniform)'); box off;

subplot(2,2,2); hold on;
for im=1:nM
    Q=Qdist{im}; Q=Q(Q<=5000); [f,x]=ksdensity(Q);
    plot(x,f,'LineWidth',1.7,'Color',colors(im,:),'DisplayName',labels{im});
end
xlabel('penetrator flow [fL/s]'); ylabel('density'); title('Baseline flow distribution');
legend('Location','northeast','FontSize',8); box off; xlim([0 5000]);

subplot(2,2,3);
bar(1:nM, [fm(:) ad(:)], 'grouped');
set(gca,'XTick',1:nM,'XTickLabel',labels,'XTickLabelRotation',20);
ylabel('%'); legend({'flow maintained','adequacy (\geq50% base)'},'Location','southeast');
title('Occlusion resilience'); box off; ylim([0 105]);

subplot(2,2,4); hold on;
for im=1:nM
    r=r_pial{im}; r=r(r<=300); [f,x]=ksdensity(r,'Bandwidth',max(2,range(r)/30));
    plot(x,f,'LineWidth',1.7,'Color',colors(im,:),'DisplayName',labels{im});
end
xline(r_dead*1e6,'k--','r_{dead}','HandleVisibility','off');
xlabel('pial radius [\mum]'); ylabel('density'); title('Radius distribution'); box off; xlim([0 200]);

sgtitle(sprintf('Perfusion & occlusion (density=%.2f /mm^2, n_{pen}=%d)', DENSITY, npen),'FontWeight','bold');
saveas(gcf, fullfile(PROJ,'FinalResults','Exp_Flow.png'));
save(fullfile(PROJ,'FinalResults','Exp_Flow_Results.mat'), ...
    'CVq','Qmean','Qdist','r_pial','flow_maint','adequacy','fm','ad','Qocc_dist', ...
    'occ_srcs','models','DENSITY','npen','Qreg','Qref');
fprintf('\nSaved: FinalResults/Exp_Flow.png and _Results.mat\n');

%% local helpers
function Qpen = solve_flows(S, r_full, Psrc, sink_idx)
    S.sourceP = Psrc;  for k=find(S.sources), S.SE(k).Ps=Psrc; end
    [S.IE.r] = vout(r_full);
    S = calcConductance2025(S); S = solvehemodyn2025(S);
    Qpen = abs([S.SE(sink_idx).Qs]);
end

function S = occlude(S, src)
    S.SE(src).Ps = S.sinkP;  S.SE(src).Gs = S.GsSink;  S.sources(src) = 0;
end
