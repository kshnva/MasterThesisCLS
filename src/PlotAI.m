% PlotAI.m
%
% Plots the autoregulation index from the saved sweep. No simulation here:
% loads AutoregulationCurve_Results.mat (produced by TestAutoregulationCurve.m).
%
% Two panels:
%   (A) AI per model: single global value (local 80-90 slope, Lassen form).
%   (B) local AI vs pressure: AI on every consecutive pressure bracket, showing
%       the regulation notch at the operating point and AI<0 away from it.

close all;
PROJ = fileparts(mfilename('fullpath'));
load(fullfile(PROJ,'FinalResults','AutoregulationCurve_Results.mat'), ...
    'Q_norm_all','pressures_mmHg','AI_all','models','labels','P_op');

nM = numel(models);
colors = [0.15 0.35 0.75;   % blue
          0.85 0.45 0.05;   % orange
          0.55 0.10 0.55;   % purple
          0.05 0.60 0.20];  % green

% Local AI at each consecutive bracket: AI = 1 - P_op * dQnorm/dP
Pmid   = (pressures_mmHg(1:end-1) + pressures_mmHg(2:end)) / 2;
AI_loc = nan(nM, numel(Pmid));
for im = 1:nM
    dQ = diff(Q_norm_all(im,:));
    dP = diff(pressures_mmHg);
    AI_loc(im,:) = 1 - P_op * (dQ ./ dP);
end

figure('Color','w','Position',[120 200 1000 400]);

% (A) global AI per model
subplot(1,2,1);
b = bar(1:nM, AI_all, 'FaceColor','flat','EdgeColor','none','BarWidth',0.6);
b.CData = colors(1:nM,:);
hold on; yline(0,'k-','LineWidth',1);
set(gca,'XTick',1:nM,'XTickLabel',labels,'XTickLabelRotation',20,'FontSize',10);
ylabel('Autoregulation Index  (local 80-90 mmHg)');
title('AI per model'); box off; grid on;
for im = 1:nM
    text(im, AI_all(im)+sign(AI_all(im))*0.02, sprintf('%+.3f',AI_all(im)), ...
        'HorizontalAlignment','center','FontSize',9);
end
text(0.55, 0.92, 'AI>0: regulating', 'Units','normalized','FontSize',8,'Color',[0.1 0.5 0.1]);
text(0.55, 0.86, 'AI<0: amplifying', 'Units','normalized','FontSize',8,'Color',[0.6 0.1 0.1]);

% (B) local AI vs pressure
subplot(1,2,2); hold on;
for im = 1:nM
    plot(Pmid, AI_loc(im,:), '-o','LineWidth',1.8,'Color',colors(im,:), ...
        'MarkerFaceColor',colors(im,:),'DisplayName',labels{im});
end
yline(0,'k-','LineWidth',1,'HandleVisibility','off');
xline(P_op,'k:','LineWidth',1,'Label','P_0=86','HandleVisibility','off');
xlabel('Arterial pressure  (mmHg)'); ylabel('local AI  (per bracket)');
title('Local AI vs pressure (regulation notch at set point)');
legend('Location','best','FontSize',8); box off; grid on;

sgtitle('Autoregulation Index','FontWeight','bold');
saveas(gcf,fullfile(PROJ,'FinalResults','AI_plot.png'));
fprintf('Saved AI_plot.png\n');

% also print the local-AI table
fprintf('\nLocal AI per bracket:\n');
fprintf('%-22s', 'bracket (mmHg)');
for k = 1:numel(Pmid)
    fprintf('%7d-%-3d', pressures_mmHg(k), pressures_mmHg(k+1));
end
fprintf('\n');
for im = 1:nM
    fprintf('%-22s', models{im});
    for k = 1:numel(Pmid), fprintf('%11.2f', AI_loc(im,k)); end
    fprintf('\n');
end
