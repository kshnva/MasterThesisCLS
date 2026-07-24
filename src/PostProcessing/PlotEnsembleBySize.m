% PlotEnsembleBySize.m
%
% Post-processing of RunGeometryEnsemble output: separate the loop-survival
% result by NETWORK SIZE (ncombsx), to expose the size trend that inflates
% FullCouplingFlowReg's variance in the pooled summary.
%
% No simulation — just loads FinalResults/GeometryEnsemble_Results.mat.
%
% Produces:
%   - a per-size table (mean loop fraction per model at each ncombsx)
%   - FinalResults/GeometryEnsemble_bySize.png  (loop fraction vs ncombsx, per model)

close all;
PROJ = fileparts(fileparts(mfilename('fullpath')));   % absolute so this works from any folder
load(fullfile(PROJ,'FinalResults','GeometryEnsemble_Results.mat'), ...
    'loops_frac','loops_abs','ncombsx_list','models','area_mm2','npen_used');

nM     = numel(models);
sizes  = unique(ncombsx_list(:)');           % distinct network sizes
colors = [0.15 0.35 0.75; 0.85 0.45 0.05; 0.55 0.10 0.55; 0.05 0.60 0.20];

% ── aggregate loop fraction by size ─────────────────────────────────────────
mean_by = nan(numel(sizes), nM);
sd_by   = nan(numel(sizes), nM);
n_by    = nan(numel(sizes),1);
area_by = nan(numel(sizes),1);
pen_by  = nan(numel(sizes),1);
for is = 1:numel(sizes)
    rows = (ncombsx_list(:)' == sizes(is));
    n_by(is)   = nnz(rows);
    area_by(is)= mean(area_mm2(rows));
    pen_by(is) = round(mean(npen_used(rows)));
    mean_by(is,:) = mean(loops_frac(rows,:), 1);
    sd_by(is,:)   = std(loops_frac(rows,:), 0, 1);
end

% ── table ───────────────────────────────────────────────────────────────────
fprintf('\n%s\nLOOP FRACTION BY NETWORK SIZE\n%s\n', repmat('=',1,72), repmat('=',1,72));
fprintf('%-8s %6s %6s %5s', 'ncombsx','area','n_pen','runs');
for im=1:nM, fprintf(' %16s', models{im}); end
fprintf('\n');
for is = 1:numel(sizes)
    fprintf('%-8d %5.0f  %5d %5d', sizes(is), area_by(is), pen_by(is), n_by(is));
    for im = 1:nM
        fprintf('   %5.2f+/-%4.2f', mean_by(is,im), sd_by(is,im));
    end
    fprintf('\n');
end
fprintf('\n(area in mm^2; loop fraction = n_loops / n_loops_total)\n');

% ── figure: loop fraction vs size, per model ────────────────────────────────
figure('Color','w','Position',[150 200 720 440]); hold on;
for im = 1:nM
    % faint individual realizations
    scatter(ncombsx_list(:) + 0.05*(im-2.5), loops_frac(:,im), 22, colors(im,:), ...
        'filled', 'MarkerFaceAlpha', 0.30, 'HandleVisibility','off');
    % mean +/- SD trend
    errorbar(sizes, mean_by(:,im), sd_by(:,im), '-o', 'LineWidth', 1.9, ...
        'Color', colors(im,:), 'MarkerFaceColor', colors(im,:), 'DisplayName', models{im});
end
xlabel('network size  (ncombsx)   \rightarrow larger area / more penetrators');
ylabel('surviving loop fraction');
title('Loop survival vs network size, by coupling model');
xticks(sizes); ylim([-0.02 1.05]); grid on; box off;
legend('Location','best','FontSize',9);
saveas(gcf, fullfile(PROJ,'FinalResults','GeometryEnsemble_bySize.png'));

% ── quick trend readout (does FlowReg rise with size, relative to others) ───
fprintf('\nSize trend (loop fraction, smallest -> largest ncombsx):\n');
for im = 1:nM
    d = mean_by(end,im) - mean_by(1,im);
    fprintf('  %-22s %.2f -> %.2f   (%+.2f)\n', models{im}, mean_by(1,im), mean_by(end,im), d);
end
fprintf('\nSaved: FinalResults/GeometryEnsemble_bySize.png\n');
