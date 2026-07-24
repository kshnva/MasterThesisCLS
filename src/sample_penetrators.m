function [sampled_penetrators] = sample_penetrators(penetrator_data, n_samples, sampling_method)
% Stochastic sampling of penetrator vessels from dataset
% penetrator_data: struct array with fields .r, .l, .R (radius, length, resistance)
% n_samples: number of penetrators to sample (e.g., 100)
% sampling_method: 'uniform', 'weighted', 'stratified'

if nargin < 3, sampling_method = 'uniform'; end

n_data = length(penetrator_data.r);  % Vector length, not struct length

switch sampling_method
    case 'uniform'
        % Simple random sampling
        sample_indices = randperm(n_data, min(n_samples, n_data));
        
    case 'weighted'
        % Weighted sampling based on some criterion (e.g., radius distribution)
        weights = penetrator_data.r / sum(penetrator_data.r);
        sample_indices = randsample(n_data, min(n_samples, n_data), true, weights);
        
    case 'stratified'
        % Stratified sampling by radius ranges
        radius_ranges = [0, 20e-6; 20e-6, 40e-6; 40e-6, 60e-6; 60e-6, inf];
        samples_per_stratum = ceil(n_samples / length(radius_ranges));
        
        sample_indices = [];
        for i = 1:length(radius_ranges)
            in_range = find(penetrator_data.r >= radius_ranges(i,1) & ...
                           penetrator_data.r < radius_ranges(i,2));
            if ~isempty(in_range)
                n_stratum = min(samples_per_stratum, length(in_range));
                stratum_samples = randsample(in_range, n_stratum);
                sample_indices = [sample_indices; stratum_samples];
            end
        end
        
        % Ensure we have exactly n_samples by adding more from available strata
        fprintf('  DEBUG PRE-WHILE: sample_indices=%d, n_samples=%d, n_data=%d\n', ...
                length(sample_indices), n_samples, n_data);
        while length(sample_indices) < n_samples && length(sample_indices) < n_data
            remaining = setdiff(1:n_data, sample_indices);
            if isempty(remaining), break; end
            n_needed = min(n_samples - length(sample_indices), length(remaining));
            extra = randsample(remaining, n_needed);
            extra = extra(:); % Ensure column vector
            sample_indices = [sample_indices; extra];
            fprintf('  DEBUG: Added %d extra samples, total now: %d\n', n_needed, length(sample_indices));
        end
        
    otherwise
        error('Unknown sampling method');
end

% Extract sampled data
sampled_penetrators.r = penetrator_data.r(sample_indices);
sampled_penetrators.l = penetrator_data.l(sample_indices);
sampled_penetrators.R = penetrator_data.R(sample_indices);
sampled_penetrators.G = 1 ./ penetrator_data.R(sample_indices); % conductance

end
