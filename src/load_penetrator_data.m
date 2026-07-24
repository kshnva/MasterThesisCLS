function penetrator_data = load_penetrator_data(csv_file)
% Load penetrator data from CSV file and convert to simulation format
% csv_file: path to penetrators.csv file

fprintf('Loading penetrator data from %s...\n', csv_file);
% Read CSV file
try
    data = readtable(csv_file);
    disp(readtable('penetrators.csv').Properties.VariableNames)
catch
    % Fallback for older MATLAB versions
    data = csvread(csv_file, 1, 0); % Skip header
    disp(readtable('penetrators.csv').Properties.VariableNames)
    headers = {'label', 'zmin', 'zmax', 'z_span', 'z_span_um', 'radius_nm', ...
              'radius_um', 'radius_source', 'length_um', 'angle_deg', 'xy_diagonal'};
    data = cell2table(num2cell(data), 'VariableNames', headers);
end

% Extract relevant columns and convert units
n_penetrators = height(data);

% Convert radius from micrometers to meters
if ismember('radius_um', data.Properties.VariableNames)
    radii_um = data.radius_um;
    radii_m = radii_um * 1e-6; % Convert μm to m
else
    error('radius_um column not found in CSV');
end

% Convert length from micrometers to meters  
if ismember('length_um', data.Properties.VariableNames)
    lengths_um = data.length_um;
    lengths_m = lengths_um * 1e-6; % Convert μm to m
else
    error('length_um column not found in CSV');
end

% Calculate resistances using Poiseuille's law
% You can replace this with your actual resistance calculations if available
viscosity = 4e-3; % Ns/m² (blood viscosity)
resistances = 8 * viscosity * lengths_m ./ (pi * radii_m.^4);

% Create output structure
penetrator_data.r = radii_m;
penetrator_data.l = lengths_m;
penetrator_data.R = resistances;
penetrator_data.G = 1 ./ resistances; % Conductance

% Store additional metadata if needed
penetrator_data.labels = data.label;
penetrator_data.angles = data.angle_deg;
penetrator_data.z_span = data.z_span_um * 1e-6; % Convert to meters

% Display statistics
fprintf('Loaded %d penetrators from dataset\n', n_penetrators);
fprintf('Radius range: %.2f - %.2f μm\n', min(radii_um), max(radii_um));
fprintf('Length range: %.2f - %.2f μm\n', min(lengths_um), max(lengths_um));
fprintf('Resistance range: %.2e - %.2e Pa·s/m³\n', min(resistances), max(resistances));

end
