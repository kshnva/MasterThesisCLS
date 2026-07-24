function S = integrate_sampled_penetrators(S, sampled_penetrators)
% Integrate stochastically sampled penetrators into the leptomeningeal network
% S: existing simulation structure
% sampled_penetrators: struct with .r, .l, .G fields from sampling function

% Get current network statistics
n_pial_vessels = sum([S.IE.pial] == 1);
n_penetrators_needed = length(sampled_penetrators.r);

% Strategy: Distribute penetrators across pial vessels
% Option 1: One penetrator per pial vessel (if enough pial vessels)
% Option 2: Multiple penetrators per pial vessel (if more penetrators than pial vessels)

if n_penetrators_needed <= n_pial_vessels
    % Option 1: One penetrator per selected pial vessel
    pial_indices = find([S.IE.pial] == 1);
    selected_pial = randsample(pial_indices, n_penetrators_needed);
    
    for i = 1:n_penetrators_needed
        ie = selected_pial(i);
        S = add_sampled_penetrator(S, ie, sampled_penetrators.r(i), ...
                                   sampled_penetrators.l(i), sampled_penetrators.G(i));
    end
    
else
    % Option 2: Multiple penetrators per pial vessel
    pial_indices = find([S.IE.pial] == 1);
    penetrators_per_vessel = ceil(n_penetrators_needed / n_pial_vessels);
    
    penetrator_count = 1;
    for ie = pial_indices
        if penetrator_count > n_penetrators_needed
            break;
        end
        
        % Add multiple penetrators to this pial vessel
        n_this_vessel = min(penetrators_per_vessel, n_penetrators_needed - penetrator_count + 1);
        
        for j = 1:n_this_vessel
            if penetrator_count <= n_penetrators_needed
                S = add_sampled_penetrator(S, ie, sampled_penetrators.r(penetrator_count), ...
                                           sampled_penetrators.l(penetrator_count), ...
                                           sampled_penetrators.G(penetrator_count));
                penetrator_count = penetrator_count + 1;
            end
        end
    end
end

% Update node table and network connectivity
[S.IN, S.nIN] = MakeNodeTable(S.IE, S.SE, S.IN);

% Update element lengths based on positions
S.IE = LengthFromPosition(S.IE, S.IN);

% Update network counts
S.nIE = length(S.IE);
S.nSE = length(S.SE);
S.nIN = length(S.IN);

fprintf('Integrated %d penetrators into the network\n', n_penetrators_needed);
fprintf('Network now has: %d internal elements, %d source elements, %d nodes\n', ...
        S.nIE, S.nSE, S.nIN);

end

function S = add_sampled_penetrator(S, pial_ie, radius, length, conductance)
% Add a single sampled penetrator to a pial vessel
% S: simulation structure
% pial_ie: index of pial vessel to connect to
% radius, length, conductance: penetrator properties

% Get the pial vessel nodes
nodes = S.IE(pial_ie).nodes;
proximal_node = nodes(1);
distal_node = nodes(2);

% Create new node at penetrator origin (along the pial vessel)
new_node1 = S.nIN + 1;
S.IN(new_node1) = S.IN(proximal_node); % Copy properties
S.IN(new_node1).pos = S.IN(proximal_node).pos + ...
                      (S.IN(distal_node).pos - S.IN(proximal_node).pos) * 0.5;

% Create new node at penetrator distal end
new_node2 = S.nIN + 2;
S.IN(new_node2) = S.IN(new_node1);
S.IN(new_node2).pos = S.IN(new_node1).pos + [0; 0; -length]; % Dive into tissue

% Split the original pial vessel
S.IE(pial_ie).nodes = [proximal_node, new_node1];

% Create new pial segment (second half of original)
new_pial_ie = S.nIE + 1;
S.IE(new_pial_ie) = S.IE(pial_ie);
S.IE(new_pial_ie).nodes = [new_node1, distal_node];

% Create penetrator vessel with sampled properties
new_penetrator_ie = S.nIE + 2;
S.IE(new_penetrator_ie) = S.IE(pial_ie);
S.IE(new_penetrator_ie).nodes = [new_node1, new_node2];
S.IE(new_penetrator_ie).r = radius;
S.IE(new_penetrator_ie).l = length;
S.IE(new_penetrator_ie).G = conductance;
S.IE(new_penetrator_ie).pial = 0; % Mark as penetrator

% Create sink element at distal end
new_se = S.nSE + 1;
S.SE(new_se).node = new_node2;
S.SE(new_se).Ps = S.sinkP; % Venous pressure
S.SE(new_se).Gs = S.GsSink; % Sink conductance

% Update counts
S.nIN = S.nIN + 2;
S.nIE = S.nIE + 2;
S.nSE = S.nSE + 1;

end
