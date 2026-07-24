function [cycle_rank, info] = sens_cycle_rank(S, seg_list)
% Cycle rank (number of independent loops) of a pial subgraph:
%   cycle_rank = m - n + c   (edges - nodes + connected components)
% Shared helper for the sensitivity-analysis scripts.
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
