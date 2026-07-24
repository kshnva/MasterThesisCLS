function [nloops, nsurv, conv, r_ss] = sens_ss(S_ref, model_name, Qreg, ECr, Psrc, pial_idx, r_dead, maxchunks, r0_pial)
% Shared steady-state runner for the sensitivity-analysis scripts.
% Mirrors the proven tr_run logic used across the Phase-2 scripts.
%
% In:
%   S_ref     - built network struct (from sens_build)
%   model_name- 'NoCoupling'|'ArtCoupling'|'FullCoupling'|'FullCouplingFlowReg'
%   Qreg      - flow-regulation target (QrefCapillary)
%   ECr       - ECresistivity (sets lambda). Default 0.1 if [].
%   Psrc      - source pressure [Pa]
%   pial_idx  - indices of pial IE segments
%   r_dead    - collapse threshold [m]
%   maxchunks - chunk cap for convergence
%   r0_pial   - (optional) initial PIAL radii: scalar or vector [m]. [] = neutral 100 um.
% Out:
%   nloops, nsurv, conv (bool), r_ss (full radius vector)

    if nargin < 9, r0_pial = []; end
    if isempty(ECr), ECr = 0.1; end

    S = S_ref;
    S.sourceP = Psrc;
    for k = find(S.sources), S.SE(k).Ps = Psrc; end
    % canonical capillary parameters (same as run_ss in the Phase-2 scripts)
    S.gCapBase = 0.05;  S.kSensQ = 0.6;
    S.GscCapillary = 2e-5;  S.kSmCapillary = 7e-6;
    S.ECresistivity = ECr;  S.QrefCapillary = Qreg;

    pen_ie = find(~[S.IE.pial]);
    S.kReg = S.kReg * ones(S.nIE,1);  S.kReg(pen_ie) = 0;

    % initial radii: pial = r0_pial (or neutral 100 um), penetrators frozen at 15 um
    r0 = 100e-6*ones(1,S.nIE);  r0(pen_ie) = 15e-6;
    if ~isempty(r0_pial)
        r0(pial_idx) = r0_pial;             % scalar broadcast or same-length vector
    end
    S.r0 = r0;
    S.KeepDynamics = 0;

    nc = any(strcmp(model_name,{'FullCoupling','FullCouplingFlowReg'}));
    S.DoCapillaryCurrent = nc;
    S = MakeElectricCircuit(S,0);
    S.DoElectricCoupling       = ~strcmp(model_name,'NoCoupling');
    S.DoCapillaryCurrent       = nc;
    S.FlowRegulatedCapillaries = strcmp(model_name,'FullCouplingFlowReg');
    S.annot = model_name;

    opts = odeset('Events',@sens_nan,'RelTol',1e-6,'AbsTol',1e-9);
    t=0; y=S.r0(:); conv=false; nch=0;
    while nch < maxchunks
        [tc,Xc,tE] = ode45(@(t,X) rdotfunCoupling(t,X,S),[t,t+1e5],y,opts);
        nch = nch+1;
        if ~isempty(tE), y=Xc(end,:)'; break; end
        if max(abs(Xc(end,:)-Xc(1,:)))/(tc(end)-tc(1)) < 1e-12, conv=true; y=Xc(end,:)'; break; end
        t=tc(end); y=Xc(end,:)';
    end
    r_ss = y';
    alive = pial_idx(r_ss(pial_idx) > r_dead);
    nloops = sens_cycle_rank(S_ref, alive);
    nsurv  = numel(alive);
end

function [v,i,d] = sens_nan(~,X)
    v = all(~isnan(X)); i = 1; d = 0;
end
